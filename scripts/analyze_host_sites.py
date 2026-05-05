#!/usr/bin/env python3
"""Rank ChiVMV genome sites whose allele distributions differ by host."""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
import shutil
import subprocess
import sys
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path


VALID_BASES = {"A", "C", "G", "T"}
MISSING_BASES = {"", "N", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "?"}
COMMON_HOST_MAP = {
    "chile pepper": "Capsicum annuum",
    "chili": "Capsicum annuum",
    "chili pepper": "Capsicum annuum",
    "chilli": "Capsicum annuum",
    "chilli pepper": "Capsicum annuum",
    "chilli prepper": "Capsicum annuum",
    "hot pepper": "Capsicum annuum",
    "pepper": "Capsicum annuum",
    "sweet pepper": "Capsicum annuum",
    "capsicum": "Capsicum spp.",
    "tobacco": "Nicotiana tabacum",
    "tomato": "Solanum lycopersicum",
}


@dataclass
class SiteResult:
    mode: str
    reference_id: str
    alignment_column: str
    reference_position: str
    reference_base: str
    n_sequences: int
    n_hosts: int
    alleles: str
    allele_counts: str
    major_by_host: str
    chi_square: float
    cramers_v: float
    empirical_p: str = "NA"
    q_value: str = "NA"
    host_specific: str = "false"


def normalize_host(host: str | None) -> str:
    if not host:
        return "unknown"
    normalized = re.sub(r"\s+", " ", host.strip())
    normalized = normalized.replace("'", "")
    if not normalized or normalized.lower() in {"na", "n/a", "unknown", "not provided"}:
        return "unknown"
    common_name = normalized.lower()
    if common_name in COMMON_HOST_MAP:
        return COMMON_HOST_MAP[common_name]
    normalized = re.sub(r"\b(cultivar|cv\.?|var\.?)\b.*$", "", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\b[Ll]\.$", "", normalized).strip()
    parts = normalized.split()
    if len(parts) >= 2 and re.fullmatch(r"[A-Z][A-Za-z-]+", parts[0]) and re.fullmatch(
        r"[a-z][A-Za-z-]+", parts[1]
    ):
        return f"{parts[0]} {parts[1]}"
    return normalized.lower()


def read_fasta(path: Path) -> OrderedDict[str, str]:
    sequences: OrderedDict[str, str] = OrderedDict()
    current_id: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(parts).upper()
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(parts).upper()
    if not sequences:
        raise ValueError(f"No FASTA sequences found in {path}")
    return sequences


def write_fasta(sequences: OrderedDict[str, str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for sequence_id, sequence in sequences.items():
            handle.write(f">{sequence_id}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            host = row.get("host_normalized") or normalize_host(row.get("host"))
            row["host_normalized"] = normalize_host(host)
            for key in (row.get("accession_version"), row.get("accession")):
                if key:
                    by_id[key] = row
    return by_id


def accession_stem(sequence_id: str) -> str:
    return sequence_id.split(".", 1)[0]


def metadata_for_sequence(sequence_id: str, metadata: dict[str, dict[str, str]]) -> dict[str, str] | None:
    return metadata.get(sequence_id) or metadata.get(accession_stem(sequence_id))


def resolve_reference(reference: str | None, sequences: OrderedDict[str, str]) -> str:
    if reference is None:
        return max(sequences, key=lambda sequence_id: len(sequences[sequence_id].replace("-", "")))
    if reference in sequences:
        return reference
    matches = [sequence_id for sequence_id in sequences if accession_stem(sequence_id) == accession_stem(reference)]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(f"Reference {reference!r} was not found in the FASTA IDs.")


def run_mafft(input_fasta: Path, output_alignment: Path, mafft_cmd: str) -> None:
    mafft_path = shutil.which(mafft_cmd)
    if not mafft_path:
        raise FileNotFoundError(f"MAFFT executable not found: {mafft_cmd}")
    with output_alignment.open("w", encoding="utf-8") as stdout_handle:
        process = subprocess.run(
            [mafft_path, "--auto", str(input_fasta)],
            stdout=stdout_handle,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if process.returncode != 0:
        raise RuntimeError(process.stderr.strip() or "MAFFT failed")


def usable_base(base: str, include_gaps: bool) -> bool:
    if base in VALID_BASES:
        return True
    if include_gaps and base == "-":
        return True
    return False


def host_sequence_ids(
    sequences: OrderedDict[str, str], metadata: dict[str, dict[str, str]], min_per_host: int
) -> tuple[list[str], dict[str, str], Counter[str]]:
    sequence_hosts: dict[str, str] = {}
    counts: Counter[str] = Counter()
    for sequence_id in sequences:
        row = metadata_for_sequence(sequence_id, metadata)
        host = normalize_host(row.get("host_normalized") if row else "")
        if host == "unknown":
            continue
        sequence_hosts[sequence_id] = host
        counts[host] += 1
    allowed_hosts = {host for host, count in counts.items() if count >= min_per_host}
    kept_ids = [sequence_id for sequence_id in sequences if sequence_hosts.get(sequence_id) in allowed_hosts]
    kept_counts = Counter(sequence_hosts[sequence_id] for sequence_id in kept_ids)
    if len(kept_counts) < 2:
        raise ValueError(
            "Need at least two host groups after filtering. Lower --min-per-host or add records with host metadata."
        )
    return kept_ids, sequence_hosts, kept_counts


def contingency_chi_square(table: dict[str, Counter[str]]) -> float:
    hosts = list(table)
    alleles = sorted({allele for counts in table.values() for allele in counts})
    row_totals = {host: sum(table[host].values()) for host in hosts}
    col_totals = {allele: sum(table[host][allele] for host in hosts) for allele in alleles}
    total = sum(row_totals.values())
    if total == 0:
        return 0.0
    chi_square = 0.0
    for host in hosts:
        for allele in alleles:
            expected = row_totals[host] * col_totals[allele] / total
            if expected > 0:
                observed = table[host][allele]
                chi_square += (observed - expected) ** 2 / expected
    return chi_square


def cramers_v(chi_square: float, total: int, n_hosts: int, n_alleles: int) -> float:
    denominator = total * min(n_hosts - 1, n_alleles - 1)
    if denominator <= 0:
        return 0.0
    return math.sqrt(chi_square / denominator)


def format_counter(counter: Counter[str]) -> str:
    return ";".join(f"{key}:{counter[key]}" for key in sorted(counter))


def site_result_from_calls(
    mode: str,
    reference_id: str,
    alignment_column: str,
    reference_position: str,
    reference_base: str,
    calls: dict[str, str],
    kept_ids: list[str],
    sequence_hosts: dict[str, str],
    include_gaps: bool,
    min_samples: int,
    min_major_frequency: float,
) -> SiteResult | None:
    table: dict[str, Counter[str]] = defaultdict(Counter)
    for sequence_id in kept_ids:
        base = calls.get(sequence_id, "N").upper()
        if not usable_base(base, include_gaps) or base in MISSING_BASES:
            continue
        table[sequence_hosts[sequence_id]][base] += 1

    table = {host: counts for host, counts in table.items() if sum(counts.values()) > 0}
    total = sum(sum(counts.values()) for counts in table.values())
    alleles = sorted({allele for counts in table.values() for allele in counts})
    if total < min_samples or len(table) < 2 or len(alleles) < 2:
        return None

    chi_square = contingency_chi_square(table)
    effect = cramers_v(chi_square, total, len(table), len(alleles))
    overall_counts: Counter[str] = Counter()
    major_parts: list[str] = []
    major_alleles: set[str] = set()
    strong_major_count = 0
    for host in sorted(table):
        counts = table[host]
        overall_counts.update(counts)
        host_total = sum(counts.values())
        major_allele, major_count = counts.most_common(1)[0]
        major_frequency = major_count / host_total
        major_alleles.add(major_allele)
        if major_frequency >= min_major_frequency:
            strong_major_count += 1
        major_parts.append(f"{host}:{major_allele}({major_count}/{host_total},{major_frequency:.3f})")

    host_specific = len(major_alleles) > 1 and strong_major_count >= 2
    return SiteResult(
        mode=mode,
        reference_id=reference_id,
        alignment_column=alignment_column,
        reference_position=reference_position,
        reference_base=reference_base,
        n_sequences=total,
        n_hosts=len(table),
        alleles=",".join(alleles),
        allele_counts=format_counter(overall_counts),
        major_by_host=";".join(major_parts),
        chi_square=chi_square,
        cramers_v=effect,
        host_specific=str(host_specific).lower(),
    )


def iter_msa_sites(
    sequences: OrderedDict[str, str], reference_id: str, skip_reference_insertions: bool = True
) -> tuple[str, list[tuple[str, str, str, str, dict[str, str]]]]:
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError("Aligned FASTA sequences must all have the same length.")
    alignment_length = lengths.pop()
    reference_sequence = sequences[reference_id]
    reference_position = 0
    sites: list[tuple[str, str, str, str, dict[str, str]]] = []
    for column_index in range(alignment_length):
        reference_base = reference_sequence[column_index]
        if reference_base != "-":
            reference_position += 1
        elif skip_reference_insertions:
            continue
        calls = {sequence_id: sequence[column_index] for sequence_id, sequence in sequences.items()}
        position_label = str(reference_position) if reference_base != "-" else "NA"
        sites.append((str(column_index + 1), position_label, reference_base, reference_base, calls))
    return "mafft_msa", sites


def reference_coordinate_calls(
    sequences: OrderedDict[str, str], reference_id: str
) -> tuple[str, list[tuple[str, str, str, str, dict[str, str]]]]:
    try:
        from Bio import Align
    except ImportError as error:
        raise RuntimeError(
            "Biopython is required for reference-coordinate fallback. Install requirements.txt or use MAFFT."
        ) from error

    reference_sequence = sequences[reference_id].replace("-", "").upper()
    mapped: dict[str, list[str]] = {
        sequence_id: ["N"] * len(reference_sequence) for sequence_id in sequences
    }
    mapped[reference_id] = list(reference_sequence)

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -1

    for sequence_id, sequence in sequences.items():
        if sequence_id == reference_id:
            continue
        query_sequence = sequence.replace("-", "").upper()
        alignment = aligner.align(reference_sequence, query_sequence)[0]
        coordinates = alignment.coordinates
        calls = ["N"] * len(reference_sequence)
        for segment_index in range(coordinates.shape[1] - 1):
            ref_start = int(coordinates[0, segment_index])
            ref_end = int(coordinates[0, segment_index + 1])
            query_start = int(coordinates[1, segment_index])
            query_end = int(coordinates[1, segment_index + 1])

            ref_advance = ref_end - ref_start
            query_advance = query_end - query_start
            if ref_advance > 0 and query_advance > 0:
                for offset in range(ref_advance):
                    calls[ref_start + offset] = query_sequence[query_start + offset]
            elif ref_advance > 0 and query_advance == 0:
                for offset in range(ref_advance):
                    calls[ref_start + offset] = "-"
        mapped[sequence_id] = calls

    sites: list[tuple[str, str, str, str, dict[str, str]]] = []
    for position_index, reference_base in enumerate(reference_sequence, start=1):
        calls = {sequence_id: mapped[sequence_id][position_index - 1] for sequence_id in sequences}
        sites.append(("NA", str(position_index), reference_base, reference_base, calls))
    return "reference_pairwise", sites


def build_site_results(
    mode: str,
    sites: list[tuple[str, str, str, str, dict[str, str]]],
    reference_id: str,
    kept_ids: list[str],
    sequence_hosts: dict[str, str],
    args: argparse.Namespace,
) -> list[SiteResult]:
    results: list[SiteResult] = []
    for alignment_column, reference_position, reference_base, _, calls in sites:
        result = site_result_from_calls(
            mode=mode,
            reference_id=reference_id,
            alignment_column=alignment_column,
            reference_position=reference_position,
            reference_base=reference_base,
            calls=calls,
            kept_ids=kept_ids,
            sequence_hosts=sequence_hosts,
            include_gaps=args.include_gaps,
            min_samples=args.min_samples,
            min_major_frequency=args.min_major_frequency,
        )
        if result is not None:
            results.append(result)
    results.sort(key=lambda item: (item.cramers_v, item.chi_square, item.n_sequences), reverse=True)
    return results


def empirical_p_value(
    result: SiteResult,
    calls: dict[str, str],
    kept_ids: list[str],
    sequence_hosts: dict[str, str],
    include_gaps: bool,
    permutations: int,
    rng: random.Random,
) -> float:
    alleles: list[str] = []
    hosts: list[str] = []
    for sequence_id in kept_ids:
        base = calls.get(sequence_id, "N").upper()
        if usable_base(base, include_gaps) and base not in MISSING_BASES:
            alleles.append(base)
            hosts.append(sequence_hosts[sequence_id])

    observed = result.chi_square
    exceedances = 0
    for _ in range(permutations):
        shuffled_hosts = hosts[:]
        rng.shuffle(shuffled_hosts)
        table: dict[str, Counter[str]] = defaultdict(Counter)
        for host, allele in zip(shuffled_hosts, alleles, strict=True):
            table[host][allele] += 1
        if contingency_chi_square(table) >= observed:
            exceedances += 1
    return (exceedances + 1) / (permutations + 1)


def add_empirical_p_values(
    results: list[SiteResult],
    sites: list[tuple[str, str, str, str, dict[str, str]]],
    kept_ids: list[str],
    sequence_hosts: dict[str, str],
    args: argparse.Namespace,
) -> None:
    if args.permutations <= 0 or not results:
        return
    by_site = {(site[0], site[1]): site[4] for site in sites}
    rng = random.Random(args.seed)
    tested = results[: args.permutation_top_sites]
    for result in tested:
        calls = by_site[(result.alignment_column, result.reference_position)]
        p_value = empirical_p_value(
            result=result,
            calls=calls,
            kept_ids=kept_ids,
            sequence_hosts=sequence_hosts,
            include_gaps=args.include_gaps,
            permutations=args.permutations,
            rng=rng,
        )
        result.empirical_p = f"{p_value:.6g}"

    p_values = [(index, float(result.empirical_p)) for index, result in enumerate(tested)]
    p_values.sort(key=lambda item: item[1])
    m_tests = len(p_values)
    adjusted: dict[int, float] = {}
    running_min = 1.0
    for rank_from_end, (index, p_value) in enumerate(reversed(p_values), start=1):
        rank = m_tests - rank_from_end + 1
        running_min = min(running_min, p_value * m_tests / rank)
        adjusted[index] = running_min
    for index, q_value in adjusted.items():
        tested[index].q_value = f"{min(q_value, 1.0):.6g}"


def write_results(results: list[SiteResult], path: Path) -> None:
    fields = [
        "mode",
        "reference_id",
        "alignment_column",
        "reference_position",
        "reference_base",
        "n_sequences",
        "n_hosts",
        "alleles",
        "allele_counts",
        "major_by_host",
        "chi_square",
        "cramers_v",
        "empirical_p",
        "q_value",
        "host_specific",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for result in results:
            row = result.__dict__.copy()
            row["chi_square"] = f"{result.chi_square:.6f}"
            row["cramers_v"] = f"{result.cramers_v:.6f}"
            writer.writerow(row)


def write_summary(
    results: list[SiteResult],
    host_counts: Counter[str],
    reference_id: str,
    mode: str,
    args: argparse.Namespace,
    path: Path,
) -> None:
    lines = [
        "# ChiVMV Host-Associated Site Summary",
        "",
        f"Mode: `{mode}`",
        f"Reference: `{reference_id}`",
        f"Candidate variable sites ranked: {len(results)}",
        f"Min records per host: {args.min_per_host}",
        f"Permutations per tested top site: {args.permutations}",
        "",
        "## Host Groups",
        "",
    ]
    for host, count in host_counts.most_common():
        lines.append(f"- {host}: {count}")
    lines.extend(["", "## Top Candidate Sites", ""])
    if not results:
        lines.append("No candidate sites passed the filters.")
    else:
        lines.append(
            "| Rank | Reference position | Ref base | Cramer's V | Empirical p | Major alleles by host |"
        )
        lines.append("| ---: | ---: | :---: | ---: | ---: | --- |")
        for rank, result in enumerate(results[:20], start=1):
            lines.append(
                f"| {rank} | {result.reference_position} | {result.reference_base} | "
                f"{result.cramers_v:.3f} | {result.empirical_p} | {result.major_by_host} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare_sites(
    sequences: OrderedDict[str, str], reference_id: str, args: argparse.Namespace
) -> tuple[str, OrderedDict[str, str], list[tuple[str, str, str, str, dict[str, str]]]]:
    if args.alignment:
        alignment_sequences = read_fasta(args.alignment)
        alignment_reference = resolve_reference(reference_id, alignment_sequences)
        mode, sites = iter_msa_sites(alignment_sequences, alignment_reference)
        return mode, alignment_sequences, sites

    if args.mode in {"auto", "mafft"} and shutil.which(args.mafft_cmd):
        args.outdir.mkdir(parents=True, exist_ok=True)
        alignment_path = args.outdir / "chivmv_aligned.fasta"
        run_mafft(args.fasta, alignment_path, args.mafft_cmd)
        alignment_sequences = read_fasta(alignment_path)
        alignment_reference = resolve_reference(reference_id, alignment_sequences)
        mode, sites = iter_msa_sites(alignment_sequences, alignment_reference)
        return mode, alignment_sequences, sites

    if args.mode == "mafft":
        raise SystemExit("MAFFT was requested but not found. Install it with `brew install mafft` or provide --alignment.")

    mode, sites = reference_coordinate_calls(sequences, reference_id)
    return mode, sequences, sites


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=Path("data/raw/chivmv_genomes.fasta"))
    parser.add_argument("--metadata", type=Path, default=Path("data/metadata/chivmv_metadata.tsv"))
    parser.add_argument("--alignment", type=Path, default=None, help="Existing aligned FASTA to analyze.")
    parser.add_argument("--outdir", type=Path, default=Path("data/results"))
    parser.add_argument("--reference", default=None, help="Reference accession/version. Default: longest sequence.")
    parser.add_argument("--mode", choices=["auto", "mafft", "reference"], default="auto")
    parser.add_argument("--mafft-cmd", default="mafft")
    parser.add_argument("--min-per-host", type=int, default=2)
    parser.add_argument("--min-samples", type=int, default=4)
    parser.add_argument("--include-gaps", action="store_true")
    parser.add_argument("--min-major-frequency", type=float, default=0.75)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--permutation-top-sites", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260505)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.fasta.exists():
        raise SystemExit(f"Missing FASTA: {args.fasta}. Run scripts/download_chivmv.py first.")
    if not args.metadata.exists():
        raise SystemExit(f"Missing metadata: {args.metadata}. Run scripts/download_chivmv.py first.")

    sequences = read_fasta(args.fasta)
    metadata = read_metadata(args.metadata)
    reference_id = resolve_reference(args.reference, sequences)
    mode, analysis_sequences, sites = prepare_sites(sequences, reference_id, args)
    analysis_reference = resolve_reference(reference_id, analysis_sequences)
    kept_ids, sequence_hosts, host_counts = host_sequence_ids(
        analysis_sequences, metadata, min_per_host=args.min_per_host
    )

    results = build_site_results(
        mode=mode,
        sites=sites,
        reference_id=analysis_reference,
        kept_ids=kept_ids,
        sequence_hosts=sequence_hosts,
        args=args,
    )
    add_empirical_p_values(results, sites, kept_ids, sequence_hosts, args)

    args.outdir.mkdir(parents=True, exist_ok=True)
    result_path = args.outdir / "host_site_associations.tsv"
    summary_path = args.outdir / "host_site_summary.md"
    write_results(results, result_path)
    write_summary(results, host_counts, analysis_reference, mode, args, summary_path)

    print(f"Analysis mode: {mode}")
    print(f"Reference: {analysis_reference}")
    print(f"Host groups retained: {dict(host_counts)}")
    print(f"Candidate sites ranked: {len(results)}")
    print(f"Wrote {result_path}")
    print(f"Wrote {summary_path}")
    if mode == "reference_pairwise":
        print(
            "Note: MAFFT was not found, so insertion-only alignment columns were not analyzed.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()