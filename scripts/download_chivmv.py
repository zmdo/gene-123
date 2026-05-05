#!/usr/bin/env python3
"""Download complete ChiVMV genomes and host metadata from NCBI."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_QUERY = (
    '(("Chilli veinal mottle virus"[Organism]) OR ChiVMV[All Fields] '
    'OR chivmv[All Fields]) AND ("complete genome"[Title] '
    'OR "complete sequence"[Title])'
)
DEFAULT_MIN_LENGTH = 9000
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
METADATA_FIELDS = [
    "accession_version",
    "accession",
    "ncbi_uid",
    "organism",
    "definition",
    "length",
    "host",
    "host_normalized",
    "isolate",
    "strain",
    "clone",
    "country",
    "collection_date",
    "mol_type",
    "topology",
]


@dataclass
class GenomeRecord:
    accession_version: str
    accession: str
    ncbi_uid: str
    organism: str
    definition: str
    length: int
    host: str
    host_normalized: str
    isolate: str
    strain: str
    clone: str
    country: str
    collection_date: str
    mol_type: str
    topology: str
    sequence: str


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


def request_text(endpoint: str, params: dict[str, str | int], retries: int = 3) -> str:
    url = f"{EUTILS_BASE}/{endpoint}?{urlencode(params)}"
    request = Request(url, headers={"User-Agent": "chivmv-host-analysis/0.1"})
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except (HTTPError, URLError, TimeoutError) as error:
            last_error = error
            if attempt == retries:
                break
            time.sleep(2 * attempt)
    curl_path = shutil.which("curl")
    if curl_path:
        process = subprocess.run(
            [
                curl_path,
                "--location",
                "--fail",
                "--silent",
                "--show-error",
                "--retry",
                str(retries),
                "--retry-delay",
                "2",
                "--max-time",
                "60",
                url,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if process.returncode == 0:
            return process.stdout
        last_error = RuntimeError(process.stderr.strip() or "curl fallback failed")
    raise RuntimeError(f"NCBI request failed after {retries} attempts: {last_error}")


def ncbi_params(args: argparse.Namespace, extra: dict[str, str | int]) -> dict[str, str | int]:
    params: dict[str, str | int] = {"tool": "chivmv_host_analysis"}
    if args.email:
        params["email"] = args.email
    if args.api_key:
        params["api_key"] = args.api_key
    params.update(extra)
    return params


def esearch(args: argparse.Namespace) -> list[str]:
    params = ncbi_params(
        args,
        {
            "db": "nuccore",
            "term": args.term,
            "retmode": "json",
            "retmax": args.retmax,
            "sort": "relevance",
        },
    )
    payload = json.loads(request_text("esearch.fcgi", params))
    result = payload["esearchresult"]
    count = int(result["count"])
    ids = result.get("idlist", [])
    if count > len(ids):
        print(
            f"Warning: query matched {count} records but retmax={args.retmax}; only {len(ids)} will be fetched.",
            file=sys.stderr,
        )
    return ids


def batched(items: list[str], batch_size: int) -> Iterable[list[str]]:
    for start in range(0, len(items), batch_size):
        yield items[start : start + batch_size]


def source_qualifiers(gbseq: ET.Element) -> dict[str, str]:
    qualifiers: dict[str, str] = {}
    for feature in gbseq.findall("./GBSeq_feature-table/GBFeature"):
        if feature.findtext("GBFeature_key") != "source":
            continue
        for qualifier in feature.findall("./GBFeature_quals/GBQualifier"):
            name = qualifier.findtext("GBQualifier_name", default="")
            value = qualifier.findtext("GBQualifier_value", default="")
            if name and value and name not in qualifiers:
                qualifiers[name] = value
    return qualifiers


def extract_uid(gbseq: ET.Element) -> str:
    for seqid in gbseq.findall("./GBSeq_other-seqids/GBSeqid"):
        text = seqid.text or ""
        if text.startswith("gi|"):
            return text.split("|", 1)[1]
    return ""


def parse_record(gbseq: ET.Element) -> GenomeRecord:
    qualifiers = source_qualifiers(gbseq)
    accession_version = gbseq.findtext("GBSeq_accession-version", default="")
    accession = gbseq.findtext("GBSeq_primary-accession", default=accession_version.split(".")[0])
    sequence = (gbseq.findtext("GBSeq_sequence", default="") or "").upper()
    host = qualifiers.get("host", "")
    return GenomeRecord(
        accession_version=accession_version,
        accession=accession,
        ncbi_uid=extract_uid(gbseq),
        organism=gbseq.findtext("GBSeq_organism", default=""),
        definition=gbseq.findtext("GBSeq_definition", default=""),
        length=int(gbseq.findtext("GBSeq_length", default=str(len(sequence))) or len(sequence)),
        host=host,
        host_normalized=normalize_host(host),
        isolate=qualifiers.get("isolate", ""),
        strain=qualifiers.get("strain", ""),
        clone=qualifiers.get("clone", ""),
        country=qualifiers.get("country", ""),
        collection_date=qualifiers.get("collection_date", ""),
        mol_type=qualifiers.get("mol_type", ""),
        topology=gbseq.findtext("GBSeq_topology", default=""),
        sequence=sequence,
    )


def fetch_records(args: argparse.Namespace, ids: list[str]) -> tuple[list[GenomeRecord], ET.Element]:
    records: list[GenomeRecord] = []
    combined_xml = ET.Element("GBSet")
    sleep_seconds = args.sleep if args.sleep is not None else (0.12 if args.api_key else 0.35)

    for batch in batched(ids, args.batch_size):
        params = ncbi_params(
            args,
            {
                "db": "nuccore",
                "id": ",".join(batch),
                "rettype": "gb",
                "retmode": "xml",
            },
        )
        xml_text = request_text("efetch.fcgi", params)
        root = ET.fromstring(xml_text)
        for gbseq in root.findall("GBSeq"):
            combined_xml.append(gbseq)
            record = parse_record(gbseq)
            if record.length < args.min_length:
                continue
            if args.require_complete_title and not re.search(
                r"complete (genome|sequence)", record.definition, flags=re.IGNORECASE
            ):
                continue
            records.append(record)
        time.sleep(sleep_seconds)
    return records, combined_xml


def fasta_header(record: GenomeRecord) -> str:
    host = record.host_normalized.replace(" ", "_")
    isolate = re.sub(r"\s+", "_", record.isolate or record.clone or "unknown")
    return f"{record.accession_version} host={host} isolate={isolate} length={record.length}"


def write_fasta(records: list[GenomeRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{fasta_header(record)}\n")
            sequence = record.sequence
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def write_metadata(records: list[GenomeRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_FIELDS, delimiter="\t")
        writer.writeheader()
        for record in records:
            row = {field: getattr(record, field) for field in METADATA_FIELDS}
            writer.writerow(row)


def write_summary(records: list[GenomeRecord], args: argparse.Namespace, path: Path) -> None:
    host_counts = Counter(record.host_normalized for record in records)
    lengths = [record.length for record in records]
    summary = {
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "query": args.term,
        "record_count": len(records),
        "records_with_host": sum(1 for record in records if record.host_normalized != "unknown"),
        "host_counts": dict(sorted(host_counts.items(), key=lambda item: (-item[1], item[0]))),
        "min_length": min(lengths) if lengths else None,
        "max_length": max(lengths) if lengths else None,
    }
    path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def save_xml(root: ET.Element, path: Path) -> None:
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(path, encoding="utf-8", xml_declaration=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--term", default=DEFAULT_QUERY, help="NCBI Nucleotide search term.")
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL"), help="Email passed to NCBI E-utilities.")
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"), help="NCBI API key.")
    parser.add_argument("--retmax", type=int, default=1000, help="Maximum NCBI UIDs to fetch.")
    parser.add_argument("--batch-size", type=int, default=50, help="EFetch batch size.")
    parser.add_argument("--min-length", type=int, default=DEFAULT_MIN_LENGTH, help="Minimum genome length.")
    parser.add_argument(
        "--no-title-filter",
        action="store_false",
        dest="require_complete_title",
        help="Do not require 'complete genome' or 'complete sequence' in the GenBank definition.",
    )
    parser.add_argument("--sleep", type=float, default=None, help="Seconds to sleep between EFetch batches.")
    parser.add_argument("--outdir", type=Path, default=Path("data"), help="Output directory root.")
    parser.set_defaults(require_complete_title=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_dir = args.outdir / "raw"
    metadata_dir = args.outdir / "metadata"
    results_dir = args.outdir / "results"
    for directory in (raw_dir, metadata_dir, results_dir):
        directory.mkdir(parents=True, exist_ok=True)

    ids = esearch(args)
    if not ids:
        raise SystemExit("No NCBI Nucleotide records matched the query.")

    records, xml_root = fetch_records(args, ids)
    if not records:
        raise SystemExit("No records passed the complete-genome filters.")

    fasta_path = raw_dir / "chivmv_genomes.fasta"
    metadata_path = metadata_dir / "chivmv_metadata.tsv"
    summary_path = results_dir / "download_summary.json"
    xml_path = raw_dir / "chivmv_records.xml"

    write_fasta(records, fasta_path)
    write_metadata(records, metadata_path)
    write_summary(records, args, summary_path)
    save_xml(xml_root, xml_path)

    host_counts = Counter(record.host_normalized for record in records)
    print(f"Downloaded {len(records)} ChiVMV complete-genome records.")
    print(f"Wrote {fasta_path}")
    print(f"Wrote {metadata_path}")
    print(f"Wrote {summary_path}")
    print("Top hosts:")
    for host, count in host_counts.most_common(10):
        print(f"  {host}: {count}")


if __name__ == "__main__":
    main()
