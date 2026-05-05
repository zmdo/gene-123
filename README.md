# ChiVMV Genome and Host Association Analysis

This project downloads complete Chilli veinal mottle virus (ChiVMV) genome
records from NCBI Nucleotide, extracts host metadata from GenBank source
features, and ranks genome sites whose allele distributions differ by host.

## Outputs

- `data/raw/chivmv_genomes.fasta`: downloaded complete genome sequences.
- `data/metadata/chivmv_metadata.tsv`: accession, length, host, isolate, country,
  collection date, and GenBank definition.
- `data/results/download_summary.json`: download counts and host totals.
- `data/results/host_site_associations.tsv`: ranked host-associated candidate
  sites.
- `data/results/host_site_summary.md`: short human-readable analysis summary.
- `data/results/host_site_manhattan.png`: Manhattan plot of permutation-tested
  candidate sites.

## Run

Create or activate a Python environment, then install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Download genomes and metadata:

```bash
python scripts/download_chivmv.py --email your.email@example.com
```

Rank candidate host-associated sites:

```bash
python scripts/analyze_host_sites.py --permutations 1000
```

Draw a Manhattan plot:

```bash
python scripts/plot_manhattan.py
```

For a formal whole-genome multiple sequence alignment, install MAFFT first:

```bash
brew install mafft
python scripts/analyze_host_sites.py --mode mafft --permutations 1000
```

Without MAFFT, the analysis falls back to pairwise reference-coordinate mapping
with Biopython. That fallback is useful for substitutions and deletions relative
to a selected reference genome, but it does not evaluate insertion-only columns.

## Notes

- NCBI recommends providing an email address and an API key for heavier use. The
  scripts accept `--email` and `--api-key`, or `NCBI_EMAIL` and `NCBI_API_KEY`.
  If Python SSL fails against NCBI but `curl` works, the downloader falls back
  to `curl` automatically.
- Host names are normalized lightly so values like `Capsicum annuum L.`,
  `Capsicum annuum cultivar ...`, `hot pepper`, `pepper`, and `chilli` group
  under `Capsicum annuum`.
- The association statistic is exploratory. Candidate sites should be checked
  against sampling bias, geography, isolate history, and recombination before any
  biological interpretation.
