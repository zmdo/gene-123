#!/usr/bin/env python3
"""Draw a Manhattan plot for ChiVMV host-site association results."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ManhattanPoint:
    position: int
    minus_log10_p: float
    empirical_p: float
    q_value: float | None
    cramers_v: float
    reference_base: str
    host_specific: bool
    major_by_host: str


def parse_float(value: str) -> float | None:
    if not value or value.upper() == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_points(path: Path, max_q: float | None = None) -> list[ManhattanPoint]:
    points: list[ManhattanPoint] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            empirical_p = parse_float(row.get("empirical_p", ""))
            position_text = row.get("reference_position", "")
            if empirical_p is None or empirical_p <= 0 or not position_text.isdigit():
                continue
            q_value = parse_float(row.get("q_value", ""))
            if max_q is not None and (q_value is None or q_value > max_q):
                continue
            points.append(
                ManhattanPoint(
                    position=int(position_text),
                    minus_log10_p=-math.log10(empirical_p),
                    empirical_p=empirical_p,
                    q_value=q_value,
                    cramers_v=float(row.get("cramers_v", "0") or 0),
                    reference_base=row.get("reference_base", ""),
                    host_specific=row.get("host_specific", "").lower() == "true",
                    major_by_host=row.get("major_by_host", ""),
                )
            )
    if not points:
        raise ValueError(f"No rows with numeric empirical_p and reference_position found in {path}")
    points.sort(key=lambda point: point.position)
    return points


def draw_plot(points: list[ManhattanPoint], output: Path, title: str, annotate_top: int) -> None:
    import matplotlib.pyplot as plt

    positions = [point.position for point in points]
    values = [point.minus_log10_p for point in points]
    colors = ["#2563eb" if (point.position // 1000) % 2 == 0 else "#0f766e" for point in points]

    significant = [point for point in points if point.q_value is not None and point.q_value < 0.05]

    fig, ax = plt.subplots(figsize=(13, 6.5), dpi=160)
    ax.scatter(positions, values, c=colors, s=32, alpha=0.78, linewidths=0)

    if significant:
        ax.scatter(
            [point.position for point in significant],
            [point.minus_log10_p for point in significant],
            c="#dc2626",
            s=58,
            alpha=0.94,
            edgecolors="white",
            linewidths=0.5,
            label="q < 0.05",
            zorder=4,
        )

    ax.axhline(-math.log10(0.05), color="#64748b", linestyle="--", linewidth=1.1, label="p = 0.05")
    ax.set_title(title, fontsize=16, pad=14)
    ax.set_xlabel("Reference genome position (OP378160.1)", fontsize=12)
    ax.set_ylabel("-log10(empirical permutation p)", fontsize=12)
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(max(1, min(positions) - 150), max(positions) + 150)

    top_points = sorted(points, key=lambda point: (point.minus_log10_p, point.cramers_v), reverse=True)[:annotate_top]
    for rank, point in enumerate(top_points, start=1):
        offset = 12 if rank % 2 else -18
        ax.annotate(
            str(point.position),
            xy=(point.position, point.minus_log10_p),
            xytext=(0, offset),
            textcoords="offset points",
            ha="center",
            va="bottom" if offset > 0 else "top",
            fontsize=8,
            color="#111827",
            arrowprops={"arrowstyle": "-", "color": "#94a3b8", "lw": 0.6},
        )

    if significant:
        ax.legend(frameon=False, loc="upper right")
    else:
        ax.legend(frameon=False, loc="upper right")

    note = (
        f"Plotted {len(points)} permutation-tested candidate sites; "
        f"red points pass Benjamini-Hochberg q < 0.05."
    )
    fig.text(0.01, 0.01, note, fontsize=9, color="#475569")
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/results/host_site_associations.tsv"))
    parser.add_argument("--png", type=Path, default=Path("data/results/host_site_manhattan.png"))
    parser.add_argument("--svg", type=Path, default=Path("data/results/host_site_manhattan.svg"))
    parser.add_argument("--title", default="ChiVMV Host-Associated Sites")
    parser.add_argument("--annotate-top", type=int, default=12)
    parser.add_argument("--max-q", type=float, default=None, help="Optional q-value filter before plotting.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    points = read_points(args.input, max_q=args.max_q)
    draw_plot(points, args.png, args.title, args.annotate_top)
    draw_plot(points, args.svg, args.title, args.annotate_top)
    print(f"Plotted {len(points)} permutation-tested sites.")
    print(f"Wrote {args.png}")
    print(f"Wrote {args.svg}")


if __name__ == "__main__":
    main()
