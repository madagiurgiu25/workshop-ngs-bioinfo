#!/usr/bin/env python3
"""
Summarize and plot Picard CollectAlignmentSummaryMetrics output.

Usage:
    python plot_alignment_summary.py dedup/HG00171_female_ERR034564.alignment_summary.txt
    python plot_alignment_summary.py metrics.txt -o my_plot.png

Output:
    - a summary of key metrics printed to the terminal
    - a PNG with metric tiles and the read-length histogram

Requires: python3, matplotlib  (conda install -c conda-forge matplotlib)
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def to_number(value):
    """Convert a Picard value to float. Handles decimal commas (German locale)."""
    if value is None or value.strip() == "":
        return None
    try:
        return float(value.replace(",", "."))
    except ValueError:
        return None


def parse_picard(path):
    """Return (metrics_by_category, histogram) from a Picard metrics file."""
    lines = Path(path).read_text().splitlines()
    metrics, histogram = {}, []
    section, header = None, None

    for line in lines:
        if line.startswith("## METRICS CLASS"):
            section, header = "metrics", None
            continue
        if line.startswith("## HISTOGRAM"):
            section, header = "histogram", None
            continue
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.split("\t")
        if header is None:
            header = fields
            continue

        if section == "metrics":
            category = fields[0]
            # Numeric columns sit between CATEGORY and SAMPLE/LIBRARY/READ_GROUP
            numeric_header = header[1:header.index("SAMPLE")] if "SAMPLE" in header else header[1:]
            numeric_values = []
            for v in fields[1:]:
                if to_number(v) is None:
                    break
                numeric_values.append(v)
            if len(numeric_values) != len(numeric_header):
                print(f"WARNING: row '{category}' has {len(numeric_values)} numeric values, "
                      f"expected {len(numeric_header)}. Columns would be shifted, so this row "
                      f"is skipped and rebuilt from R1 + R2 where possible.", file=sys.stderr)
                continue
            metrics[category] = {h: to_number(v) for h, v in zip(numeric_header, numeric_values)}
        elif section == "histogram":
            values = [to_number(v) or 0 for v in fields]
            histogram.append(values)

    return metrics, histogram


def pair_from_reads(metrics):
    """Rebuild the PAIR row from FIRST_OF_PAIR and SECOND_OF_PAIR if needed."""
    if "PAIR" in metrics:
        return metrics["PAIR"]
    r1, r2 = metrics.get("FIRST_OF_PAIR"), metrics.get("SECOND_OF_PAIR")
    if not (r1 and r2):
        sys.exit("ERROR: no PAIR row and no FIRST/SECOND_OF_PAIR rows to rebuild it from.")
    total = r1["TOTAL_READS"] + r2["TOTAL_READS"]
    aligned = r1["PF_READS_ALIGNED"] + r2["PF_READS_ALIGNED"]
    weight = lambda key: (r1[key] * r1["PF_READS_ALIGNED"] + r2[key] * r2["PF_READS_ALIGNED"]) / aligned
    return {
        "TOTAL_READS": total,
        "PF_READS_ALIGNED": aligned,
        "PCT_PF_READS_ALIGNED": aligned / total,
        "PF_ALIGNED_BASES": r1["PF_ALIGNED_BASES"] + r2["PF_ALIGNED_BASES"],
        "PCT_READS_ALIGNED_IN_PAIRS": (r1["READS_ALIGNED_IN_PAIRS"] + r2["READS_ALIGNED_IN_PAIRS"]) / aligned,
        "PCT_PF_READS_IMPROPER_PAIRS": (r1["PF_READS_IMPROPER_PAIRS"] + r2["PF_READS_IMPROPER_PAIRS"]) / aligned,
        "PF_MISMATCH_RATE": weight("PF_MISMATCH_RATE"),
        "PF_INDEL_RATE": weight("PF_INDEL_RATE"),
        "STRAND_BALANCE": weight("STRAND_BALANCE"),
        "PCT_CHIMERAS": weight("PCT_CHIMERAS"),
        "PCT_ADAPTER": (r1["PCT_ADAPTER"] + r2["PCT_ADAPTER"]) / 2,
        "PCT_SOFTCLIP": weight("PCT_SOFTCLIP"),
        "MEAN_READ_LENGTH": (r1["MEAN_READ_LENGTH"] * r1["TOTAL_READS"]
                             + r2["MEAN_READ_LENGTH"] * r2["TOTAL_READS"]) / total,
    }


def summarize(metrics):
    pair = pair_from_reads(metrics)
    r1, r2 = metrics.get("FIRST_OF_PAIR", {}), metrics.get("SECOND_OF_PAIR", {})
    total, aligned = pair["TOTAL_READS"], pair["PF_READS_ALIGNED"]
    return [
        ("Total reads", f"{total:,.0f}",
         f"{r1.get('TOTAL_READS', 0):,.0f} R1 · {r2.get('TOTAL_READS', 0):,.0f} R2"),
        ("Reads aligned", f"{pair['PCT_PF_READS_ALIGNED']:.1%}",
         f"{aligned:,.0f} aligned · {total - aligned:,.0f} unaligned"),
        ("Aligned in pairs", f"{pair['PCT_READS_ALIGNED_IN_PAIRS']:.1%}",
         f"{pair['PCT_PF_READS_IMPROPER_PAIRS']:.1%} improper pairs"),
        ("Mismatch rate", f"{pair['PF_MISMATCH_RATE']:.2%}",
         f"R1 {r1.get('PF_MISMATCH_RATE', 0):.2%} · R2 {r2.get('PF_MISMATCH_RATE', 0):.2%}"),
        ("Indel rate", f"{pair['PF_INDEL_RATE']:.3%}", ""),
        ("Strand balance", f"{pair['STRAND_BALANCE']:.2f}", "0.50 = no strand bias"),
        ("Chimeras", f"{pair['PCT_CHIMERAS']:.2%}", ""),
        ("Adapter · soft-clip", f"{pair['PCT_ADAPTER']:.2%} · {pair['PCT_SOFTCLIP']:.1%}", ""),
        ("Mean read length", f"{pair['MEAN_READ_LENGTH']:.1f} bp",
         f"aligned reads only: {pair['PF_ALIGNED_BASES'] / aligned:.1f} bp aligned"),
    ]


def plot(summary, histogram, sample, out_png):
    fig = plt.figure(figsize=(11, 8))
    grid = fig.add_gridspec(3, 3, height_ratios=[1, 1, 2.4], hspace=0.45, wspace=0.15)

    for i, (label, value, note) in enumerate(summary[:6]):
        ax = fig.add_subplot(grid[i // 3, i % 3])
        ax.axis("off")
        ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, color="#f1efe8"))
        ax.text(0.05, 0.75, label, fontsize=10, color="#5f5e5a", transform=ax.transAxes)
        ax.text(0.05, 0.38, value, fontsize=18, color="#2c2c2a", transform=ax.transAxes)
        ax.text(0.05, 0.12, note, fontsize=8.5, color="#888780", transform=ax.transAxes)

    ax = fig.add_subplot(grid[2, :])
    if histogram:
        lengths = [row[0] for row in histogram]
        total = [row[1] for row in histogram]
        aligned = [row[2] for row in histogram]
        ax.bar(lengths, total, width=0.8, color="#2a78d6", label="Read length")
        ax.plot(lengths, aligned, "o-", color="#eb6834", markersize=3, linewidth=1.5,
                label="Aligned length")
        ax.set_yscale("log")
        ax.set_xlabel("Length (bp)")
        ax.set_ylabel("Reads (log scale)")
        ax.legend(frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        unaligned = next((row[2] for row in histogram if row[0] == 0), 0)
        if unaligned:
            ax.annotate(f"{unaligned:,.0f} unaligned reads", xy=(0, unaligned),
                        xytext=(12, unaligned), fontsize=9, color="#993c1d",
                        arrowprops=dict(arrowstyle="->", color="#993c1d"))
    else:
        ax.text(0.5, 0.5, "No histogram found", ha="center", transform=ax.transAxes)

    fig.suptitle(f"Picard alignment summary: {sample}", fontsize=14, x=0.125, ha="left")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to: {out_png}")


def main():
    parser = argparse.ArgumentParser(description="Summarize Picard alignment summary metrics.")
    parser.add_argument("metrics_file", help="Output of picard CollectAlignmentSummaryMetrics")
    parser.add_argument("-o", "--output", help="Output PNG (default: <input>.png)")
    args = parser.parse_args()

    sample = Path(args.metrics_file).name.replace(".alignment_summary.txt", "")
    out_png = args.output or str(Path(args.metrics_file).with_suffix(".png"))

    metrics, histogram = parse_picard(args.metrics_file)
    summary = summarize(metrics)

    print(f"Alignment summary: {sample}\n")
    for label, value, note in summary:
        print(f"  {label:<22}{value:<18}{note}")

    plot(summary, histogram, sample, out_png)


if __name__ == "__main__":
    main()
