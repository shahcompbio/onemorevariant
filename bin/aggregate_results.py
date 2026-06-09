#!/usr/bin/env python3
"""Aggregate stratified variant calls into benchmark summary statistics.

Reads one or more stratified_calls.tsv files and computes precision, recall,
and F1 per sample × caller × category, plus genome-wide totals.
"""

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", nargs="+", required=True, help="Input stratified_calls.tsv files"
    )
    parser.add_argument(
        "--variant-type", required=True, help="Variant type (snv/indel)"
    )
    parser.add_argument(
        "--output-summary", required=True, help="Output benchmark summary TSV"
    )
    parser.add_argument(
        "--output-calls", required=True, help="Output concatenated stratified calls TSV"
    )
    args = parser.parse_args()

    # Read and concatenate all input files
    df = pd.concat([pd.read_csv(f, sep="\t") for f in args.input], ignore_index=True)

    # Write concatenated calls
    df.to_csv(args.output_calls, sep="\t", index=False)

    # Compute per-category summary: group by sample × caller × category
    per_category = (
        df.groupby(["sample", "caller", "category", "classification"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["TP", "FP", "FN"], fill_value=0)
        .reset_index()
        .groupby(["sample", "caller", "category"])[["TP", "FP", "FN"]]
        .sum()
        .reset_index()
    )

    # Compute genome-wide summary: deduplicate variants first (a variant in
    # multiple categories should only be counted once genome-wide)
    gw = (
        df.drop_duplicates(
            subset=["sample", "caller", "chrom", "pos", "ref", "alt", "classification"]
        )
        .groupby(["sample", "caller", "classification"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["TP", "FP", "FN"], fill_value=0)
        .reset_index()
    )
    gw["category"] = "genome_wide"

    # Combine
    summary = pd.concat([per_category, gw], ignore_index=True)
    summary["variant_type"] = args.variant_type
    summary["precision"] = summary["TP"] / (summary["TP"] + summary["FP"]).replace(
        0, float("nan")
    )
    summary["recall"] = summary["TP"] / (summary["TP"] + summary["FN"]).replace(
        0, float("nan")
    )
    summary["f1"] = (
        2
        * summary["precision"]
        * summary["recall"]
        / (summary["precision"] + summary["recall"]).replace(0, float("nan"))
    )
    summary[["precision", "recall", "f1"]] = (
        summary[["precision", "recall", "f1"]].fillna(0.0).round(4)
    )

    # Write summary
    col_order = [
        "sample",
        "caller",
        "variant_type",
        "category",
        "TP",
        "FP",
        "FN",
        "precision",
        "recall",
        "f1",
    ]
    summary[col_order].sort_values(["sample", "caller", "category"]).to_csv(
        args.output_summary, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
