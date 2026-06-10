#!/usr/bin/env python3
"""Stratify TP/FP/FN variants by genomic region using bedtools intersect.

Produces a long-format TSV with one row per variant × category overlap,
annotated with VCF quality metrics from both the query and truth VCFs.
"""

import argparse
import csv
import gzip
import os
import subprocess
import tempfile


def parse_vcf_with_metrics(vcf_path):
    """Parse VCF and return dict keyed by (chrom, pos, ref, alt) with quality metrics.

    Extracts:
      - qual: QUAL column
      - gq: FORMAT/GQ
      - dp: FORMAT/DP
      - af: FORMAT/AF
      - naf: FORMAT/NAF
      - ndp: FORMAT/NDP
      - haplotype_support: INFO/H flag (1 if present, 0 otherwise)
    """
    records = {}
    format_fields_of_interest = {"GQ", "DP", "AF", "NAF", "NDP"}

    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5] if fields[5] != "." else ""

            # Parse INFO for H flag
            info = fields[7] if len(fields) > 7 else ""
            info_fields = info.split(";") if info else []
            haplotype_support = "1" if "H" in info_fields else "0"

            # Parse FORMAT fields
            metrics = {f: "" for f in format_fields_of_interest}
            if len(fields) > 9:
                fmt_keys = fields[8].split(":")
                fmt_vals = fields[9].split(":")
                for key, val in zip(fmt_keys, fmt_vals):
                    if key in format_fields_of_interest:
                        metrics[key] = val if val != "." else ""

            records[(chrom, pos, ref, alt)] = {
                "qual": qual,
                "gq": metrics.get("GQ", ""),
                "dp": metrics.get("DP", ""),
                "af": metrics.get("AF", ""),
                "naf": metrics.get("NAF", ""),
                "ndp": metrics.get("NDP", ""),
                "haplotype_support": haplotype_support,
            }

    return records


def vcf_to_bed(vcf_path, bed_path):
    """Convert VCF to BED format for bedtools intersect."""
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh, open(bed_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            # BED is 0-based half-open
            start = pos - 1
            end = pos - 1 + len(ref)
            out.write(f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\n")


def intersect_bed(variant_bed, region_bed):
    """Return set of (chrom, start, end, ref, alt) tuples that overlap the region BED."""
    result = set()
    try:
        proc = subprocess.run(
            ["bedtools", "intersect", "-u", "-a", variant_bed, "-b", region_bed],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError:
        return result

    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        result.add((fields[0], fields[1], fields[2], fields[3], fields[4]))
    return result


def classify_variants(
    vcf_path,
    classification,
    sample,
    caller,
    variant_type,
    stratification_beds,
    tmpdir,
    query_metrics,
    truth_metrics,
):
    """Classify variants from a VCF by stratification category, annotating with quality metrics."""
    rows = []

    # Convert VCF to BED
    variant_bed = os.path.join(tmpdir, f"{classification}.bed")
    vcf_to_bed(vcf_path, variant_bed)

    # Check if empty
    if os.path.getsize(variant_bed) == 0:
        return rows

    # Get all variants
    all_variants = set()
    with open(variant_bed) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            all_variants.add((fields[0], fields[1], fields[2], fields[3], fields[4]))

    def make_row(chrom, start, ref, alt, category):
        pos = str(int(start) + 1)  # Convert back to 1-based
        row = {
            "sample": sample,
            "caller": caller,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_type": variant_type,
            "classification": classification,
            "category": category,
        }
        # Add query metrics (ClairS) for TP and FP
        key = (chrom, pos, ref, alt)
        if classification in ("TP", "FP") and query_metrics:
            m = query_metrics.get(key, {})
            row["qual"] = m.get("qual", "")
            row["gq"] = m.get("gq", "")
            row["dp"] = m.get("dp", "")
            row["af"] = m.get("af", "")
            row["naf"] = m.get("naf", "")
            row["ndp"] = m.get("ndp", "")
            row["haplotype_support"] = m.get("haplotype_support", "")
        else:
            row["qual"] = ""
            row["gq"] = ""
            row["dp"] = ""
            row["af"] = ""
            row["naf"] = ""
            row["ndp"] = ""
            row["haplotype_support"] = ""
        # Add truth QUAL for TP and FN
        if classification in ("TP", "FN") and truth_metrics:
            m = truth_metrics.get(key, {})
            row["truth_qual"] = m.get("qual", "")
        else:
            row["truth_qual"] = ""
        return row

    # For each stratification BED, find overlapping variants
    categorized = set()
    for category, bed_path in stratification_beds:
        overlapping = intersect_bed(variant_bed, bed_path)
        for chrom, start, end, ref, alt in overlapping:
            rows.append(make_row(chrom, start, ref, alt, category))
            categorized.add((chrom, start, end, ref, alt))

    # Variants not in any category get "Unclassified"
    for chrom, start, end, ref, alt in all_variants - categorized:
        rows.append(make_row(chrom, start, ref, alt, "Unclassified"))

    return rows


def parse_stratification_manifest(manifest_path):
    """Parse CSV manifest of category,bed pairs."""
    beds = []
    with open(manifest_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            bed_path = row["bed"]
            # If the path doesn't exist, try finding by basename in current directory
            if not os.path.exists(bed_path):
                basename = os.path.basename(bed_path)
                if os.path.exists(basename):
                    bed_path = basename
                else:
                    raise FileNotFoundError(
                        f"BED file not found: {row['bed']} (also tried {basename})"
                    )
            beds.append((row["category"], bed_path))
    return beds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tp", required=True, help="TP VCF (truth perspective, isec 0002)"
    )
    parser.add_argument(
        "--tp-query", default=None, help="TP VCF (query perspective, isec 0003)"
    )
    parser.add_argument("--fp", required=True, help="FP VCF (query-private, isec 0001)")
    parser.add_argument("--fn", required=True, help="FN VCF (truth-private, isec 0000)")
    parser.add_argument(
        "--stratification", required=True, help="Stratification manifest CSV"
    )
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--caller", required=True, help="Query caller name")
    parser.add_argument(
        "--variant-type", required=True, help="Variant type (snv/indel)"
    )
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    stratification_beds = parse_stratification_manifest(args.stratification)

    # Build quality metric lookups
    # Query metrics: from tp_query (0003, ClairS perspective of TPs) and fp (0001, ClairS FPs)
    query_metrics = {}
    if args.tp_query:
        query_metrics.update(parse_vcf_with_metrics(args.tp_query))
    query_metrics.update(parse_vcf_with_metrics(args.fp))

    # Truth metrics: from tp (0002, truth perspective of TPs) and fn (0000, truth FNs)
    truth_metrics = {}
    truth_metrics.update(parse_vcf_with_metrics(args.tp))
    truth_metrics.update(parse_vcf_with_metrics(args.fn))

    all_rows = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for vcf_path, classification in [
            (args.tp, "TP"),
            (args.fp, "FP"),
            (args.fn, "FN"),
        ]:
            rows = classify_variants(
                vcf_path,
                classification,
                args.sample,
                args.caller,
                args.variant_type,
                stratification_beds,
                tmpdir,
                query_metrics,
                truth_metrics,
            )
            all_rows.extend(rows)

    # Write output
    fieldnames = [
        "sample",
        "caller",
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_type",
        "classification",
        "category",
        "qual",
        "gq",
        "dp",
        "af",
        "naf",
        "ndp",
        "haplotype_support",
        "truth_qual",
    ]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)


if __name__ == "__main__":
    main()
