#!/usr/bin/env python3
"""Stratify TP/FP/FN variants by genomic region using bedtools intersect.

Produces a long-format TSV with one row per variant × category overlap.
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile


def parse_vcf_records(vcf_path):
    """Parse VCF and yield (chrom, pos, ref, alt) tuples."""
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alt = (
                fields[0],
                fields[1],
                fields[2],
                fields[3],
                fields[4],
            )
            yield chrom, pos, ref, alt


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
    vcf_path, classification, sample, caller, variant_type, stratification_beds, tmpdir
):
    """Classify variants from a VCF by stratification category."""
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

    # For each stratification BED, find overlapping variants
    categorized = set()
    for category, bed_path in stratification_beds:
        overlapping = intersect_bed(variant_bed, bed_path)
        for chrom, start, end, ref, alt in overlapping:
            pos = str(int(start) + 1)  # Convert back to 1-based
            rows.append(
                {
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
            )
            categorized.add((chrom, start, end, ref, alt))

    # Variants not in any category get "Unclassified"
    for chrom, start, end, ref, alt in all_variants - categorized:
        pos = str(int(start) + 1)
        rows.append(
            {
                "sample": sample,
                "caller": caller,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "variant_type": variant_type,
                "classification": classification,
                "category": "Unclassified",
            }
        )

    return rows


def parse_stratification_manifest(manifest_path):
    """Parse CSV manifest of category,bed pairs."""
    beds = []
    with open(manifest_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            beds.append((row["category"], row["bed"]))
    return beds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tp", required=True, help="TP VCF file")
    parser.add_argument("--fp", required=True, help="FP VCF file")
    parser.add_argument("--fn", required=True, help="FN VCF file")
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
    ]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)


if __name__ == "__main__":
    main()
