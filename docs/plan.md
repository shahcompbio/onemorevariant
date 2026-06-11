# Plan: shahcompbio/onemorevariant — Somatic SNV/Indel Benchmarking Pipeline

## TL;DR

Build a Nextflow DSL2 pipeline (nf-core style) that benchmarks somatic SNV/indel calls from one or more query callers against a constructed truth set, using bcftools isec for TP/FP/FN classification and bedtools intersect for region stratification. Designed for ONT vs Illumina comparison but generic enough for any caller combination.

## Pipeline Overview

**Input:** Samplesheet (sample, variant, caller, category, vcf) + reference FASTA + stratification BED manifest + truth construction params.

**Core logic:**

1. Preprocess: PASS filter, contig filter, normalize (split multiallelic + left-align)
2. Build truth: Combine truth-category VCFs per sample × variant_type using configurable strategy (union/intersect/single caller)
3. Benchmark: bcftools isec (truth vs each query caller) → classify TP/FP/FN
4. Stratify: bedtools intersect TP/FP/FN against each BED → per-variant category annotations
5. Aggregate: summary stats + stratified_calls tables, separate for SNVs and indels

## Steps

### Phase 1: Input & Samplesheet

1. **Define samplesheet schema** (`assets/schema_input.json`). Columns:
   - `sample` (string, required) — sample ID
   - `variant` (string, required, enum: "snv", "indel") — variant type
   - `caller` (string, required) — caller name (e.g., "mutect2", "strelka", "clairs", "deepsomatic")
   - `category` (string, required, enum: "query", "truth") — whether this VCF is a query or truth input
   - `vcf` (file-path, required) — path to VCF file

   Example samplesheet:

   ```csv
   sample,variant,caller,category,vcf
   SHAH_H002190_T02,snv,mutect2,truth,/path/to/mutect2_snv.vcf.gz
   SHAH_H002190_T02,snv,strelka,truth,/path/to/strelka_snv.vcf.gz
   SHAH_H002190_T02,snv,clairs,query,/path/to/clairs_snv.vcf.gz
   SHAH_H002190_T02,indel,strelka,truth,/path/to/strelka_indel.vcf.gz
   SHAH_H002190_T02,indel,clairs,query,/path/to/clairs_indel.vcf.gz
   ```

2. **Define stratification BED manifest schema** (`assets/schema_stratification.json`). Columns:
   - `category` (string) — stratification category name (e.g., "segdup", "tandem_repeat")
   - `bed` (file-path) — path to BED file (.bed or .bed.gz)

3. **Add pipeline params** to `nextflow.config` / `nextflow_schema.json`:
   - `--input` — samplesheet CSV
   - `--fasta` — reference FASTA (required)
   - `--stratification` — path to stratification BED manifest CSV (optional)
   - `--contigs` — comma-separated contig list (default: chr1,...,chr22,chrX,chrY)
   - `--snv_truth` — truth construction strategy for SNVs (default: "union"; options: "union", "intersect", or a caller name)
   - `--indel_truth` — truth construction strategy for indels (default: "strelka"; options: "union", "intersect", or a caller name)

### Phase 2: Preprocessing (prefer nf-core modules)

4. **PREPROCESS_VCF** — applied to every samplesheet row (truth and query alike)
   - Use nf-core modules where available: `bcftools/view`, `bcftools/norm`, `tabix/tabix`
   - Logic: `bcftools view -f PASS -t ${contigs}` → `bcftools norm -m - -f ${fasta}` → `tabix -p vcf`
   - Output: preprocessed VCF + index

### Phase 3: Truth Construction

5. **BUILD_TRUTH** — groups preprocessed truth VCFs by sample × variant_type, then applies strategy:
   - Use nf-core modules where available: `bcftools/concat`, `bcftools/sort`, `bcftools/norm`, `bcftools/isec`
   - Logic:
     - `--snv_truth union` → `bcftools concat -a` → sort → `bcftools norm -d exact` (dedup) → index
     - `--snv_truth intersect` → `bcftools isec -n =N` (keep positions in ALL N truth VCFs) → index
     - `--snv_truth <caller_name>` → pass through that caller's VCF only
   - Output: `truth_{snv|indel}.vcf.gz` + `.tbi`

### Phase 4: Benchmarking

6. **BCFTOOLS_ISEC** — each query caller benchmarked independently against the same truth
   - Use nf-core `bcftools/isec` module if available, otherwise local module
   - Logic: `bcftools isec -p ${prefix}_isec ${truth_vcf} ${query_vcf}`
     - `0000.vcf` = truth-only = FN
     - `0001.vcf` = query-only = FP
     - `0002.vcf` = shared (truth side) = TP
   - Output: `tp.vcf.gz`, `fp.vcf.gz`, `fn.vcf.gz` (bgzipped + indexed)

### Phase 5: Stratification

7. **STRATIFY_VARIANTS** — local module (no nf-core equivalent)
   - Input: meta (sample, variant_type, query_caller) + tp_vcf + fp_vcf + fn_vcf + stratification_beds (from manifest)
   - Logic: Python/bash script that:
     1. For each variant in TP/FP/FN, intersect against each stratification BED using `bedtools intersect -u`
     2. Produce `stratified_calls.tsv` in long format:
        - Columns: `sample`, `caller`, `chrom`, `pos`, `ref`, `alt`, `variant_type`, `classification` (TP/FP/FN), `category`
        - Variant with no category overlaps → `category = Unclassified`
        - Variant overlapping multiple BEDs → multiple rows (one per category)
   - Output: `stratified_calls.tsv`
   - Container: bcftools + bedtools + python3

### Phase 6: Aggregation

8. **AGGREGATE_RESULTS** — local module, run separately for SNVs and indels
   - Input: All per-sample stratified_calls.tsv files for a given variant_type (collected)
   - Logic: Python script that:
     1. Concatenates all stratified_calls.tsv
     2. Computes per sample × caller × category: TP, FP, FN, precision, recall, F1
     3. Also computes genome-wide (category = "genome_wide") stats per sample × caller
   - Output (separate files for SNV and indel):
     - `snv_benchmark_summary.tsv` / `indel_benchmark_summary.tsv` — one row per sample × caller × region
     - `snv_stratified_calls.tsv` / `indel_stratified_calls.tsv` — concatenated per-variant table
   - Container: python3 + pandas

### Phase 7: Workflow Wiring

9. **Wire modules in `workflows/onemorevariant.nf`**
   - Parse samplesheet → channel grouped by (sample, variant_type, category)
   - Parse stratification manifest → channel of (category, bed_path) tuples
   - Run PREPROCESS_VCF on every samplesheet row (parallel)
   - Group preprocessed truth VCFs by sample × variant_type → BUILD_TRUTH (depends on preprocess)
   - Cross each query VCF with its corresponding truth → BCFTOOLS_ISEC (parallel per query caller)
   - Run STRATIFY_VARIANTS per sample × variant_type × query_caller (depends on isec)
   - Branch by variant_type → AGGREGATE_RESULTS separately for SNVs and indels (depends on all stratify)
   - Collect versions for MultiQC

10. **Update `conf/modules.config`** with publishDir for each module.

11. **Update `nextflow_schema.json`** with new params (fasta, stratification, contigs, snv_truth, indel_truth).

### Phase 8: Containers & Testing

12. **Container strategy**: Use nf-core module containers for bcftools steps. Create a lightweight custom container (or conda env) for the Python stratification + aggregation scripts (pandas + bedtools).

13. **Test data**: Create minimal test VCFs (a few variants on chr21) + a small stratification BED for CI testing under `tests/`.

14. **nf-test**: Write basic pipeline test in `tests/default.nf.test`.

## Relevant Files

- `workflows/onemorevariant.nf` — main workflow logic (rewrite from boilerplate)
- `main.nf` — entry point (already wired to workflow)
- `assets/schema_input.json` — samplesheet validation (rewrite from fastq-based)
- `assets/samplesheet.csv` — example samplesheet (rewrite)
- `nextflow.config` — pipeline params
- `nextflow_schema.json` — param schema for nf-core
- `conf/modules.config` — per-module publishDir config
- `conf/base.config` — resource labels
- `modules/` — nf-core modules to install + local modules to create
- `bin/` — custom scripts: `stratify_variants.py`, `aggregate_results.py`

**Reference implementations:**

- `/Users/preskaa/VSCodeProjects/onemoresv/modules/local/stratify/main.nf` — STRATIFY module pattern
- `/Users/preskaa/VSCodeProjects/onemoresv/bin/minda_stratify.py` — stratified_calls.tsv output format
- `/Users/preskaa/PycharmProjects/SarcAtlas/analyses/260605_APS066_ont_v_ill_snv/scripts/APS066_generate_truth_vcfs.sh` — existing bash logic

## Verification

1. Run `nextflow run . -profile test,docker --outdir results` with minimal test data
2. Verify `snv_benchmark_summary.tsv` and `indel_benchmark_summary.tsv` match expected precision/recall
3. Verify `snv_stratified_calls.tsv` has correct variant × category assignments
4. Test with multiple query callers — ensure each gets independent benchmarking
5. Test truth strategies: verify "union" produces superset, "intersect" produces subset, caller name selects correctly
6. Run on one real sample (SHAH_H002190_T02) and compare to lab notebook pilot stats

## Decisions

| Decision               | Choice                                                                          |
| ---------------------- | ------------------------------------------------------------------------------- |
| Benchmarking tool      | `bcftools isec`                                                                 |
| Samplesheet format     | Generic long (sample, variant, caller, category, vcf) — future-proof            |
| Truth construction     | Parameterized: `--snv_truth` (default union), `--indel_truth` (default strelka) |
| Multiple query callers | Supported — each benchmarked independently                                      |
| Stratification         | Post-hoc bedtools intersect; CSV manifest of BEDs                               |
| Output format          | Long-format stratified_calls + summary, separate SNV/indel files                |
| Normalization          | `bcftools norm -m - -f`                                                         |
| Tumor ID               | Not needed (position-based matching)                                            |
| Module preference      | Use nf-core modules over custom where available                                 |

## Further Considerations

1. **Deduplication in truth union**: After `bcftools concat` of multiple truth callers, the same variant position may appear with identical REF/ALT from both callers. `bcftools norm -d exact` in BUILD_TRUTH handles this.
