# Plan: Add `consensus` truth mode to onemorevariant

## Summary

Add a `consensus` truth construction strategy to onemorevariant that builds truth from variants called by **≥2 of N** input truth callers (majority vote). This enables orthogonal-tools benchmarking as described in the DeepSomatic paper (Samarakoon et al. 2025), where cross-platform caller agreement defines the truth set.

## Motivation

The current `union` and `intersect` strategies have limitations for cross-platform benchmarking:
- **`union`** (any caller) is too permissive — includes low-confidence single-caller variants
- **`intersect`** (all callers) is too strict — requires agreement from every caller, which is rare across platforms with different error profiles

**`consensus` (≥2/N)** is the sweet spot: a variant needs support from at least two independent callers/platforms to enter the truth set. This:
- Reduces platform bias (ONT-confirmed variants get credit if one Illumina caller also sees them)
- Avoids circularity (a single caller can't place variants into truth alone)
- Matches established benchmarking practice in the field

## Use Case (APS066)

Benchmark ClairS (ONT) somatic SNV/indel calls using a consensus truth:
- **SNV truth:** ≥2 of {Mutect2, Strelka2, ClairS} agree
- **Indel truth:** ≥2 of {Strelka2, ClairS} agree (equivalent to Strelka2 ∩ ClairS for N=2)

**Samplesheet:** ClairS VCFs are listed as **both** `truth` and `query`:

```csv
sample,variant,caller,category,vcf
SHAH_H002190_T02,snv,mutect2,truth,/path/to/mutect2_snv.vcf.gz
SHAH_H002190_T02,snv,strelka,truth,/path/to/strelka_snv.vcf.gz
SHAH_H002190_T02,snv,clairs,truth,/path/to/clairs_snv.vcf.gz
SHAH_H002190_T02,snv,clairs,query,/path/to/clairs_snv.vcf.gz
SHAH_H002190_T02,indel,strelka,truth,/path/to/strelka_indel.vcf.gz
SHAH_H002190_T02,indel,clairs,truth,/path/to/clairs_indel.vcf.gz
SHAH_H002190_T02,indel,clairs,query,/path/to/clairs_indel.vcf.gz
```

Launch: `--snv_truth consensus --indel_truth consensus`

## Implementation

### Changes to `workflows/onemorevariant.nf`

In the truth grouping logic (the `.map` that applies the strategy), add a `consensus` branch:

```groovy
if (strategy == 'consensus') {
    // Pass all truth VCFs through to a new CONSENSUS_TRUTH process
    return [[id: group_key.id, variant: group_key.variant, _consensus: true], vcfs, tbis]
}
```

### New Process: `CONSENSUS_TRUTH`

**Location:** `modules/local/consensus_truth/main.nf`

**Input:** meta + list of N normalized truth VCFs + their indices

**Logic (bash/bcftools):**
1. Compute all pairwise intersections using `bcftools isec -n=2 -w1`:
   - For N truth VCFs, this is N×(N-1)/2 pairwise comparisons
   - Each produces variants from VCF_i that are also in VCF_j (positional match)
2. Union all pairwise intersection VCFs: `bcftools concat --allow-overlaps --remove-duplicates`
3. Sort: `bcftools sort`
4. Index: `bcftools index -t`

**For N=2:** Consensus = intersect (just one pairwise comparison). This is fine — `consensus` with 2 callers degenerates to `intersect`, which is the expected behavior for the indel case (Strelka2 ∩ ClairS).

**For N=3:** Three pairwise intersections → union = ≥2/3.

**Output:** `consensus_truth.vcf.gz` + `.tbi`

**Script sketch:**
```bash
#!/bin/bash
set -euo pipefail

VCFS=($@)  # all input VCFs passed as positional args
N=${#VCFS[@]}
ISEC_DIR="pairwise_isec"
mkdir -p ${ISEC_DIR}

# Compute all pairwise intersections
PAIR_IDX=0
ISEC_FILES=()
for ((i=0; i<N-1; i++)); do
    for ((j=i+1; j<N; j++)); do
        OUT="${ISEC_DIR}/pair_${PAIR_IDX}.vcf.gz"
        bcftools isec -n=2 -w1 -Oz -o ${OUT} ${VCFS[$i]} ${VCFS[$j]}
        bcftools index -t ${OUT}
        ISEC_FILES+=("${OUT}")
        PAIR_IDX=$((PAIR_IDX + 1))
    done
done

# Union of pairwise intersections = ≥2/N consensus
bcftools concat --allow-overlaps --remove-duplicates "${ISEC_FILES[@]}" | \
    bcftools sort -Oz -o consensus_truth.vcf.gz
bcftools index -t consensus_truth.vcf.gz
```

### Workflow Wiring

After the strategy routing, add:

```groovy
// Filter consensus-strategy truth groups
def ch_truth_consensus = ch_truth_grouped
    .filter { meta, vcfs, _tbis -> meta.containsKey('_consensus') }
    .map { meta, vcfs, tbis -> [meta.subMap('id', 'variant'), vcfs, tbis] }

CONSENSUS_TRUTH(ch_truth_consensus)

// Merge consensus truth output into the unified truth channel
// (alongside union/intersect/single-caller results)
```

### Config Changes

Update `nextflow.config` params docs and `nextflow_schema.json`:
- `--snv_truth` options: `"union"`, `"intersect"`, `"consensus"`, or a caller name
- `--indel_truth` options: same

### Testing

1. **Unit test (nf-test):** 3 small VCFs with known overlaps on chr21. Verify consensus output has exactly the ≥2/3 variants.
2. **Integration test:** Run on SHAH_H002190_T02 with the consensus samplesheet. Compare variant counts to manually computed values.
3. **Edge case (N=2):** Verify consensus with 2 truth callers produces same result as intersect.

## Files to Modify/Create

| File | Action |
|------|--------|
| `workflows/onemorevariant.nf` | Add `consensus` strategy branch + wire CONSENSUS_TRUTH process |
| `modules/local/consensus_truth/main.nf` | New process (pairwise isec → union) |
| `nextflow.config` | Update param docs |
| `nextflow_schema.json` | Add "consensus" to enum |
| `docs/usage.md` | Document consensus mode with example samplesheet |
| `conf/test_consensus.config` | Test profile for consensus mode |

## Verification Criteria

- [ ] `--snv_truth consensus` with 3 truth callers produces ≥2/3 consensus truth
- [ ] `--indel_truth consensus` with 2 truth callers produces their intersection
- [ ] Query caller listed as both truth and query works correctly (same VCF used in both roles)
- [ ] Pipeline produces correct precision/recall when ClairS contributes to truth (precision should increase vs Illumina-only truth)
- [ ] Existing `union`/`intersect`/single-caller modes still work (no regression)
