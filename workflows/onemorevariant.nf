/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { BCFTOOLS_VIEW                        } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_NORM                        } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_CONCAT                              } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_CONSENSUS } from '../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT                                } from '../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_CONSENSUS     } from '../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_ISEC                                } from '../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_ISEC as BCFTOOLS_ISEC_CONSENSUS     } from '../modules/nf-core/bcftools/isec/main'
include { BCFTOOLS_REHEADER                            } from '../modules/nf-core/bcftools/reheader/main'
include { STRATIFY_VARIANTS                            } from '../modules/local/stratify_variants/main'
include { AGGREGATE_RESULTS as AGGREGATE_SNV           } from '../modules/local/aggregate_results/main'
include { AGGREGATE_RESULTS as AGGREGATE_INDEL         } from '../modules/local/aggregate_results/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_onemorevariant_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONEMOREVARIANT {
    take:
    ch_samplesheet // channel: [ [meta], vcf ]
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_multiqc_files = channel.empty()

    // Prepare reference fasta channel
    def ch_fasta = channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])

    //
    // STEP 1: PREPROCESS — bcftools view (PASS filter + contig filter) on all VCFs
    //
    // Input channel: [ meta, vcf, index ] — no index available at this stage
    def ch_view_input = ch_samplesheet.map { meta, vcf ->
        [meta, vcf, []]
    }

    BCFTOOLS_VIEW(ch_view_input, [], [], [])

    //
    // STEP 2: NORMALIZE — bcftools norm (split multiallelic + left-align) on filtered VCFs
    //
    def ch_norm_input = BCFTOOLS_VIEW.out.vcf.map { meta, vcf ->
        [meta, vcf, []]
    }

    BCFTOOLS_NORM(ch_norm_input, ch_fasta)

    //
    // STEP 2b: REHEADER — standardize sample names to meta.id
    //
    def ch_reheader_input = BCFTOOLS_NORM.out.vcf.map { meta, vcf ->
        def samples_file = file("${workDir}/reheader/${meta.id}_${meta.caller}_${meta.variant}_samples.txt")
        samples_file.parent.mkdirs()
        samples_file.text = "${meta.id}\n"
        [meta, vcf, [], samples_file]
    }

    BCFTOOLS_REHEADER(ch_reheader_input, [[], []])

    //
    // STEP 3: BUILD TRUTH — group truth VCFs by sample × variant_type, apply strategy
    //
    // Separate truth and query channels — join vcf with index
    def ch_preprocessed = BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index, by: [0])

    def ch_truth = ch_preprocessed.filter { meta, _vcf, _tbi -> meta.category == 'truth' }

    def ch_query = ch_preprocessed.filter { meta, _vcf, _tbi -> meta.category == 'query' }

    // Group truth VCFs by sample × variant_type for truth construction
    // Key: [sample, variant] -> list of [meta, vcf] entries
    def ch_truth_grouped = ch_truth
        .map { meta, vcf, tbi -> [[id: meta.id, variant: meta.variant], meta, vcf, tbi] }
        .groupTuple(by: 0)
        .map { group_key, metas, vcfs, tbis ->
            // Determine strategy based on variant type
            def strategy = group_key.variant == 'snv' ? params.snv_truth : params.indel_truth
            def callers = metas.collect { it.caller }

            if (strategy == 'intersect') {
                // Intersect — will be handled by isec downstream
                return [[id: group_key.id, variant: group_key.variant, _intersect: true], vcfs, tbis]
            }
            else if (strategy == 'consensus') {
                // Consensus — >=2/N agreement via pairwise intersections
                return [[id: group_key.id, variant: group_key.variant, _consensus: true], vcfs, tbis]
            }
            else if (strategy != 'union' && strategy in callers) {
                // Single caller strategy — find matching VCF
                def idx = callers.indexOf(strategy)
                return [[id: group_key.id, variant: group_key.variant], [vcfs[idx]], [tbis[idx]]]
            }
            else {
                // Union strategy (default) — also used as fallback if caller name not found
                return [[id: group_key.id, variant: group_key.variant], vcfs, tbis]
            }
        }

    // For consensus strategy: bcftools isec -n+2 → concat → sort
    def ch_truth_consensus = ch_truth_grouped
        .filter { meta, _vcfs, _tbis -> meta.containsKey('_consensus') }
        .map { meta, vcfs, tbis -> [meta.subMap('id', 'variant'), vcfs, tbis, [], [], []] }

    BCFTOOLS_ISEC_CONSENSUS(ch_truth_consensus)

    // Extract numbered VCFs from isec output directory
    def ch_consensus_isec_vcfs = BCFTOOLS_ISEC_CONSENSUS.out.results.map { meta, dir ->
        def vcfs = file("${dir}/*.vcf.gz").sort()
        def tbis = vcfs.collect { vcf -> file("${vcf}.tbi") }
        [meta, vcfs, tbis]
    }

    BCFTOOLS_CONCAT_CONSENSUS(ch_consensus_isec_vcfs)
    BCFTOOLS_SORT_CONSENSUS(BCFTOOLS_CONCAT_CONSENSUS.out.vcf)

    def ch_consensus_truth_out = BCFTOOLS_SORT_CONSENSUS.out.vcf.join(BCFTOOLS_SORT_CONSENSUS.out.index, by: [0])

    // For union strategy: concat (with --remove-duplicates) → sort
    def ch_truth_to_concat = ch_truth_grouped
        .filter { meta, vcfs, _tbis -> !meta.containsKey('_intersect') && !meta.containsKey('_consensus') && vcfs.size() > 1 }
        .map { meta, vcfs, tbis -> [meta, vcfs, tbis] }

    def ch_truth_single = ch_truth_grouped
        .filter { meta, vcfs, _tbis -> !meta.containsKey('_intersect') && !meta.containsKey('_consensus') && vcfs.size() == 1 }
        .map { meta, vcfs, tbis -> [meta, vcfs[0], tbis[0]] }

    BCFTOOLS_CONCAT(ch_truth_to_concat)
    BCFTOOLS_SORT(BCFTOOLS_CONCAT.out.vcf)

    // Combine single-VCF truth with concat+sort multi-VCF truth (with indexes)
    def ch_truth_sort_with_idx = BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.index, by: [0])

    def ch_truth_final = ch_truth_single.mix(ch_truth_sort_with_idx).mix(ch_consensus_truth_out)

    //
    // STEP 4: BENCHMARK — bcftools isec (truth vs query)
    //
    // Join query VCFs with their corresponding truth by [sample, variant]
    def ch_query_keyed = ch_query.map { meta, vcf, tbi ->
        [[id: meta.id, variant: meta.variant], meta, vcf, tbi]
    }

    def ch_truth_keyed = ch_truth_final.map { meta, vcf, tbi ->
        [[id: meta.id, variant: meta.variant], vcf, tbi]
    }

    // Combine: each query gets paired with its truth
    def ch_isec_input = ch_query_keyed
        .combine(ch_truth_keyed, by: 0)
        .map { _group_key, query_meta, query_vcf, query_tbi, truth_vcf, truth_tbi ->
            // bcftools isec input: [meta, [vcfs], [tbis], file_list, targets, regions]
            def meta = [
                id: "${query_meta.id}_${query_meta.caller}_${query_meta.variant}",
                sample: query_meta.id,
                caller: query_meta.caller,
                variant: query_meta.variant,
            ]
            [meta, [truth_vcf, query_vcf], [truth_tbi, query_tbi], [], [], []]
        }

    BCFTOOLS_ISEC(ch_isec_input)

    //
    // STEP 5: STRATIFY — bedtools intersect TP/FP/FN against stratification BEDs
    //
    if (params.stratification) {
        // Prepare stratification manifest and stage BED files
        def manifest_file = file(params.stratification, checkIfExists: true)
        def ch_stratification = channel.value(manifest_file)

        // Parse the manifest CSV to collect BED file paths for staging
        def ch_stratification_beds = channel.value(
            manifest_file.readLines().drop(1).collect { line ->
                file(line.split(',')[1].trim(), checkIfExists: true)
            }
        )

        // Extract TP, FP, FN VCFs from isec results directory
        def ch_isec_results = BCFTOOLS_ISEC.out.results.map { meta, results_dir ->
            def tp = file("${results_dir}/0002.vcf.gz").exists()
                ? file("${results_dir}/0002.vcf.gz")
                : file("${results_dir}/0002.vcf")
            def tp_query = file("${results_dir}/0003.vcf.gz").exists()
                ? file("${results_dir}/0003.vcf.gz")
                : file("${results_dir}/0003.vcf")
            def fp = file("${results_dir}/0001.vcf.gz").exists()
                ? file("${results_dir}/0001.vcf.gz")
                : file("${results_dir}/0001.vcf")
            def fn = file("${results_dir}/0000.vcf.gz").exists()
                ? file("${results_dir}/0000.vcf.gz")
                : file("${results_dir}/0000.vcf")
            [meta, tp, tp_query, fp, fn]
        }

        STRATIFY_VARIANTS(ch_isec_results, ch_stratification, ch_stratification_beds)

        //
        // STEP 6: AGGREGATE — collect per variant_type and compute summary stats
        //
        def ch_stratified_snv = STRATIFY_VARIANTS.out.stratified_calls
            .filter { meta, _calls -> meta.variant == 'snv' }
            .map { _meta, calls -> calls }
            .collect()
            .filter { calls -> calls.size() > 0 }
            .map { calls -> ['snv', calls] }

        def ch_stratified_indel = STRATIFY_VARIANTS.out.stratified_calls
            .filter { meta, _calls -> meta.variant == 'indel' }
            .map { _meta, calls -> calls }
            .collect()
            .filter { calls -> calls.size() > 0 }
            .map { calls -> ['indel', calls] }

        AGGREGATE_SNV(ch_stratified_snv)
        AGGREGATE_INDEL(ch_stratified_indel)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'onemorevariant_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'onemorevariant'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
}
