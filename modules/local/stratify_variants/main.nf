process STRATIFY_VARIANTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.11 bioconda::bedtools=2.31.1 bioconda::htslib=1.21"

    input:
    tuple val(meta), path(tp_vcf), path(fp_vcf), path(fn_vcf)
    path stratification_manifest
    path stratification_beds

    output:
    tuple val(meta), path("${prefix}_stratified_calls.tsv"), emit: stratified_calls
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed 's/bedtools v//'"), topic: versions, emit: versions_bedtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    stratify_variants.py \\
        --tp ${tp_vcf} \\
        --fp ${fp_vcf} \\
        --fn ${fn_vcf} \\
        --stratification ${stratification_manifest} \\
        --sample ${meta.sample} \\
        --caller ${meta.caller} \\
        --variant-type ${meta.variant} \\
        --output ${prefix}_stratified_calls.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stratified_calls.tsv
    """
}
