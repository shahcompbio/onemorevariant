process AGGREGATE_RESULTS {
    tag "${variant_type}"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.2"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(variant_type), path(stratified_calls)

    output:
    tuple val(variant_type), path("${variant_type}_benchmark_summary.tsv"), emit: summary
    tuple val(variant_type), path("${variant_type}_stratified_calls.tsv"), emit: stratified_calls
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    aggregate_results.py \\
        --input ${stratified_calls} \\
        --variant-type ${variant_type} \\
        --output-summary ${variant_type}_benchmark_summary.tsv \\
        --output-calls ${variant_type}_stratified_calls.tsv
    """

    stub:
    """
    touch ${variant_type}_benchmark_summary.tsv
    touch ${variant_type}_stratified_calls.tsv
    """
}
