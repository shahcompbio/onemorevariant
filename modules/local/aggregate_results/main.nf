process AGGREGATE_RESULTS {
    tag "${variant_type}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/onemorevariant-aggregate:python-3.11_pandas-2.2--0f9e990626d9afe0"

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
