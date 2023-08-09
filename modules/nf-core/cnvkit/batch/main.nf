process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::cnvkit=0.9.10 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-780d630a9bb6a0ff2e7b6f730906fd703e40e98f:c94363856059151a2974dc501fb07a0360cc60a3-0' :
        'biocontainers/mulled-v2-780d630a9bb6a0ff2e7b6f730906fd703e40e98f:c94363856059151a2974dc501fb07a0360cc60a3-0' }"

    input:
    tuple val(meta), path(normal)
    tuple val(meta2), path(fasta)
    path  targets

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn, optional: true
    tuple val(meta), path("*.cnr"), emit: cnr, optional: true
    tuple val(meta), path("*.cns"), emit: cns, optional: true
    tuple val(meta), path("*.pdf"), emit: pdf, optional: true
    tuple val(meta), path("*.png"), emit: png, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // PON needs --normal, target files if they exists
    fasta_args = fasta ? "--fasta ${fasta}" : ""
    target_args = targets ? "--targets ${targets}" : ""

    pon_input = normal.collect().join(' ')
    normal_args = "--normal $pon_input"

    def versions =
        "cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')"
    """

    cnvkit.py \\
        batch \\
        $normal_args \\
        $fasta_args \\
        $target_args \\
        --processes $task.cpus \\
        $args \\
        --output-reference ${prefix}.cnn

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${versions}
    END_VERSIONS
    """
}
