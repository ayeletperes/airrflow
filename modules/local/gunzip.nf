process GUNZIP {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(R1, stageAs: 'input/*'), path(R2, stageAs: 'input/*')

    output:
    tuple val(meta), path("${R1.simpleName}*"), path("${R2.simpleName}*")   , emit: reads
    tuple val("${task.process}"), val('gunzip'), eval('gunzip --version 2>&1 | head -n 1 | sed \'s/^.*(gzip) //; s/ Copyright.*$//\''), emit: versions_gunzip, topic: versions

    script:
    def r1_out = R1.name.endsWith('.gz') ? R1.baseName : "${R1.baseName}.fastq"
    def r2_out = R2.name.endsWith('.gz') ? R2.baseName : "${R2.baseName}.fastq"
    """
    gzip -cdf "${R1}" > "${r1_out}"
    gzip -cdf "${R2}" > "${r2_out}"

    """
}
