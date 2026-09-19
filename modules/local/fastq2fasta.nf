process FASTQ2FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta

    script:
    """
    sed -n '1~4s/^@/>/p;2~4p' $reads > ${meta.id}.fasta
    """
}
