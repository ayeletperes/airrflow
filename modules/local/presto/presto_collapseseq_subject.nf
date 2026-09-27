process PRESTO_COLLAPSESEQ_SUBJECT {
    tag "$meta.id"
    label "process_medium"

    conda "bioconda::presto=0.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/103c49b8078f59cf606995618535a988c1055c13f06d060bdb5f642c6b217fc6/data' :
        'biocontainers/presto:0.7.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_collapse-unique.fastq"), emit: reads
    path("*_command_log.txt"), emit: logs
    path("*_table.tab")
    tuple val("${task.process}"), val('presto'), eval('CollapseSeq.py --version | grep -o "[0-9][0-9.]*" | head -n 1'), emit: versions_presto, topic: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def countf = task.ext.count_field ?: 'READCOUNT'
    // The read count is carried under ANOTHER NAME across the collapse and renamed
    // back, because CollapseSeq cannot sum DUPCOUNT into itself. It generates
    // DUPCOUNT as the size of the duplicate set AFTER the copy fields are combined,
    // so '--cf DUPCOUNT --act sum' is overwritten by the group size -- and so are
    // --act min, max and set. Measured on presto 0.7.8 with two identical reads
    // carrying DUPCOUNT 5 and 3, and a third carrying 7:
    //
    //   --cf DUPCOUNT --act sum   ->  DUPCOUNT=2   and   DUPCOUNT=1
    //   rename, sum, rename back  ->  DUPCOUNT=8   and   DUPCOUNT=7
    //
    // The second is what the reads actually are. Without this a subject's counts are
    // replaced by how many of its runs happened to agree, and a locus sequenced once
    // has every count reset to 1.
    //
    // DUPCOUNT means the same thing afterwards as it did before -- reads behind this
    // sequence -- so a sample that skipped this step (the `single` branch, one run)
    // and one that did not are read the same way downstream.
    """
    cat ${reads.join(' ')} > ${meta.id}_runs.fastq
    ParseHeaders.py rename -s ${meta.id}_runs.fastq -f DUPCOUNT -k ${countf} \\
        --outname ${meta.id}_in > /dev/null
    CollapseSeq.py -s ${meta.id}_in_reheader.fastq $args \\
        --outname ${meta.id}_col --log ${meta.id}.log > ${meta.id}_command_log.txt
    ParseHeaders.py delete -s ${meta.id}_col_collapse-unique.fastq -f DUPCOUNT \\
        --outname ${meta.id}_del > /dev/null
    ParseHeaders.py rename -s ${meta.id}_del_reheader.fastq -f ${countf} -k DUPCOUNT \\
        --outname ${meta.id}_out > /dev/null
    mv ${meta.id}_out_reheader.fastq ${meta.id}_collapse-unique.fastq
    ParseLog.py -l ${meta.id}.log $args2
    rm -f ${meta.id}_runs.fastq ${meta.id}_in_reheader.fastq \\
          ${meta.id}_col_collapse-unique.fastq ${meta.id}_del_reheader.fastq
    """
}
