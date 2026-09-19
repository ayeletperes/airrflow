// Include statements

include { GUNZIP              as GUNZIP_UMI_HEADER              } from '../../modules/local/gunzip'
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_UMI_HEADER } from '../../modules/local/fastqc_postassembly'
include { FASTQ2FASTA         as FASTQ2FASTA_UMI_HEADER         } from '../../modules/local/fastq2fasta'

//PRESTO
include { PRESTO_ASSEMBLEPAIRS          as PRESTO_ASSEMBLEPAIRS_UMI_HEADER          } from '../../modules/local/presto/presto_assemblepairs'
include { PRESTO_FILTERSEQ_POSTASSEMBLY as PRESTO_FILTERSEQ_POSTASSEMBLY_UMI_HEADER } from '../../modules/local/presto/presto_filterseq_postassembly'


workflow PRESTO_UMI_HEADER {
    take:
    ch_reads    // channel: [ val(meta), [ reads ] ]
    filterseq_q

    main:
    ch_versions = channel.empty()

    // Reads arrive demultiplexed and trimmed, with the UMI and the trimmed prefix carried
    // as UMI= / TRIM= annotations in the R1 header, so no adapter trimming is needed.
    GUNZIP_UMI_HEADER ( ch_reads.map{ meta, reads -> [meta, reads[0], reads[1]] } )

    // Assemble read pairs, propagating the R1 header annotations onto the assembled read
    PRESTO_ASSEMBLEPAIRS_UMI_HEADER (
        GUNZIP_UMI_HEADER.out.reads
    )

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_POSTASSEMBLY_UMI_HEADER (
        PRESTO_ASSEMBLEPAIRS_UMI_HEADER.out.reads,
        filterseq_q
    )

    FASTQC_POSTASSEMBLY_UMI_HEADER (
        PRESTO_FILTERSEQ_POSTASSEMBLY_UMI_HEADER.out.reads
    )

    // ponytail: no PRESTO_PARSEHEADERS_METADATA here - the DSL1 RM protocol carries sample
    // metadata in the samplesheet, not in the read headers. Revisit if downstream steps
    // start reading it off the FASTA description.
    FASTQ2FASTA_UMI_HEADER (
        PRESTO_FILTERSEQ_POSTASSEMBLY_UMI_HEADER.out.reads
    )

    emit:
    fasta = FASTQ2FASTA_UMI_HEADER.out.fasta
    versions = ch_versions
    fastqc_postassembly_gz = FASTQC_POSTASSEMBLY_UMI_HEADER.out.zip
    presto_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_UMI_HEADER.out.logs.collect()
    presto_filterseq_logs = PRESTO_FILTERSEQ_POSTASSEMBLY_UMI_HEADER.out.logs
}
