// Include statements

include { GUNZIP              as GUNZIP_TS_HEADER              } from '../../modules/local/gunzip'
include { FASTQC_POSTASSEMBLY as FASTQC_POSTASSEMBLY_TS_HEADER } from '../../modules/local/fastqc_postassembly'
include { FASTQ2FASTA         as FASTQ2FASTA_TS_HEADER         } from '../../modules/local/fastq2fasta'

//PRESTO
include { PRESTO_ASSEMBLEPAIRS          as PRESTO_ASSEMBLEPAIRS_TS_HEADER          } from '../../modules/local/presto/presto_assemblepairs'
include { PRESTO_FILTERSEQ_POSTASSEMBLY as PRESTO_FILTERSEQ_POSTASSEMBLY_TS_HEADER } from '../../modules/local/presto/presto_filterseq_postassembly'
include { PRESTO_COLLAPSESEQ            as PRESTO_COLLAPSESEQ_TS_HEADER            } from '../../modules/local/presto/presto_collapseseq'
include { PRESTO_COLLAPSESEQ_SUBJECT    as PRESTO_COLLAPSESEQ_SUBJECT_TS_HEADER    } from '../../modules/local/presto/presto_collapseseq_subject'


workflow PRESTO_TS_HEADER {
    take:
    ch_reads    // channel: [ val(meta), [ reads ] ]
    filterseq_q

    main:
    ch_versions = channel.empty()

    // Reads arrive demultiplexed and trimmed, with the UMI and the trimmed prefix carried
    // as UMI= / TRIM= annotations in the R1 header, so no adapter trimming is needed.
    GUNZIP_TS_HEADER ( ch_reads.map{ meta, reads -> [meta, reads[0], reads[1]] } )

    // Assemble read pairs, propagating the R1 header annotations onto the assembled read
    PRESTO_ASSEMBLEPAIRS_TS_HEADER (
        GUNZIP_TS_HEADER.out.reads
    )

    // Filter sequences by quality score
    PRESTO_FILTERSEQ_POSTASSEMBLY_TS_HEADER (
        PRESTO_ASSEMBLEPAIRS_TS_HEADER.out.reads,
        filterseq_q
    )

    FASTQC_POSTASSEMBLY_TS_HEADER (
        PRESTO_FILTERSEQ_POSTASSEMBLY_TS_HEADER.out.reads
    )

    // Collapse identical sequences before anything has to align them. On the macaque
    // libraries this removes a third of the reads, and V(D)J annotation is the most
    // expensive step in the pipeline by a wide margin. The count of collapsed reads
    // travels on as DUPCOUNT.
    PRESTO_COLLAPSESEQ_TS_HEADER (
        PRESTO_FILTERSEQ_POSTASSEMBLY_TS_HEADER.out.reads
    )

    // One subject's locus can be sequenced several times. Assembly and quality filtering
    // have to stay per run, or a failed run disappears into the averages of the merged
    // one, but by here nothing is run specific: those are the same molecules read twice.
    // Collapse them together before they become FASTA, so annotation does not pay for
    // the duplicates and clone sizes and genotype depth do not count them twice.
    // DUPCOUNT is summed, so a sequence read 3 times in one run and 5 in another stands
    // for 8 reads.
    PRESTO_COLLAPSESEQ_TS_HEADER.out.reads
        .map { meta, reads -> [ [ meta.subject_id, meta.locus ], meta, reads ] }
        .groupTuple()
        .branch { _key, metas, _reads ->
            single: metas.size() == 1
            runs:   metas.size() > 1
        }
        .set { ch_by_subject_locus }

    PRESTO_COLLAPSESEQ_SUBJECT_TS_HEADER (
        ch_by_subject_locus.runs.map { key, metas, reads ->
            [ metas[0] + [ id: "${key[0]}_${key[1]}" ], reads ]
        }
    )

    ch_collapsed = ch_by_subject_locus.single
        .map { _key, metas, reads -> [ metas[0], reads[0] ] }
        .mix( PRESTO_COLLAPSESEQ_SUBJECT_TS_HEADER.out.reads )

    // ponytail: no PRESTO_PARSEHEADERS_METADATA here - the DSL1 RM protocol carries sample
    // metadata in the samplesheet, not in the read headers. Revisit if downstream steps
    // start reading it off the FASTA description.
    FASTQ2FASTA_TS_HEADER ( ch_collapsed )

    emit:
    fasta = FASTQ2FASTA_TS_HEADER.out.fasta
    versions = ch_versions
    fastqc_postassembly_gz = FASTQC_POSTASSEMBLY_TS_HEADER.out.zip
    presto_assemblepairs_logs = PRESTO_ASSEMBLEPAIRS_TS_HEADER.out.logs.collect()
    presto_filterseq_logs = PRESTO_FILTERSEQ_POSTASSEMBLY_TS_HEADER.out.logs
    presto_collapseseq_logs = PRESTO_COLLAPSESEQ_TS_HEADER.out.logs
}
