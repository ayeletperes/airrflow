include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'
include { VALIDATE_IGBLAST_DB } from '../../modules/local/validate_igblast_db'

workflow DATABASES {

    take:

    main:
    ch_versions = Channel.empty()
    ch_airrdb_source = Channel.value(params.airrdb_source)

    if( !params.fetch_database ){
        if (params.reference_igblast.endsWith(".zip")) {
            Channel.fromPath("${params.reference_igblast}")
                    .ifEmpty{ error "IGBLAST DB not found: ${params.reference_igblast}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_IGBLAST.out.versions)
        } else {
            Channel.fromPath("${params.reference_igblast}")
                .ifEmpty { error "IGBLAST DB not found: ${params.reference_igblast}" }
                .set { ch_igblast }
        }
    }

    if( !params.fetch_database ){
        if (params.reference_fasta.endsWith(".zip")) {
            Channel.fromPath("${params.reference_fasta}")
                    .ifEmpty{ error "IMGTDB not found: ${params.reference_fasta}" }
                    .set { ch_reference_fasta_zipped }
            UNZIP_REFERENCE_FASTA( ch_reference_fasta_zipped.collect() )
            ch_reference_fasta = UNZIP_REFERENCE_FASTA.out.unzipped
            ch_versions = ch_versions.mix(UNZIP_REFERENCE_FASTA.out.versions)
        } else {
            Channel.fromPath("${params.reference_fasta}")
                .ifEmpty { error "IMGT DB not found: ${params.reference_fasta}" }
                .set { ch_reference_fasta }
        }

        VALIDATE_IGBLAST_DB(ch_igblast, ch_reference_fasta, ch_airrdb_source)
        ch_igblast = VALIDATE_IGBLAST_DB.out.igblast
        ch_versions = ch_versions.mix(VALIDATE_IGBLAST_DB.out.versions)
    }

    if (params.fetch_database == "imgt" || params.fetch_database == "airrc-imgt") {
        FETCH_DATABASES(Channel.value(params.fetch_database), ch_airrdb_source)
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_reference_fasta = FETCH_DATABASES.out.reference_fasta
        ch_versions = ch_versions.mix(FETCH_DATABASES.out.versions)
    }

    emit:
    versions = ch_versions
    reference_fasta = ch_reference_fasta
    igblast = ch_igblast
}
