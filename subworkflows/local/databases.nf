include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'
include { VALIDATE_IGBLAST_DB } from '../../modules/local/validate_igblast_db'
include { MAKE_IGBLAST_AUX } from '../../modules/local/make_igblast_aux'
include { BUILD_GGS_REFERENCE } from '../../modules/local/build_ggs_reference'
include { UNZIP_DB as UNZIP_GGS } from '../../modules/local/unzip_db'

workflow DATABASES {

    take:
    fetch_germlines
    reference_igblast
    reference_fasta
    generate_igblast_aux
    ch_meta // channel: [ val(meta) ] one per sample

    main:

    if( fetch_germlines == "none"){
        if (reference_igblast.endsWith(".zip")) {
            channel.fromPath("${reference_igblast}")
                    .ifEmpty{ error "IGBLAST DB not found: ${reference_igblast}" }
                    .set { ch_igblast_zipped }
            UNZIP_IGBLAST( ch_igblast_zipped.collect() )
            ch_igblast = UNZIP_IGBLAST.out.unzipped
        } else {
            channel.fromPath("${reference_igblast}")
                .ifEmpty { error "IGBLAST DB not found: ${reference_igblast}" }
                .set { ch_igblast }
        }
    }

    if( fetch_germlines == "none" ){
        if (reference_fasta.endsWith(".zip")) {
            channel.fromPath("${reference_fasta}")
                    .ifEmpty{ error "IMGTDB not found: ${reference_fasta}" }
                    .set { ch_reference_fasta_zipped }
            UNZIP_REFERENCE_FASTA( ch_reference_fasta_zipped.collect() )
            ch_reference_fasta = UNZIP_REFERENCE_FASTA.out.unzipped
        } else {
            channel.fromPath("${reference_fasta}")
                .ifEmpty { error "IMGT DB not found: ${reference_fasta}" }
                .set { ch_reference_fasta }
        }

        VALIDATE_IGBLAST_DB(ch_igblast, ch_reference_fasta)
        ch_igblast = VALIDATE_IGBLAST_DB.out.igblast
    }

    if (fetch_germlines in ["imgt", "airrc-imgt", "ogrdb"]) {
        FETCH_DATABASES(channel.value(fetch_germlines))
        ch_igblast = FETCH_DATABASES.out.igblast
        ch_reference_fasta = FETCH_DATABASES.out.reference_fasta
    }

    // Replace the stock NCBI .ndm / .aux files with ones derived from the
    // germline reference actually in use. The boundaries and reading frames they
    // encode are germline-set specific, so the shipped files mis-annotate a
    // custom reference.
    if (generate_igblast_aux) {
        MAKE_IGBLAST_AUX(ch_igblast, ch_reference_fasta)
        ch_igblast = MAKE_IGBLAST_AUX.out.igblast
    }

    ch_generic = ch_igblast.combine( ch_reference_fasta )

    // One build per subject with a personal germline set; a .zip is unpacked first.
    ch_meta
        .filter { meta -> meta.ggs_path }
        .unique { meta -> meta.subject_id }
        .map { meta -> [ [ id: "ggs_${meta.subject_id}", subject: meta.subject_id, species: meta.species, required_loci: meta.ggs_loci ], meta.ggs_path ] }
        .branch { _meta, ggs ->
            zipped: ggs.endsWith('.zip')
            dir: true
        }
        .set { ch_ggs }

    UNZIP_GGS( ch_ggs.zipped.map { _meta, ggs -> file(ggs) } )

    BUILD_GGS_REFERENCE(
        ch_ggs.dir.map { meta, ggs -> [ meta, file(ggs) ] }
            .mix( ch_ggs.zipped
                .map { meta, ggs -> [ file(ggs).simpleName, meta ] }
                .combine( UNZIP_GGS.out.unzipped.map { dir -> [ dir.name, dir ] }, by: 0 )
                .map { _name, meta, dir -> [ meta, dir ] } )
            .combine( ch_generic )
    )

    emit:
    igblast = ch_igblast
    // channel: [ val(key), path(igblast_base), path(reference_base) ]
    reference_by_key = ch_generic.map { igblast, reference -> [ 'generic', igblast, reference ] }
        .mix( BUILD_GGS_REFERENCE.out.reference.map { meta, igblast, reference -> [ "ggs:${meta.subject}".toString(), igblast, reference ] } )
}

// Samples of a subject with a personal germline set use its build; the rest share the generic reference.
def germlineKey(meta) {
    return meta.ggs_path ? "ggs:${meta.subject_id}".toString() : 'generic'
}
