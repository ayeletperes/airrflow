include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'
include { VALIDATE_IGBLAST_DB } from '../../modules/local/validate_igblast_db'
include { MAKE_IGBLAST_AUX } from '../../modules/local/make_igblast_aux'
include { BUILD_GERMLINE_REFERENCE } from '../../modules/local/build_germline_reference'
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

    if (fetch_germlines == "imgt" || fetch_germlines == "airrc-imgt") {
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

    // One build per reference key: a personal germline set, a locus restriction, or both.
    // A .zip germline set is unpacked first.
    ch_meta
        .filter { meta -> meta.ggs_path || meta.locus_restriction }
        .unique { meta -> germlineKey(meta) }
        .map { meta -> [ [ id: germlineKey(meta).replaceAll(':', '_'), key: germlineKey(meta), subject: meta.subject_id,
                            species: meta.species, locus: meta.locus_restriction, required_loci: meta.ggs_loci ?: [] ],
                            meta.ggs_path ?: '' ] }
        .branch { _meta, ggs ->
            zipped: ggs.endsWith('.zip')
            dir: true
        }
        .set { ch_ggs }

    UNZIP_GGS( ch_ggs.zipped.map { _meta, ggs -> file(ggs) } )

    BUILD_GERMLINE_REFERENCE(
        ch_ggs.dir.map { meta, ggs -> [ meta, ggs ? file(ggs) : [] ] }
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
        .mix( BUILD_GERMLINE_REFERENCE.out.reference.map { meta, igblast, reference -> [ meta.key, igblast, reference ] } )
}

// A subject with a personal germline set gets its own build; a restricted locus keys a
// build of its own, shared by every sample restricted to that locus. The rest share the
// generic reference.
def germlineKey(meta) {
    def key = meta.ggs_path ? "ggs:${meta.subject_id}" : 'generic'
    return (meta.locus_restriction ? "${key}:${meta.locus_restriction}" : key).toString()
}
