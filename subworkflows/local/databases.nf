include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'
include { VALIDATE_IGBLAST_DB } from '../../modules/local/validate_igblast_db'
include { MAKE_IGBLAST_AUX } from '../../modules/local/make_igblast_aux'
include { RESTRICT_REFERENCE } from '../../modules/local/restrict_reference'
include { UNZIP_DB as UNZIP_GGS } from '../../modules/local/unzip_db'

workflow DATABASES {

    take:
    fetch_germlines
    reference_igblast
    reference_fasta
    generate_igblast_aux
    ch_meta // channel: [ val(meta) ] one entry per sample in the run

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

    // One reference per distinct key. A sample needing neither a restriction nor
    // a personal set keys to the generic reference and nothing is built.
    ch_generic = ch_igblast.combine( ch_reference_fasta )

    ch_meta
        .map { meta -> [ germlineKey(meta), meta.species, meta.locus, meta.locus_restriction,
                            meta.ggs_path, meta.ggs_subject, meta.required_loci ] }
        .unique { row -> row[0] }
        .branch { _key, _species, _locus, restriction, ggs, _subject, _required ->
            build: restriction != null || ggs != null
            generic: true
        }
        .set { ch_key_lanes }

    ch_key_lanes.generic
        .combine( ch_generic )
        .map { key, _species, _locus, _restriction, _ggs, _subject, _required, igblast, reference ->
            [ key, igblast, reference ]
        }
        .set { ch_generic_by_key }

    ch_build = ch_key_lanes.build
        .map { key, species, locus, restriction, ggs, subject, required ->
            def meta = [ id: key.replaceAll(':', '_'), key: key, species: species,
                            locus: locus, locus_restriction: restriction ]
            if (ggs) {
                meta.ggs_subject = subject
                meta.required_loci = required
            }
            [ meta, ggs ?: [] ]
        }

    // A germline set given as a directory is staged as it stands. A .zip is
    // unpacked first -- Nextflow cannot list a directory over https, so an
    // archive is the only form a remote germline set can take. UNZIP_DB names
    // its output after the archive, which is also how each build finds its own.
    if (params.ggs_input) {
        ch_build
            .branch { _meta, ggs ->
                zipped: ggs.toString().endsWith('.zip')
                ready: true
            }
            .set { ch_ggs }

        UNZIP_GGS( ch_ggs.zipped.map { _meta, ggs -> file(ggs) }.unique() )

        ch_build = ch_ggs.ready
            .map { meta, ggs -> [ meta, ggs ? file(ggs) : [] ] }
            .mix(
                ch_ggs.zipped
                    .map { meta, ggs -> [ file(ggs).simpleName, meta ] }
                    .combine( UNZIP_GGS.out.unzipped.map { dir -> [ dir.name, dir ] }, by: 0 )
                    .map { _name, meta, dir -> [ meta, dir ] }
            )
    }

    RESTRICT_REFERENCE(
        ch_build
            .combine( ch_generic )
            .map { meta, ggs, igblast, reference -> [ meta, reference, igblast, ggs ] },
        fetch_germlines == 'airrc-imgt' ? 'airrc-imgt' : 'imgt'
    )

    RESTRICT_REFERENCE.out.reference
        .map { meta, igblast, reference -> [ meta.key, igblast, reference ] }
        .set { ch_built_by_key }

    emit:
    igblast = ch_igblast
    // channel: [ val(key), path(igblast_base), path(reference_base) ]
    reference_by_key = ch_generic_by_key.mix( ch_built_by_key )
}


// A restricted reference depends only on (species, locus), so two subjects both
// restricted to IGH share one build. A personal set is per subject.
def germlineKey(meta) {
    def locus = meta.locus_restriction ?: meta.locus.toUpperCase()
    return meta.ggs_subject
        ? "ggs:${meta.subject_id}:${locus}"
        : "gen:${meta.species}:${locus}"
}
