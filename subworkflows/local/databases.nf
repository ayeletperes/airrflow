include { FETCH_DATABASES } from '../../modules/local/fetch_databases'
include { UNZIP_DB as UNZIP_IGBLAST } from '../../modules/local/unzip_db'
include { UNZIP_DB as UNZIP_REFERENCE_FASTA } from '../../modules/local/unzip_db'
include { VALIDATE_IGBLAST_DB } from '../../modules/local/validate_igblast_db'
include { MAKE_IGBLAST_AUX } from '../../modules/local/make_igblast_aux'
include { RESTRICT_REFERENCE } from '../../modules/local/restrict_reference'

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

    //
    // Resolve which germline reference each sample needs.
    //
    // A sample that asks for no locus restriction and has no personal germline
    // set resolves to the generic reference above and nothing is built, so a run
    // using neither feature behaves exactly as it did before. The rest are keyed,
    // and one reference is built per distinct key -- every IGH-restricted sample
    // in the run shares a single restricted database.
    //
    ch_generic = ch_igblast.combine( ch_reference_fasta )

    ch_meta
        .map { meta -> [ germlineKey(meta), meta.species, meta.locus, meta.locus_restriction ] }
        .unique { row -> row[0] }
        .branch { _key, _species, _locus, restriction ->
            build: restriction != null
            generic: true
        }
        .set { ch_key_lanes }

    ch_key_lanes.generic
        .combine( ch_generic )
        .map { key, _species, _locus, _restriction, igblast, reference ->
            [ key, igblast, reference ]
        }
        .set { ch_generic_by_key }

    RESTRICT_REFERENCE(
        ch_key_lanes.build
            .map { key, species, locus, restriction ->
                [ [ id: key.replaceAll(':', '_'), key: key, species: species,
                    locus: locus, locus_restriction: restriction ] ]
            }
            .combine( ch_generic )
            .map { meta, igblast, reference -> [ meta, reference, igblast ] },
        fetch_germlines == 'airrc-imgt' ? 'airrc-imgt' : 'imgt'
    )

    RESTRICT_REFERENCE.out.reference
        .map { meta, igblast, reference -> [ meta.key, igblast, reference ] }
        .set { ch_built_by_key }

    emit:
    reference_fasta = ch_reference_fasta
    igblast = ch_igblast
    // channel: [ val(key), path(igblast_base), path(reference_base) ]
    reference_by_key = ch_generic_by_key.mix( ch_built_by_key )
}


//
// The key identifying which germline reference a sample needs.
//
// A locus-restricted reference depends only on (species, locus) -- two subjects
// both restricted to IGH share one reference, so it is built once. A personal
// germline set is per subject, so it keys on the subject as well.
//
def germlineKey(meta) {
    def locus = meta.locus_restriction ?: meta.locus.toUpperCase()
    return meta.ggs_subject
        ? "ggs:${meta.subject_id}:${locus}"
        : "gen:${meta.species}:${locus}"
}

//
// Attach the germline reference(s) a task needs to each item of a [meta, ...]
// channel. `kinds` selects which to append and in what order, e.g. ['igblast']
// or ['igblast', 'reference_fasta'].
//
// `combine(by: 0)` and not `join`: join is strictly 1:1 on the key, so it would
// silently drop every sample after the first that shares a reference key. The
// item is nested under its key first so this works for any tuple arity --
// [meta, tab], [meta, R1, R2] and [meta] alike.
//
def withGermline(ch_items, ch_by_key, kinds) {
    def index = [ 'igblast': 1, 'reference_fasta': 2 ]

    return ch_items
        .map { it -> [ germlineKey(it[0]), it ] }
        .combine( ch_by_key, by: 0 )
        .map { entry ->
            def item = entry[1]
            def refs = entry[2..-1]
            item + kinds.collect { kind -> refs[ index[kind] - 1 ] }
        }
}

