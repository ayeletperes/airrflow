include { AMULETY_TRANSLATE  } from '../../modules/nf-core/amulety/translate/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ANTIBERTY} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ANTIBERTA2} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ESM2} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_BALMPAIRED} from '../../modules/nf-core/amulety/embed/main'
include { germlineKey } from './databases'

workflow TRANSLATE_EMBED {
    take:
    ch_repertoire
    ch_reference_by_key // channel: [ val(key), path(igblast_base), path(reference_base) ]
    embeddings
    embedding_chain

    main:

    // AMULETY_TRANSLATE is an nf-core module and takes the reference as a separate
    // broadcast input, so split the per-sample reference back out of the tuple.
    ch_repertoire
        .map { meta, tab -> [ germlineKey(meta), meta, tab ] }
        .combine( ch_reference_by_key, by: 0 )
        .multiMap { _key, meta, tab, igblast, _reference ->
            repertoire: [ meta, tab ]
            igblast: igblast
        }
        .set { ch_translate }

    AMULETY_TRANSLATE(
        ch_translate.repertoire,
        ch_translate.igblast
    )

    if (embeddings && embeddings.split(',').contains('antiberty') ){
        AMULETY_EMBED_ANTIBERTY(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "antiberty"
        )
    }

    if (embeddings && embeddings.split(',').contains('antiberta2') ){
        AMULETY_EMBED_ANTIBERTA2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "antiberta2"
        )
    }

    if (embeddings && embeddings.split(',').contains('esm2') ){
        AMULETY_EMBED_ESM2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "esm2"
        )
    }

    if (embeddings && embeddings.split(',').contains('balmpaired') ){
        AMULETY_EMBED_BALMPAIRED(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "balm-paired"
        )
    }

}
