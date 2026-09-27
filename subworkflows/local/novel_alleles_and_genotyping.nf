include { NOVEL_ALLELE_INFERENCE } from '../../modules/local/enchantr/novel_allele_inference'
include { GENOTYPE_INFERENCE } from '../../modules/local/enchantr/genotype_inference'
include { REASSIGN_ALLELES as REASSIGN_ALLELES_NOVEL; REASSIGN_ALLELES as REASSIGN_ALLELES_GENOTYPE} from '../../modules/local/enchantr/reassign_alleles'
include { CLONAL_ANALYSIS } from './clonal_analysis'
include { CLONAL_ASSIGNMENT as CLONAL_ASSIGNMENT_GENOTYPING } from '../../modules/local/enchantr/clonal_assignment'

workflow NOVEL_ALLELES_AND_GENOTYPING {
    take:
    ch_repertoire // channel: [ val(meta), path(tab), path(reference_fasta) ]
    ch_validated_samplesheet
    ch_logo
    genotypeby
    genotype_method
    allele_thresholds_db
    novel_allele_inference
    single_clone_representative
    genotyping_clonal_threshold
    cloneby
    singlecell

    main:
    ch_logs = channel.empty()

    if (allele_thresholds_db) {
        ch_allele_thresholds_db = channel.fromPath(allele_thresholds_db, checkIfExists: true)
    } else {
        ch_allele_thresholds_db = []
    }

    // merge all repertoires by genotypeby metadata field
    ch_repertoire
        .map{ it ->
                def meta = it[0]
                def rep = it[1]
                def ref = it[2]
                def genotypeby_field = genotypeby=="sample_id" ? "id" : genotypeby
                [ meta[genotypeby_field],
                                    meta.id,
                                    meta.sample_id,
                                    meta.subject_id,
                                    meta.species,
                                    meta.single_cell,
                                    meta.locus,
                                    rep,
                                    ref ] }
                    .groupTuple()
                    .map{ get_meta_tabs(it) }
                    .set{ ch_grouped_repertoires }

    // infer novel alleles
    if (novel_allele_inference) {
        NOVEL_ALLELE_INFERENCE (
            ch_grouped_repertoires
        )

        // reassign novel alleles (we can skip this step if no novel alleles were inferred)
        ch_grouped_repertoires
            .join(NOVEL_ALLELE_INFERENCE.out.reference)
            .map { it ->
                def meta = it[0]
                def reps = it[1]
                def new_ref = it[3]
                [ meta, reps, new_ref ]
            }
            .set{ ch_reassign_alleles }

        REASSIGN_ALLELES_NOVEL (
            ch_reassign_alleles,
            ["v"],
            genotypeby //TODO: @ayeletperes check if this is correct
        )

        REASSIGN_ALLELES_NOVEL.out.tab.dump(tag: "reassign alleles novel")

        REASSIGN_ALLELES_NOVEL.out.tab
            .join(NOVEL_ALLELE_INFERENCE.out.reference)
            .set{ ch_repertoire_reference }

    } else {
        ch_repertoire_reference = ch_grouped_repertoires
    }

    if (single_clone_representative) {
        // Fork so the process and the join each get their own copy, and key the join on
        // the id rather than the whole meta.
        ch_repertoire_reference.dump(tag: 'dbg_reference')
        ch_repertoire_reference
            .multiMap { meta, tabs, ref ->
                to_clones: [ meta, tabs, ref ]
                refs:      [ meta.id, ref ]
            }
            .set { ch_genotype_fork }
        ch_genotype_fork.refs.dump(tag: 'dbg_refs')

        CLONAL_ASSIGNMENT_GENOTYPING(
            ch_genotype_fork.to_clones,
            [genotyping_clonal_threshold],
            [],
            cloneby,
            singlecell
        )
        CLONAL_ASSIGNMENT_GENOTYPING.out.tab.dump(tag: 'dbg_clone_out')
        CLONAL_ASSIGNMENT_GENOTYPING.out.tab
            .map { meta, tab -> [ meta.id, meta, tab ] }
            .join(ch_genotype_fork.refs)
            .map { _id, meta, tab, ref -> [ meta, tab, ref ] }
            .set{ ch_for_genotyping }
        ch_for_genotyping.dump(tag: 'dbg_for_genotype')
    } else {
        ch_for_genotyping = ch_repertoire_reference
    }

    // infer genotype. One module, one enchantr report; `genotype_method` selects between
    // bayesian, fraction and allele_based.
    GENOTYPE_INFERENCE (
        ch_for_genotyping,
        genotypeby,
        single_clone_representative,
        genotype_method,
        ch_allele_thresholds_db
    )
    ch_genotype_reference = GENOTYPE_INFERENCE.out.reference

    ch_grouped_repertoires
        .map{ it -> [it[0], it[1]] }
        .join(ch_genotype_reference)
        .set{ ch_for_reassign }

    ch_genotype_reference.dump(tag: "genotype inference out ref")


    // reassign genotypes
    REASSIGN_ALLELES_GENOTYPE (
        ch_for_reassign,
        ["auto"],
        genotypeby
    )

    REASSIGN_ALLELES_GENOTYPE.out.tab.dump(tag: "reassign alleles genotype out tab")

    ch_repertoire_reference = REASSIGN_ALLELES_GENOTYPE.out.tab.join(ch_genotype_reference)
    ch_repertoire_reference.dump(tag: "ch_repertoire_reference_genotyping")


    emit:
    repertoire_reference = ch_repertoire_reference
    logs = ch_logs
}

// Function to map
def get_meta_tabs(arr) {
    if (arr[3].unique().size() > 1) {
        error "Multiple subject IDs found for ${arr[0]} (${arr[3].join(', ')}). It is not possible to perform joint genotyping of samples from different subjects. Please check the 'genotypeby' parameter."
    }

    def meta = [:]
    meta.id            = [arr[0]].unique().join("")
    // Sorted, so the map is stable as a join key.
    meta.sample_id          = arr[2].flatten().sort()
    meta.subject_id         = arr[3].unique().join("")
    meta.species            = arr[4].unique().join("")
    meta.single_cell        = arr[5].unique().join("")
    meta.locus              = arr[6].unique().join("")
    def array = []

    array = [ meta, arr[7].flatten(), arr[8].unique() ]
    if (arr[8].size() > 1) {
        error "Multiple reference fasta files found for ${meta.id}."
    }


    return array
}
