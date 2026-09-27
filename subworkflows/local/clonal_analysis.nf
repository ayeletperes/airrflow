include { FIND_THRESHOLD as FIND_CLONAL_THRESHOLD } from '../../modules/local/enchantr/find_threshold'
include { FIND_THRESHOLD as REPORT_THRESHOLD } from '../../modules/local/enchantr/find_threshold'
include { CLONAL_ASSIGNMENT } from '../../modules/local/enchantr/clonal_assignment'
include { REPERTOIRE_ANALYSIS} from '../../modules/local/enchantr/repertoire_analysis'
include { DOWSER_LINEAGES } from '../../modules/local/enchantr/dowser_lineages'

workflow CLONAL_ANALYSIS {
    take:
    ch_repertoire_reference
    ch_logo
    clonal_threshold
    clonal_threshold_fallback
    skip_report_threshold
    cloneby
    skip_all_clones_report
    lineage_trees
    genotypeby
    crossby
    singlecell
    lineage_tree_builder
    lineage_tree_exec

    main:
    ch_logs = channel.empty()

    if (clonal_threshold == 'auto') {

        ch_find_threshold = ch_repertoire_reference.map{ it -> it[1] }
                                        .collect()
        ch_find_threshold_samplesheet =  ch_find_threshold
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .collectFile(name: 'find_threshold_samplesheet.txt', newLine: true)

        FIND_CLONAL_THRESHOLD (
            ch_find_threshold,
            ch_logo,
            ch_find_threshold_samplesheet,
            cloneby,
            crossby,
            singlecell
        )
        def ch_threshold = FIND_CLONAL_THRESHOLD.out.mean_threshold

        // Collect raw threshold values into a single list so we can distinguish
        // between (A) no values at all (likely upstream failure), and
        // (B) values present but all invalid thresholds ('' / 'NA' / 'NaN').
        def raw_list = ch_threshold
            .splitText( limit:1 ) { it.trim().toString() }
            .map { it -> it.trim() }
            .collect()

        // Process the collected list to identify when no valid thresholds were found
        clone_threshold = raw_list
            .map { list ->
                def valid = (list ?: []).findAll { it != '' && it != 'NA' && it != 'NaN' }
                if (valid.size() > 0) {
                    return valid
                }
                // Nothing usable from inference: use the stated fallback rather than
                // stopping, and say which value was used.
                if (clonal_threshold_fallback == null || clonal_threshold_fallback == '') {
                    error "Automatic clone_threshold detection produced no usable value. Set --clonal_threshold to a number, or --clonal_threshold_fallback to the value to use when inference fails."
                }
                log.warn "Automatic clone_threshold detection produced no usable value; using --clonal_threshold_fallback ${clonal_threshold_fallback}. Clones in this run are defined at that threshold, not an inferred one."
                return [ clonal_threshold_fallback.toString() ]
            }
            .flatten()

    } else {
        clone_threshold = clonal_threshold

        ch_find_threshold = ch_repertoire_reference.map{ it -> it[1] }
                                        .collect()
        ch_find_threshold_samplesheet =  ch_find_threshold
                        .flatten()
                        .map{ it -> it.getName().toString() }
                        .collectFile(name: 'find_threshold_samplesheet.txt', newLine: true)

        if (!skip_report_threshold){
            REPORT_THRESHOLD (
                ch_find_threshold,
                ch_logo,
                ch_find_threshold_samplesheet,
                cloneby,
                crossby,
                singlecell
            )
        }
    }

    // merge all repertoires by cloneby metadata field
    ch_repertoire_reference.map{ it -> [ it[0][cloneby],
                                it[0].id,
                                it[0].sample_id,
                                it[0].subject_id,
                                it[0].species,
                                it[0].single_cell,
                                it[0].locus,
                                it[1],
                                it[2] ] }
                .groupTuple()
                .map{ get_meta_tabs(it, genotypeby, cloneby) }
                .set{ ch_repertoire_grouped }

    ch_repertoire_grouped.dump(tag: "ch_repertoire_grouped")

    CLONAL_ASSIGNMENT(
        ch_repertoire_grouped,
        clone_threshold.collect(),
        [],
        cloneby,
        singlecell
    )


    // prepare ch for define clones all samples report
    CLONAL_ASSIGNMENT.out.tab
            .map { it -> it[1]}
            .collect()
            .map { it -> [ [id:'all_reps'], it ] }
            .set{ch_all_repertoires_cloned}

    if (!skip_all_clones_report){

        ch_all_repertoires_cloned_samplesheet = ch_all_repertoires_cloned.map{ it -> it[1] }
                                        .collect()
                                        .flatten()
                                        .map{ it -> it.getName().toString() }
                                        .collectFile(name: 'all_repertoires_cloned_samplesheet.txt', newLine: true)

        REPERTOIRE_ANALYSIS(
            ch_all_repertoires_cloned,
            ch_all_repertoires_cloned_samplesheet,
            cloneby
        )
    }

    if (lineage_trees){
        DOWSER_LINEAGES(
            CLONAL_ASSIGNMENT.out.tab,
            lineage_tree_builder,
            lineage_tree_exec
        )
    }

    emit:
    repertoire = REPERTOIRE_ANALYSIS.out.tab
    logs = ch_logs
}

// Function to map
def get_meta_tabs(arr, genotypeby, cloneby) {
    if (arr[3].unique().size() > 1) {
            error "Multiple subject_id found for ${arr[0]} (${arr[3].join(', ')}). Please check your input parameters and ensure that all samples with the same 'cloneby' value have the same 'subject_id' value."
    }

    def meta = [:]
    meta.id                 = [arr[0]].unique().join("")
    // Sorted, so the map is stable as a join key.
    meta.sample_id          = arr[2].flatten().sort()
    meta.subject_id         = arr[3].unique().join("")
    meta.species            = arr[4].unique().join("")
    meta.single_cell        = arr[5].unique().join("")
    meta.locus              = arr[6].unique().join("")

    def array = []

        array = [ meta, arr[7].flatten(), arr[8].unique() ]
        if (arr[8].size() > 1) {
            error "Multiple reference fasta files found for ${meta.id}. Please check your input parameters and ensure that all samples with the same ${genotypeby} value (parameter 'genotype_by') have the same ${cloneby} value (parameter 'clone_by')."
        }
    return array
}
