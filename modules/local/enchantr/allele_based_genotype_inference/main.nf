def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process ALLELE_BASED_GENOTYPE_INFERENCE {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'

    // ponytail: this pin is the reason this module is separate from
    // BAYESIAN_GENOTYPE_INFERENCE. Bump it here when the enchantr image that knows
    // 'piglet_genotype' is built; the bayesian module keeps 5.1.0.
    container "docker.io/immcantation/airrflow:5.1.0"

    input:
    tuple val(meta), path(tabs), path(reference_fasta) // meta, sequence tsv in AIRR format
    val single_clone_representative
    path allele_thresholds_db

    output:
    // ponytail: this is the shape REASSIGN_ALLELES_GENOTYPE and the ch_repertoire_reference join
    // require, not a shape the piglet report currently produces: its index.Rmd writes a flat
    // `<outname>_personal_reference.fasta` into outdir and never calls generate_genotyped_reference,
    // so there is no references/<id>/db_genotype directory. Left as the downstream contract rather
    // than guessed at; the report has to emit this, or the join needs rewriting.
    tuple val(meta), path("*_report/references/*/db_genotype"), emit: reference // reference folder
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    tuple val("${task.process}"), val('piglet'), eval("Rscript -e \"library(piglet); cat(as.character(packageVersion('piglet')))\""), emit: versions_piglet, topic: versions


    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def input = tabs.join(',')
    // ponytail: piglet_genotype_project_files/index.Rmd names the germline `germline`, not
    // `imgt_db`, and declares no `genotypeby`, so this is not the bayesian parameter list.
    // Its remaining parameters (find_unmutated, single_assignments, default_threshold,
    // translate_to_asc, asc_annotation, reassign, call/seq/germline/clone_id columns) have no
    // airrflow option to map onto and are left at the template defaults.
    """
    Rscript -e "enchantr::enchantr_report('piglet_genotype', \\
                                        report_params=list('input'='${input}', \\
                                        'germline'='${reference_fasta}', \\
                                        'allele_thresholds_db'='${allele_thresholds_db}', \\
                                        'single_clone_representative'='${single_clone_representative}', \\
                                        'outdir'=getwd(), \\
                                        'log'='${meta.id}_allele_based_genotype_inference_command_log' ${args}))"

    cp -r enchantr ${meta.id}_allele_based_genotype_inference_report && rm -rf enchantr

    """
}
