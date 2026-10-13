def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process GENOTYPE_INFERENCE {
    tag "${meta.id}"

    label 'process_long_parallelized'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.1dev"

    input:
    tuple val(meta), path(tabs), path(reference_fasta) // meta, compressed sequence tsv in AIRR format
    val genotypeby
    val single_clone_representative
    val method                  // bayesian or allele_based
    path allele_thresholds_db   // allele_based only, [] otherwise

    output:
    tuple val(meta), path("*_report/references/*/db_genotype"), emit: reference // reference folder
    path("*/*_command_log.txt"), emit: logs //process logs
    path "*_report/*"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    tuple val("${task.process}"), val('tigger'), eval("Rscript -e \"library(tigger); cat(as.character(packageVersion('tigger')))\""), emit: versions_tigger, topic: versions
    tuple val("${task.process}"), val('piglet'), eval("Rscript -e \"library(piglet); cat(as.character(packageVersion('piglet')))\""), emit: versions_piglet, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def input = tabs.join(',')
    def thresholds = allele_thresholds_db ? ",'allele_thresholds_db'='${allele_thresholds_db}'" : ''
    """
    Rscript -e "enchantr::enchantr_report('genotype', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'genotypeby'='${genotypeby}', \\
                                        'method'='${method}', \\
                                        'single_clone_representative'='${single_clone_representative}', \\
                                        'outdir'=getwd(), \\
                                        'log'='${meta.id}_genotype_inference_command_log' ${thresholds} ${args}))"

    mv enchantr ${meta.id}_genotype_inference_report

    """
}
