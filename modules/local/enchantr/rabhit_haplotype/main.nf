def asString (args) {
    if (args.size() == 0 || args[0] == 'none') return ""
    return args.keySet().sort().collect { param ->
        def value = args[param] instanceof Boolean ? args[param].toString().toUpperCase() : args[param].toString().isNumber() ? args[param].toString() : "'${args[param]}'"
        ",'${param}'=${value}"
    }.join('')
}

process RABHIT_HAPLOTYPE {
    tag "${meta.id}"

    label 'process_medium'
    label 'immcantation'
    label 'immcantation_container'

    container "docker.io/immcantation/airrflow:5.2.2dev"

    input:
    tuple val(meta), path(tabs), path(reference_fasta) // genotyped repertoires of one subject, and their reference
    val anchor                                         // heterozygous anchor gene, e.g. IGHJ6

    output:
    tuple val(meta), path("*_report/haplotypes/*_haplotype.tsv"), emit: haplotype
    path("*_report/haplotypes/*_deletions.tsv"), optional: true, emit: deletions
    path("*/*_command_log.txt"), emit: logs
    path "*_report"
    tuple val("${task.process}"), val('enchantr'), eval('Rscript -e "library(enchantr); cat(as.character(packageVersion(\'enchantr\')))"'), emit: versions_enchantr, topic: versions
    tuple val("${task.process}"), val('rabhit'), eval("Rscript -e \"library(rabhit); cat(as.character(packageVersion('rabhit')))\""), emit: versions_rabhit, topic: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "nf-core/airrflow currently does not support Conda. Please use a container profile instead."
    }
    def args = task.ext.args ? asString(task.ext.args) : ''
    def input = tabs.join(',')
    // Heavy chain only: the report keeps rows with a D call, so light chains in an IG group drop out.
    """
    Rscript -e "enchantr::enchantr_report('rabhit_haplotype', \\
                                        report_params=list('input'='${input}', \\
                                        'imgt_db'='${reference_fasta}', \\
                                        'chain'='IGH', \\
                                        'hapBy'='${anchor}', \\
                                        'outname'='${meta.id}', \\
                                        'outdir'=getwd(), \\
                                        'log'='${meta.id}_haplotype_command_log' ${args}))"

    cp -r enchantr ${meta.id}_haplotype_report && rm -rf enchantr
    """
}
