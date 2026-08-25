process MAKE_IGBLAST_AUX {
    tag "make_igblast_aux"
    label 'process_low'

    // No conda directive: receptor_utils is pip-only, it is not on bioconda.
    container "community.wave.seqera.io/library/pip_receptor-utils:21ada4567a12185a"

    input:
    path(igblast_dir, stageAs: "input_igblast_base")
    path(reference_fasta_dir)

    output:
    path("igblast_base"), emit: igblast
    tuple val("${task.process}"), val('receptor_utils'), eval('python3 -c \'from importlib.metadata import version; print(version("receptor_utils"))\''), emit: versions_receptor_utils, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    make_igblast_aux.py -i input_igblast_base -r "${reference_fasta_dir}" -o igblast_base ${args}
    """
}
