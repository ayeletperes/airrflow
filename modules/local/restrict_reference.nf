process RESTRICT_REFERENCE {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(reference_base, stageAs: 'source_reference'), path(igblast_base, stageAs: 'source_igblast'), path(ggs_base)
    val database_type

    output:
    tuple val(meta), path("igblast_base"), path("reference_base"), emit: reference
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | head -1 | grep -o "[0-9][0-9.]*"'), emit: versions_igblastn, topic: versions

    script:
    def restriction = meta.locus_restriction ? "-l ${meta.locus_restriction}" : ''
    // --require-loci re-checks coverage here: a zipped set cannot be inspected
    // when the samplesheet is read.
    def ggs = ggs_base ? "-g ${ggs_base} --subject '${meta.ggs_subject}' --require-loci ${meta.required_loci}" : ''
    """
    restrict_reference.py -r source_reference -s "${meta.species}" ${restriction} ${ggs} -o reference_base

    # Rebuild from what survived, so nothing outside the locus and no replaced
    # generic allele is left in a database.
    ref2igblast.sh -i reference_base -o igblast_base -d "${database_type}"

    # ref2igblast.sh makes no support trees; hardlink them rather than duplicate
    # 1.4 MB per reference. MAKE_IGBLAST_AUX replaces the ndm/aux after.
    cp -al source_igblast/internal_data source_igblast/optional_file igblast_base/

    """
}
