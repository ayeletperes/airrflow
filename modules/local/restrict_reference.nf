process RESTRICT_REFERENCE {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(reference_base, stageAs: 'source_reference'), path(igblast_base, stageAs: 'source_igblast')
    val database_type

    output:
    tuple val(meta), path("igblast_base"), path("reference_base"), emit: reference
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | head -1 | grep -o "[0-9][0-9.]*"'), emit: versions_igblastn, topic: versions

    script:
    """
    restrict_reference.py -r source_reference -l "${meta.locus_restriction}" -s "${meta.species}" -o reference_base

    # Rebuild the canonical FASTAs and BLAST databases from what survived, so no
    # gene outside the locus is left in a database.
    ref2igblast.sh -i reference_base -o igblast_base -d "${database_type}"

    # ref2igblast.sh does not produce the NCBI support trees; carry them over as
    # hardlinks so each restricted reference does not duplicate 1.4 MB of them.
    # MAKE_IGBLAST_AUX replaces the germline-derived ndm/aux afterwards.
    cp -al source_igblast/internal_data source_igblast/optional_file igblast_base/

    """
}
