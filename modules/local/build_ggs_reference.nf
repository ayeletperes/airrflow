process BUILD_GGS_REFERENCE {
    tag "$meta.id"
    label 'process_low'
    label 'immcantation'

    conda "bioconda::changeo=1.3.4 bioconda::igblast=1.22.0 conda-forge::wget=1.25.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/9632e731611070a8bc090d51443514959a3179ad6a32be77f0dea1b64f0c12c5/data' :
        'community.wave.seqera.io/library/changeo_igblast_wget:192e77f3b68daa50' }"

    input:
    tuple val(meta), path(ggs), path(igblast_base, stageAs: 'source_igblast'), path(reference_base, stageAs: 'source_reference')

    output:
    tuple val(meta), path("igblast_base"), path("reference_base"), emit: reference
    tuple val("${task.process}"), val('igblastn'), eval('igblastn -version | head -1 | grep -o "[0-9][0-9.]*"'), emit: versions_igblastn, topic: versions

    script:
    def species = meta.species.toLowerCase()
    """
    for locus in ${meta.required_loci.join(' ')}; do
        [ -d ${ggs}/\$locus ] || { echo "ERROR: germline set for subject ${meta.subject} has no \$locus directory" >&2; exit 1; }
    done

    # The subject's set replaces the generic V, D and J of every locus it provides.
    mkdir reference_base && cp -rL source_reference/${species} reference_base/
    for dir in ${ggs}/*/; do
        locus=\$(basename \$dir)
        for pair in V:V_gapped_asc D:D_asc J:J_asc; do
            [ -f \$dir/\${pair#*:}.fasta ] || continue
            rm -f reference_base/${species}/vdj*/*_\${locus}\${pair%%:*}.fasta
            cp \$dir/\${pair#*:}.fasta reference_base/${species}/vdj/imgt_${species}_\${locus}\${pair%%:*}.fasta
        done
    done

    ref2igblast.sh -i reference_base -o igblast_base -d imgt
    cp -rL source_igblast/internal_data source_igblast/optional_file igblast_base/
    """
}
