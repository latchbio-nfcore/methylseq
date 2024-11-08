process ARIOC_GENOMEPREPARATION {
    tag "$fasta"
    label 'process_high'

    container "812206152185.dkr.ecr.us-west-2.amazonaws.com/arioc:cuda_12_1_8_V100"

    input:
    path fasta, stageAs: "AriocE_Index/*"

    output:
    path "AriocE_Index" , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python /app/AriocE_Index.py ${fasta} AriocE_Index
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AriocP: 0.0.0
    END_VERSIONS
    """
}
