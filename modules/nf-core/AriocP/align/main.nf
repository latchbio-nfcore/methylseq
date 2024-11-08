process ARIOC_ALIGN {
    accelerator 4, type : "nvidia-v100"
    disk '1500 GB'

    tag "$meta.id"


    conda "bioconda::bismark=0.24.0"
    container "812206152185.dkr.ecr.us-west-2.amazonaws.com/arioc:cuda_12_1_8_V100"

    input:
    tuple val(meta), path(reads)
    path index
    val vt
    val match_score
    val mismatch_penalty
    val gap_open_penalty
    val gap_extend_penalty
    val seedDepth
    val batchsize
    val max_j

    output:
    tuple val(meta), path("*bam")       , emit: bam
    tuple val(meta), path("*report.txt"), emit: report
    tuple val(meta), path("*fq.gz")     , optional:true, emit: unmapped
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if(task.ext.prefix){
        args += " --prefix ${task.ext.prefix}"
    }
    def fastq = meta.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"


    """
    python /app/AriocP_Align.py ${reads[0]} ${reads[1]} ${index} \
                                --vt ${vt} --match ${match_score} --gap_open ${gap_open_penalty} \
                                --mismatch ${mismatch_penalty} --gap_extend ${gap_extend_penalty} \
                                --batchsize ${batchsize} --seed_depth ${seedDepth} --max_j ${max_j}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AriocP: 0.0.0
    END_VERSIONS

    """
}
