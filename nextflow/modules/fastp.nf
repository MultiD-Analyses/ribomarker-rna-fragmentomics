// nextflow/modules/fastp.nf
process FASTP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trimmed.fq.gz"), emit: reads
    path "${meta.id}.fastp.json", optional: true, emit: json
    path "${meta.id}.fastp.html", optional: true, emit: html

    script:
    def clip_5p   = params.circularization_clip_5prime_nt ?: 1
    def adapter_3p = params.adapter_3prime_seq ?: "TGGAATTCTCGGGTGCCAAGG"
    def adapter_r2 = params.adapter_r2_seq ?: "TGATCGTCGGACT"
    def min_len    = params.min_length ?: 15

    if (meta.single_end) {
        """
        cutadapt -O 1 -j ${task.cpus} -u ${clip_5p} -m ${min_len} -Z -a ${adapter_3p} \\
            -o ${meta.id}.trimmed.fq.gz ${reads} > ${meta.id}.trim.log
        """
    } else {
        """
        cutadapt -O 1 -j ${task.cpus} -u ${clip_5p} -m ${min_len} -Z \\
            -a ${adapter_3p} -A ${adapter_r2} \\
            -o ${meta.id}.R1.trimmed.fq.gz -p ${meta.id}.R2.trimmed.fq.gz \\
            ${reads[0]} ${reads[1]} > ${meta.id}.trim.log

        fastp \\
            -i ${meta.id}.R1.trimmed.fq.gz \\
            -I ${meta.id}.R2.trimmed.fq.gz \\
            --merged_out ${meta.id}.trimmed.fq.gz \\
            --merge \\
            --thread ${task.cpus} \\
            --overlap_len_require 10 \\
            --overlap_diff_limit 2 \\
            --disable_adapter_trimming \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html
        """
    }
}
