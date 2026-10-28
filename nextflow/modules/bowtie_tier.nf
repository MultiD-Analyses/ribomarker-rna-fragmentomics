// nextflow/modules/bowtie_tier.nf
process BOWTIE_TIER {
    tag "$meta.id - $tier_name"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    val tier_name
    path index_files
    val index_prefix
    val aligner

    output:
    tuple val(meta), path("${meta.id}.${tier_name}.bam"), emit: bam
    tuple val(meta), path("${meta.id}.${tier_name}.unmapped.fq.gz"), optional: true, emit: unmapped
    path "${meta.id}.${tier_name}.log", emit: log

    script:
    if (aligner == 'bowtie') {
        """
        bowtie \\
            -p ${task.cpus} \\
            --norc \\
            -v 0 \\
            -S \\
            ${index_prefix} \\
            ${reads} \\
            --un ${meta.id}.${tier_name}.unmapped.fq \\
            2> ${meta.id}.${tier_name}.log | \\
            samtools view -bS - > ${meta.id}.${tier_name}.bam

        gzip -f ${meta.id}.${tier_name}.unmapped.fq || true
        """
    } else if (aligner == 'ribospikebt2') {
        """
        bowtie2 \\
            --sam-no-qname-trunc \\
            --no-unal \\
            --no-1mm-upfront \\
            --mp 6 \\
            --rdg 6,3 \\
            --rfg 6,3 \\
            --norc \\
            --n-ceil L,5,0.15 \\
            -k 1 \\
            -D 20 \\
            -R 2 \\
            -N 0 \\
            -L 10 \\
            -i S,1,0.40 \\
            --score-min L,0,-0.6 \\
            --np 1 \\
            --ignore-quals \\
            --threads ${task.cpus} \\
            --un-gz ${meta.id}.${tier_name}.unmapped.fq.gz \\
            -x ${index_prefix} \\
            -U ${reads} \\
            2> ${meta.id}.${tier_name}.log | \\
            samtools view -bS - > ${meta.id}.${tier_name}.bam
        """
    } else {
        """
        bowtie2 \\
            --sam-no-qname-trunc \\
            --no-unal \\
            --np 5 \\
            --very-sensitive \\
            --ignore-quals \\
            --threads ${task.cpus} \\
            -k ${params.multimap ?: 100} \\
            --un-gz ${meta.id}.${tier_name}.unmapped.fq.gz \\
            -x ${index_prefix} \\
            -U ${reads} \\
            2> ${meta.id}.${tier_name}.log | \\
            samtools view -bS - > ${meta.id}.${tier_name}.bam
        """
    }
}
