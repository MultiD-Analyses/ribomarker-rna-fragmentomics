// nextflow/modules/parse_and_count.nf
process PARSE_AND_COUNT {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}/${params.run_id ?: 'ribomarker'}", mode: 'copy'

    input:
    tuple val(meta), path(bams)
    path bedfile

    output:
    path "${meta.id}.bam", emit: final_bam
    path "bedgraph/${meta.id}.bedgraph", emit: bedgraph
    path "counts/${meta.id}-counts.json", emit: json_counts

    script:
    def py_script = "${projectDir}/../docker/scripts/parseMappings.py"
    """
    mkdir -p bedgraph counts alignments

    # Merge tier BAMs into a single coordinate-sorted BAM
    samtools merge -f alignments/${meta.id}.merged.bam ${bams}
    samtools sort -@ ${task.cpus} -o ${meta.id}.bam alignments/${meta.id}.merged.bam
    samtools index ${meta.id}.bam

    # Parse mappings and generate bedgraph coverage
    python3 ${py_script} \\
        --sample ${meta.id} \\
        --bam alignments/${meta.id}.merged.bam \\
        --bed ${bedfile} \\
        --output-bedgraph bedgraph/${meta.id}.bedgraph \\
        --output-counts counts/${meta.id}-counts.json
    """
}
