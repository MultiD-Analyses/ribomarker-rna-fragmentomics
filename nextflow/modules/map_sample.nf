// nextflow/modules/map_sample.nf
process MAP_SAMPLE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)
    val assets_dir
    val species
    val align_mode

    output:
    tuple val(meta), path("alignments/${meta.id}.merged.bam"), emit: bam
    path "logs/map_${meta.id}.log", emit: log

    script:
    def py_script = "${projectDir}/../docker/scripts/main.py"
    def multimap = params.multimap ?: 100
    def dbpaths = "${assets_dir}/${species}.dbpaths.csv"
    def mode = align_mode ?: 'standard'
    """
    mkdir -p alignments logs

    # Write per-sample sheet for mapReads
    echo "${meta.id},RiboMarker,${reads}" > sample_${meta.id}.csv

    # Execute tiered mapping
    python3 ${py_script} --script mapReads -- \\
        -s sample_${meta.id}.csv \\
        -@ ${task.cpus} \\
        -o . \\
        -d ${dbpaths} \\
        -m ${multimap} \\
        --align-mode ${mode} > logs/map_${meta.id}.log 2>&1
    """
}
