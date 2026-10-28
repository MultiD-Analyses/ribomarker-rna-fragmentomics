// nextflow/modules/trim_adapters.nf
process TRIM_READS {
    label 'process_medium'
    publishDir "${params.outdir}/${params.run_id ?: 'ribomarker'}", mode: 'copy'

    input:
    path input_sheet
    val fastq_dir

    output:
    path "results/trimmed/*.fq.gz", emit: trimmed_fastqs
    path "results/*-trimstats.csv", emit: trimstats
    path "results/samples.trimmed.csv", emit: trimmed_sheet

    script:
    def py_script = "${projectDir}/../docker/scripts/main.py"
    def run_id = params.run_id ?: 'ribomarker'
    """
    mkdir -p results/trimmed

    python3 ${py_script} --script trimAdapters -- \\
        --samples ${input_sheet} \\
        --threads ${task.cpus} \\
        --output results \\
        --fastqdir ${fastq_dir}

    if [ ! -f "results/${run_id}-trimstats.csv" ]; then
        cp -f results/results-trimstats.csv "results/${run_id}-trimstats.csv" || true
    fi
    """
}
