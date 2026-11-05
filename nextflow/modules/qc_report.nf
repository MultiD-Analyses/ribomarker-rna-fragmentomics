// nextflow/modules/qc_report.nf
process QC_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/${params.run_id ?: 'ribomarker'}", mode: 'copy'

    input:
    path results_dir
    path tables_dir
    path sample_sheet
    path trimmed_sheet
    path trim_logs_dir
    path align_logs
    path map_logs
    path metrics_dir
    path ribospike_dir
    val assets_dir
    val species
    val paired_end

    output:
    path "qc/report_technical.html", emit: html_report
    path "qc/qc_metrics.json", emit: qc_metrics
    path "qc/metrics_long.tsv", emit: metrics_long
    path "qc/consistency_checks.tsv", emit: consistency_checks
    path "qc/run_config.yaml", emit: run_config
    path "qc/mqc/*", optional: true, emit: mqc_files
    path "qc", emit: qc_dir

    script:
    def py_script = "${projectDir}/../docker/scripts/main.py"
    def run_id = params.run_id ?: 'ribomarker'
    def pe_bool = paired_end ? 'true' : 'false'
    def multimap = params.multimap ?: 100
    """
    # 1. Structure workspace to match qcCollect & qcReport requirements
    mkdir -p qc mqc tables trimmed/logs alignments/logs logs metrics ribospike/results results

    # Copy tables if needed
    if [ "${tables_dir}" != "tables" ]; then
        cp -rf ${tables_dir}/* tables/ 2>/dev/null || true
    fi

    # Copy trimming logs
    if [ "${trim_logs_dir}" != "trimmed/logs" ]; then
        cp -rf ${trim_logs_dir}/* trimmed/logs/ 2>/dev/null || true
    fi

    # Copy alignment tier logs
    for f in ${align_logs}; do
        [ "\$f" = "alignments/logs/\$(basename \$f)" ] || cp -f "\$f" alignments/logs/
    done

    # Consolidate map.log
    for f in ${map_logs}; do
        cp -f "\$f" logs/map.log
        break
    done

    # Copy sample metrics logs
    if [ "${metrics_dir}" != "metrics" ]; then
        cp -rf ${metrics_dir}/* metrics/ 2>/dev/null || true
    fi

    # Copy count tables into results/ and current directory
    if [ "${results_dir}" != "results" ]; then
        cp -rf ${results_dir}/* results/ 2>/dev/null || true
    fi
    for f in results/*.csv; do
        [ -e "\$f" ] || continue
        base=\$(basename "\$f")
        [ -e "\$base" ] || cp -f "\$f" "\$base"
    done

    # Copy ribospike files
    if [ "${ribospike_dir}" != "ribospike" ]; then
        cp -rf ${ribospike_dir}/* ribospike/ 2>/dev/null || true
    fi
    for f in ribospike/*.csv; do
        [ -e "\$f" ] || continue
        base=\$(basename "\$f")
        [ -e "ribospike/results/\$base" ] || cp -f "\$f" "ribospike/results/\$base"
    done

    # Setup sample sheets
    [ -f "sample_sheet.csv" ] || cp -f "${sample_sheet}" sample_sheet.csv
    [ -f "samples.csv" ] || cp -f "${sample_sheet}" samples.csv
    [ -f "samples.trimmed.csv" ] || cp -f "${trimmed_sheet}" samples.trimmed.csv

    # Create self-referencing symlink so ${run_id}/results/ resolves to results/
    [ -e "${run_id}" ] || ln -s . "${run_id}"

    # 2. Write qc/run_config.yaml
    cat <<EOF > qc/run_config.yaml
dbdir: "${assets_dir}/"
input_dir: "${run_id}"
output_prefix: "${run_id}"
samples: "samples.csv"
sample_sheet: "sample_sheet.csv"
species: "${species}"
threads: ${task.cpus}
trim_fastq: true
trim_only: false
paired_end: ${pe_bool}
count_ribospike: true
qc_report: true
keep_classes: "all"
pos_threshold: 10
remove_multimap: true
remove_antisense: true
collapse: false
multimap: ${multimap}
subsample_reads: 0
rna_classes: "tables/rna_classes.tsv"
tiers: "tables/tiers.csv"
gene_symbols: "${projectDir}/../docker/tables/${species}.gene_symbols.tsv"
report_title: "${run_id}"
EOF

    # 3. Execute QC collector
    python3 ${py_script} --script qcCollect -- \\
        --output "${run_id}" \\
        --config qc/run_config.yaml \\
        --qc-dir "${projectDir}/../docker/qc" \\
        --version-file "${projectDir}/../docker/VERSION"

    # 4. Execute MultiQC technical report renderer
    python3 ${py_script} --script qcReport -- \\
        --output "${run_id}" \\
        --qc-dir "${projectDir}/../docker/qc" \\
        --title "${run_id}"
    """
}
