// nextflow/modules/quantification.nf
process QUANTIFICATION {
    label 'process_high'
    publishDir "${params.outdir}/${params.run_id ?: 'ribomarker'}", mode: 'copy'

    input:
    path merged_bams
    path sample_sheet
    path bedfile
    path trimstats
    path dbpaths_file
    val assets_dir
    val species

    output:
    path "results/*-genecounts.csv", emit: genecounts
    path "results/*-readcounts.csv", emit: readcounts
    path "results/*-classcounts.csv", emit: classcounts
    path "results/*-classlengths.csv", optional: true, emit: classlengths
    path "results/*-mappingstats.csv", emit: mappingstats
    path "results/*-featureclasses.csv", emit: featureclasses
    path "results/*-ribospike.csv", optional: true, emit: ribospike
    path "results/*-spikelengths.csv", optional: true, emit: spikelengths
    path "results/*-readorigin.csv", optional: true, emit: readorigin
    path "results/*-trimstats.csv", optional: true, emit: trimstats
    path "tables/*", optional: true, emit: tables
    path "bedgraph/*.bedgraph", emit: bedgraphs
    path "alignments/*.bam", emit: parsed_bams

    script:
    def py_script = "${projectDir}/../docker/scripts/main.py"
    def run_id = params.run_id ?: 'ribomarker'
    """
    mkdir -p alignments results/logs bedgraph ribospike tables

    # Copy input BAMs into alignments/
    for b in ${merged_bams}; do
        cp -f "\$b" alignments/
    done

    # Preserve trimstats
    for t in ${trimstats}; do
        cp -f "\$t" "results/${run_id}-trimstats.csv"
        break
    done

    # Create 3-column downstream sample sheet: sample,group,fastq
    awk -F, 'NF>=1 && \$1!="sample" {print \$1 ",RiboMarker," \$1 ".merged.bam"}' ${sample_sheet} > samples.downstream.csv

    # 0. Generate contig-to-tier map and resolve RNA classes table
    python3 -c "
import sys, os, shutil
sys.path.insert(0, '${projectDir}/../docker/scripts')
import dbPaths, refTables

# Generate contig_tiers.tsv
tiers = dbPaths.read_dbpaths('${dbpaths_file}', '${assets_dir}')
dbPaths.write_contig_tiers(tiers, '${species}', 'tables/contig_tiers.tsv')

# Find and copy RNA classes
rna_table = refTables.find_table('${species}', '${assets_dir}')
shutil.copyfile(rna_table, 'tables/rna_classes.tsv')

# Copy dbpaths for provenance
shutil.copyfile('${dbpaths_file}', 'tables/tiers.csv')
"

    # Determine symbols flag if gene_symbols table exists
    SYMBOLS_ARG=""
    if [ -f "${projectDir}/../docker/tables/${species}.gene_symbols.tsv" ]; then
        SYMBOLS_ARG="--gene-symbols ${projectDir}/../docker/tables/${species}.gene_symbols.tsv"
    fi

    # 1. Parse mappings and generate read counts for feature counting
    python3 ${py_script} --script parseMappings -- \\
        --samples samples.downstream.csv \\
        --threads ${task.cpus} \\
        --output results \\
        --bed ${bedfile} \\
        --rna-classes tables/rna_classes.tsv \\
        --contig-tiers tables/contig_tiers.tsv \\
        \${SYMBOLS_ARG} \\
        --bampath alignments \\
        -q \\
        --remove-multimap \\
        --remove-antisense

    # Rename output prefix if needed to match run_id
    for f in results/results-*.csv; do
        [ -e "\$f" ] || continue
        base=\$(basename "\$f" | sed "s/^results-/${run_id}-/")
        mv "\$f" "results/\$base"
    done

    # 2. Count RiboSpike reads for additional metrics
    python3 ${py_script} --script countRiboSpike -- \\
        -s samples.downstream.csv \\
        -o results \\
        --threads ${task.cpus} \\
        --bampath alignments || true

    for prefix in results ${run_id}; do
        if [ -f "results/ribospike/\${prefix}-ribospike.csv" ]; then
            cp -f "results/ribospike/\${prefix}-ribospike.csv" "results/${run_id}-ribospike.csv"
        fi
        if [ -f "results/ribospike/\${prefix}-spikelengths.csv" ]; then
            cp -f "results/ribospike/\${prefix}-spikelengths.csv" "results/${run_id}-spikelengths.csv"
        fi
        if [ -f "results/ribospike/\${prefix}-readorigin.csv" ]; then
            cp -f "results/ribospike/\${prefix}-readorigin.csv" "results/${run_id}-readorigin.csv"
        fi
    done

    # 3. Collapse per-feature counts to a gene-level, mature-exonic matrix
    python3 ${py_script} --script collapseGenes -- \\
        --readcounts "results/${run_id}-readcounts.csv" \\
        --featureclasses "results/${run_id}-featureclasses.csv" \\
        --rna-classes tables/rna_classes.tsv \\
        --output "results/${run_id}"

    # 4. Make BedGraph for each parsed sample BAM
    for b in alignments/*.bam; do
        [ -e "\$b" ] || continue
        case "\$b" in
            *.merged.bam) continue ;;
        esac
        sample=\$(basename "\$b" .bam)
        python3 ${py_script} --script makeBedgraph -- \\
            --input "\$b" \\
            --output "bedgraph/\${sample}.bedgraph" \\
            --featureclasses "results/${run_id}-featureclasses.csv" \\
            --keep-classes all
    done
    """
}
