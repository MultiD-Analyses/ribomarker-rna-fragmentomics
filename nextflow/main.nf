#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
    RiboMarker: small-RNA and fragmentomics analysis pipeline
================================================================================
    Nextflow DSL2 execution of tiered mapping, feature counting,
    gene collapsing, and coverage generation.
================================================================================
*/

params.help = false
params.input = null
params.outdir = 'results'
params.assets_dir = null
params.species = 'human'
params.run_id = 'ribomarker'
params.multimap = 100
params.adapter_r2_seq = "TGATCGTCGGACT"
params.align_mode = 'standard'
params.skip_qc = false
params.paired_end = null

include { TRIM_READS }     from './modules/trim_adapters'
include { MAP_SAMPLE }     from './modules/map_sample'
include { QUANTIFICATION } from './modules/quantification'
include { QC_REPORT }      from './modules/qc_report'

def helpMessage() {
    log.info"""
    ============================================================================
    RiboMarker Nextflow Pipeline (Release 0.3.0 Compatible)
    ============================================================================
    Usage:
      nextflow run nextflow/main.nf --input <samples.csv> --outdir <results> --assets_dir <references> [-profile <profiles>]

    Parameters:
      --input [file]           Path to sample sheet CSV (columns: sample, fastq_file, fastq_file2)
      --outdir [dir]           Path to output directory (default: results)
      --assets_dir [dir]       Path to references directory
      --species [str]          Species: human or mouse (default: human)
      --run_id [str]           Run identifier folder (default: ribomarker)
      --multimap [int]         Maximum alignments reported per read (default: 100)
      --align_mode [str]       Alignment mode: standard (single-pass) or iterative (fallback for degraded cfRNA)
      --paired_end [bool]      Override library layout (true for PE, false for SE; default: inferred)
      --skip_qc [bool]         Skip QC metrics collection and MultiQC report generation (default: false)

    HPC & Container Options:
      -profile [str]           Execution profiles (comma-separated):
                                 Executors  : local, slurm, pbs, pbspro, sge, lsf
                                 Containers : docker, singularity, apptainer
      --cluster_queue [str]    Scheduler queue / partition name (optional)
      --cluster_options [str]  Additional scheduler flags (e.g. '--account=myproj --qos=normal')
      --max_cpus [int]         Maximum CPUs per task (default: 32)
      --max_memory [str]       Maximum memory per task (default: 128.GB)
      --max_time [str]         Maximum walltime per task (default: 48.h)
      --container_image [str]  Container image URI (default: dockermultid/ribomarker-rna-fragmentomics:0.3.0)
    ============================================================================
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

workflow {
    if (!params.input) {
        error "Input sample sheet not specified! Use --input <samples.csv>"
    }

    def input_file = file(params.input, checkIfExists: true)
    def assets_dir = file(params.assets_dir ?: "${projectDir}/../docker/assets", checkIfExists: true)
    def species = params.species ?: 'human'
    def bedfile = file("${assets_dir}/${species}.bed", checkIfExists: true)
    def dbpaths_file = file("${assets_dir}/${species}.dbpaths.csv", checkIfExists: true)

    log.info """
    ============================================================================
    RiboMarker Nextflow Pipeline Execution
    ============================================================================
    Input Sheet       : ${params.input}
    Output Directory  : ${params.outdir}
    References        : ${assets_dir}
    Species           : ${species}
    Run ID            : ${params.run_id ?: 'ribomarker'}
    Align Mode        : ${params.align_mode ?: 'standard'}
    Profile           : ${workflow.profile}
    Container Engine  : ${workflow.containerEngine ?: 'none (host)'}
    ============================================================================
    """

    // 1. Adapter trimming and pair merging
    def fastq_dir = input_file.parent.toString()
    TRIM_READS(input_file, fastq_dir)

    // 2. Channel mapping for trimmed reads (prefer merged for PE, trimmed for SE)
    ch_trimmed_reads = TRIM_READS.out.trimmed_fastqs.flatten()
        .filter { fq ->
            fq.name.endsWith('.merged.fq.gz') || (!fq.name.contains('.R1.') && !fq.name.contains('.R2.') && fq.name.endsWith('.trimmed.fq.gz'))
        }
        .map { fq ->
            def sample_id = fq.name.replaceAll(/(\.merged|\.trimmed)?\.fq\.gz$/, '')
            def meta = [id: sample_id]
            tuple(meta, fq)
        }

    // 3. Parallel tiered mapping across reference tiers
    def align_mode = params.align_mode ?: 'standard'
    MAP_SAMPLE(ch_trimmed_reads, assets_dir.toString(), species, align_mode)

    // 4. Feature parsing, quantification, and BedGraph coverage generation
    ch_all_bams = MAP_SAMPLE.out.bam.map { meta, bam -> bam }.collect()
    QUANTIFICATION(ch_all_bams, input_file, bedfile, TRIM_READS.out.trimstats, dbpaths_file, assets_dir.toString(), species)

    // 5. Technical QC collector and MultiQC report generation
    if (!params.skip_qc) {
        def first_line = input_file.readLines()[0].trim()
        def is_pe = first_line.contains('fastq2') || first_line.contains('fastq_file2') || (first_line.split(',').size() >= 3 && !first_line.contains('sample,fastq_file'))
        def paired_end = params.paired_end != null ? params.paired_end : is_pe

        QC_REPORT(
            QUANTIFICATION.out.results_dir,
            QUANTIFICATION.out.tables_dir,
            input_file,
            TRIM_READS.out.trimmed_sheet,
            TRIM_READS.out.trim_logs_dir,
            MAP_SAMPLE.out.align_logs.collect(),
            MAP_SAMPLE.out.map_log.collect(),
            QUANTIFICATION.out.metrics_dir,
            QUANTIFICATION.out.ribospike_dir,
            assets_dir.toString(),
            species,
            paired_end
        )
    }
}
