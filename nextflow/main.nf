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

include { TRIM_READS }     from './modules/trim_adapters'
include { MAP_SAMPLE }     from './modules/map_sample'
include { QUANTIFICATION } from './modules/quantification'

def helpMessage() {
    log.info"""
    ============================================================================
    RiboMarker Nextflow Pipeline (Release 0.3.0 Compatible)
    ============================================================================
    Usage:
      nextflow run nextflow/main.nf --input <samples.csv> --outdir <results> --assets_dir <references>

    Parameters:
      --input [file]           Path to sample sheet CSV (columns: sample, fastq_file, fastq_file2)
      --outdir [dir]           Path to output directory (default: results)
      --assets_dir [dir]       Path to references directory
      --species [str]          Species: human or mouse (default: human)
      --run_id [str]           Run identifier folder (default: ribomarker)
      --multimap [int]         Maximum alignments reported per read (default: 100)
      --align_mode [str]       Alignment mode: standard (single-pass) or iterative (fallback for degraded cfRNA)
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
}
