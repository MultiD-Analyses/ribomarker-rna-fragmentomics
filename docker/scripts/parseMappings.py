#!/usr/bin/env python3
#Aidan Manning 7.9.24
#Script to map small RNA reads to the genome, ncRNAs, and other RNAs using various mapping tools
import argparse
import hashlib
import sys
import pandas as pd
import os
from datetime import datetime
from multiprocessing import Pool
from toolRun import run  # stops on a failed tool
import refTables  # the RNA class table (release 0.3.0)
import pysam
from collections import defaultdict, deque, Counter
import tqdm
import numpy as np
import re
from numba import jit

FEATURE_SUFFIXES = ('_multimap', '_antisense', '_5p', '_3p', '_internal', '_whole')


def base_feature(name, known):
    """The BED feature a counted feature belongs to: parseMappings suffixes feature names (_5p/_3p/_internal/_whole,
    _multimap, _antisense). A suffix is removed only while the name is not a BED feature itself: MirGeneDB arms are
    named ..._5p / ..._3p in the BED and keep that suffix (before 0.3.0 it was stripped, so 286 of 561 miRNA features
    of sample S3 had no length; release 0.3.0). known: the BED feature names."""
    def strip(name):
        changed = True
        while changed and name not in known:
            changed = False
            for suffix in FEATURE_SUFFIXES:
                if name.endswith(suffix):
                    name = name[:-len(suffix)]
                    changed = True
                    break
        return name
    # MirGeneDB arms with identical sequences are one BED feature named 'A_5p|B_5p|...': the whole name first,
    # then (as before 0.3.0) the part after the last '|'
    whole = strip(name)
    return whole if whole in known else strip(name.split('|')[-1])


def add_gene_length(featuretypes, bed_path):
    """Attach a length_nt column (gene_length before 0.3.0): the length of the base feature in the BED."""
    bed = pd.read_csv(bed_path, sep='\t', header=None,
                      usecols=[1, 2, 3], names=['start', 'end', 'feature'])
    base_length = (bed['end'] - bed['start']).groupby(bed['feature']).max().to_dict()
    featuretypes = featuretypes.copy()
    featuretypes['length_nt'] = featuretypes['feature'].map(lambda f: base_length.get(base_feature(f, base_length)))
    return featuretypes


def read_gene_symbols(path):
    """A gene-symbol table (buildDB/make_gene_symbols.py): (value columns, {feature: values}, BED sha256)."""
    bed_sha256, lines = None, []
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('# references_bed_sha256:'):
                bed_sha256 = line.split(':', 1)[1].strip()
            elif not line.startswith('#') and line.strip():
                lines.append(line.rstrip('\n').split('\t'))
    header, rows = lines[0], lines[1:]
    return header[1:], {row[0]: row[1:] + [''] * (len(header) - len(row)) for row in rows}, bed_sha256


def add_gene_symbols(featuretypes, bed_path, symbols_path):
    """Attach the gene-symbol columns (release 0.3.0) by base feature. The table is used only if it was built
    for this BED (sha256 in its header); otherwise the columns stay empty and a warning is printed."""
    columns, symbols, bed_sha256 = read_gene_symbols(symbols_path)
    with open(bed_path, 'rb') as f:
        matches = hashlib.sha256(f.read()).hexdigest() == bed_sha256
    if not matches:
        print(f'WARNING: {symbols_path} was built for another BED than {bed_path}; gene symbol columns left empty',
              file=sys.stderr)
        symbols = {}
    known = set(pd.read_csv(bed_path, sep='\t', header=None, usecols=[3], names=['feature'])['feature'])
    featuretypes = featuretypes.copy()
    bases = featuretypes['feature'].map(lambda f: base_feature(f, known))
    for number, column in enumerate(columns):
        featuretypes[column] = bases.map(lambda base: symbols.get(base, [''] * len(columns))[number])
    return featuretypes


def position_suffix(read_start, read_end, start, end, strand):
    """_5p / _3p / _internal / _whole: where a read lies in its feature (0-based, end exclusive, as the BED and
    pysam give them). Near an end = within 10 nt of it (Clark et al., RNA 32:1156, 2026, Methods); _whole: near
    both ends (not in the paper). The 5' end is the feature's start on the + strand and its end on the - strand
    (before 0.3.0 the strand was ignored, which swapped 5p and 3p of features on the - strand of a reference)."""
    near_start, near_end = read_start <= start + 10, read_end >= end - 10
    near_5p, near_3p = (near_end, near_start) if strand == '-' else (near_start, near_end)
    if near_5p and near_3p:
        return '_whole'
    if near_5p:
        return '_5p'
    if near_3p:
        return '_3p'
    return '_internal'


def read_univec_contigs(path):
    """The contigs of the UniVec tier(s) in a contig-to-tier map (contig<TAB>tier, with a header line)."""
    with open(path) as f:
        next(f, None)
        return {contig for contig, tier in (line.rstrip('\n').split('\t') for line in f if line.strip())
                if tier.startswith('univec')}


def to_classes(class_of, counts, lengths, distinct):
    """The per-label tallies of one sample relabelled to RNA class ids (class table, release 0.3.0); labels of
    one class are merged. counts {label: reads}; lengths {label: {length: reads}}; distinct {label: {sequence:
    reads}}: the union of the sequences, so a distinct sequence is counted once per class."""
    new_counts, new_lengths, new_distinct = {}, {}, {}
    for label, n in counts.items():
        rna_class = class_of(label)
        new_counts[rna_class] = new_counts.get(rna_class, 0) + n
    for label, histogram in lengths.items():
        target = new_lengths.setdefault(class_of(label), {length: 0 for length in histogram})
        for length, n in histogram.items():
            target[length] = target.get(length, 0) + n
    for label, sequences in distinct.items():
        target = new_distinct.setdefault(class_of(label), {})
        for sequence, n in sequences.items():
            target[sequence] = target.get(sequence, 0) + n
    return new_counts, new_lengths, new_distinct

@jit(nopython=True)
def find_region(starts, ends, read_start, read_end):
    for i in range(len(starts)):
        # Check for any overlap between the intervals
        if read_start <= ends[i] and read_end >= starts[i]:
            return i  # Return index of the first overlapping region
    return -1  # Indicate no overlap was found


def sort_by_rank(bed, genetype_hierarchy):
    """The BED rows from the highest-ranked biotype down (genetype_hierarchy: biotype -> rank, unknown 0).
    Rows of equal rank keep their order in the BED file (stable sort), so the same feature wins a tie on every
    computer: numpy's default quicksort orders ties differently on different CPUs (found for release 0.2.0)."""
    order = bed['biotype'].map(genetype_hierarchy).fillna(0)
    return bed.assign(order=order).sort_values(by='order', ascending=False, kind='stable').drop(columns='order')


class ParseMappings(object):

    def __init__(self, args):
        '''
        This class is used to parse the mappings from the various mapping tools.
        '''

        self.args = args
        # ranks of overlapping features and the labels with positional suffixes come from the class table
        self.classes = refTables.ClassTable(args.rna_classes or refTables.default_table('default'), 'default')
        self.position_labels = self.classes.position_labels()
        # reads without a feature on a UniVec contig are the class UniVec (release 0.3.0); the contigs come from
        # the contig-to-tier map of the run (Snakefile rule run_tables, dbPaths.write_contig_tiers)
        self.univec_contigs = read_univec_contigs(args.contig_tiers) if args.contig_tiers else set()

        self.sampleinfo = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq'])
        self.samplelist = self.sampleinfo['sample'].tolist()

        self.logdata = {}
        self.logdata['module'] = 'parseMappings'
        self.logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.logdata['samples'] = self.samplelist
        self.logdata['threads'] = self.args.threads
        self.logdata['output'] = self.args.output
        self.logdata['mode'] = 'bed'
        self.logdata['bed'] = self.args.bed
        self.logdata['tools'] = []
        self.logdata['commands'] = []
        self.logdata['input_bams'] = {}

        return

    def pick_bed_alignment(self, read_mappings_dict, read_name, gD):
        '''
        This function picks the best mapping for a read based on the alignment score, strandedness, and biotype.
        '''
        
        read_data = read_mappings_dict[read_name]
        alignments = read_data['alignments']
        features = read_data['annotations']
        alignmentscores = read_data['alignmentscores']
        biotypes = read_data['biotypes']
        strands = read_data['strand']
        barcodes = read_data['barcodes']
        biotype_scores = [gD[x] for x in biotypes]

        def assign_strandedness_tag(read, strand, barcode, status):
            strandedness = 'sense' if strand in [1, -2] else 'antisense'
            read.set_tag('XF', f"{barcode}|{strandedness}|{status}")
            return strandedness

        if len(alignments) == 1:
            # Handle single alignment case
            read = alignments[0]
            strand = strands[0]
            barcode = barcodes[0]
            strandedness = assign_strandedness_tag(read, strand, barcode, 'Unique')
            return read, 'Unique', barcode, strandedness

        elif len(alignments) > 1:
            # Handle multiple alignment cases
            # Get the index of the highest alignment score
            max_alignscore = max(alignmentscores)
            highest_alignscore_indices = [i for i, score in enumerate(alignmentscores) if score == max_alignscore]
            
            if len(highest_alignscore_indices) == 1:
                # If only one alignment has the highest score, return it
                best_index = highest_alignscore_indices[0]
                mappability = 'Unique'
                read = alignments[best_index]
                strand = strands[best_index]
                barcode = barcodes[best_index]
                strandedness = assign_strandedness_tag(read, strand, barcode, mappability)
                return read, mappability, barcode, strandedness
            else:   # several alignments with the best score
                # Resolve ties deterministically and independently of alignment order:
                # score, then strand, then biotype rank, then genomic position, then barcode.
                def tie_key(i):
                    r = alignments[i]
                    return (alignmentscores[i], strands[i], biotype_scores[i],
                            r.reference_id, r.reference_start, r.reference_end, barcodes[i])
                best_index = max(highest_alignscore_indices, key=tie_key)

                mappability = 'Multi'
                read = alignments[best_index]
                strand = strands[best_index]
                barcode = barcodes[best_index]
                strandedness = assign_strandedness_tag(read, strand, barcode, mappability)
                return read, mappability, barcode, strandedness

    def generate_metrics_bed_distinct(self, read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict):

        if barcode != 'unannotated':
            bed_info = bed_barcode_dict.get(barcode)
            feature = bed_info.get('feature')
            biotype = bed_info.get('biotype')
            start = bed_info.get('start')
            end = bed_info.get('end')
            strand = bed_info.get('strand')
            feature_length = end - start 
            read_start = read.reference_start
            read_end = read.reference_end

            if read.is_reverse:
                query_sequence = read.query_sequence[::-1].translate(str.maketrans('ATCGN', 'TAGCN'))
            else:
                query_sequence = read.query_sequence

            read_length = len(query_sequence)
            
            # an antisense read keeps its class and its feature gets _antisense, for every class (release 0.4.0;
            # before: class antisense, and protein_coding / pre-mRNA not tagged)
            if strandedness == 'antisense':
                feature = feature + '_antisense'

            if biotype in self.position_labels:     # tRNA, Mt_tRNA, snRNA, snoRNA (class table: positions = yes)
                feature = feature + position_suffix(read_start, read_end, start, end, strand)

            # Flag every multimapping read, for all biotypes — so multimappers
            # can be separated or removed reliably downstream.
            if read.get_tag('XF').split('|')[-1] == 'Multi':
                feature = feature + '_multimap'
                    
            if featuretypeD.get(feature) == None:
                featuretypeD[feature] = biotype

            read_biotype_dict[biotype] += count_value
            read_length_dict[biotype][read_length] += count_value
            read_feature_dict[feature] += count_value
            read.set_tag('CS', f"{feature}")

            read_feature_distinctseqs_dict[feature][query_sequence] += count_value

            read_biotypes_distinctseqs_dict[biotype][query_sequence] += count_value

            return read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict

        else:

            # no feature: a vector contaminant (UniVec contig) or a genome read without annotation
            biotype = 'UniVec' if read.reference_name in self.univec_contigs else 'unannotated'
            read_length = len(read.query_sequence)
            query_sequence = read.query_sequence

            read_biotype_dict[biotype] += count_value
            read_length_dict[biotype][read_length] += count_value
            read.set_tag('CS', f"other")

            read_biotypes_distinctseqs_dict[biotype][query_sequence] += count_value

            return read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict

    def parse_mapping_bed_distinct(self, sample):
        '''
        Assigns reads to bed features.
        '''

        bampath = self.args.bampath.rstrip('/') + '/'      # main() sets the default

        try:
            bed = pd.read_csv(self.args.bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'feature', 'biotype', 'strand'])
        except pd.errors.EmptyDataError:
            print(f'{self.args.bed} is empty. Check the file and rerun the script.')
            sys.exit(1)
        except FileNotFoundError:
            print(f'{self.args.bed} does not exist. Check the path and rerun the script.')
            sys.exit(1)

        ## creates a unique identifier for each feature
        bed['barcode'] = bed['feature'].astype(str) + '::' + bed['chrom'].astype(str) + ':' + bed['start'].astype(str) + '-' + bed['end'].astype(str) + ':' + bed['strand'].astype(str)
        ## get the unique biotypes in the bed file
        biotypes = bed['biotype'].drop_duplicates().tolist()

        ## make a barcodes dictionary
        bed_barcode_dict = bed.set_index('barcode').T.to_dict('dict')

        uniquefeatures = bed.drop_duplicates(subset=['feature'], keep='first')
        featuretypeD = dict(zip(uniquefeatures['feature'], uniquefeatures['biotype']))
        del uniquefeatures

        ## rank of each label from the class table (labels missing from it get rank 10, as before 0.3.0)
        genetype_hierarchy = self.classes.rank_map(biotypes + ['unannotated'])

        ## sort the bed file by rank; equal ranks keep the BED order
        bed = sort_by_rank(bed, genetype_hierarchy)

        ## group the bed file by chromosome
        bed = bed.groupby('chrom')

        ## the annotations of each chromosome (both strands, in rank order), for a quick lookup
        bed_dict = {}
        for chrom, group_df in bed:
            bed_dict[chrom] = {
                'start': group_df['start'].to_numpy(),
                'end': group_df['end'].to_numpy(),
                'feature': group_df['feature'].tolist(),
                'biotype': group_df['biotype'].tolist(),
                'strand': group_df['strand'].tolist(),
                'barcode': group_df['barcode'].tolist()
            }

        ## add the other biotypes to the biotypes list for intergenic mapping reads
        biotypes = biotypes + ['unannotated']

        if os.path.exists(f'{bampath}{sample}.merged.bam') is True:
            bam = pysam.AlignmentFile(f'{bampath}{sample}.merged.bam', 'rb')
            sample_logdata = {}
            sample_logdata['sample'] = sample
            sample_logdata['input_bam'] = f'{bampath}{sample}.merged.bam'
            sample_logdata['output_bam'] = f'{bampath}{sample}.bam'
            sample_logdata['mode'] = 'bed'
            sample_logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            sample_logdata['bed'] = self.args.bed
        else:
            print(f'An error occurred while trying to read {bampath}{sample}.merged.bam. Check the path and rerun the script.')
            sys.exit(1)
 
        ## delete output bam file if it already exists
        if os.path.exists(f'{bampath}{sample}.unsorted.bam'):
            os.remove(f'{bampath}{sample}.unsorted.bam')

        ##create new bam file to write out reads that pass filters
        try:
            output_bam = pysam.AlignmentFile(f'{bampath}{sample}.unsorted.bam', 'wb', template=bam)
        except FileNotFoundError:
            print(f'Could not write to {bampath}{sample}.unsorted.bam. Check the path and rerun the script.')
            sys.exit(1)

    
        ## define spike in regex
        ribospike_regex = r"SR_(\d+)_4N-T(\d+(\.\d+)?)"

        ## create a dictionary to store the read mappings

        # Antisense / AntisenseMulti: Unique / Multi reads whose (picked) alignment is antisense (release 0.4.0)
        read_mappability_dict = {'Mapped':0, 'Unmapped':0, 'Unique':0, 'Multi':0, 'Antisense':0, 'AntisenseMulti':0}
        read_mappings_dict = defaultdict(lambda: {'alignments': deque(), 'alignmentscores': deque(), 'annotations': deque(), 'biotypes': deque(), 'strand': deque(), 'barcodes': deque()})

        full_cov = defaultdict(lambda: {'depth': Counter(), 'mm': Counter(), 'biotype': None})

        read_feature_dict = defaultdict(lambda: 0)
        read_biotype_dict = defaultdict(lambda: 0)
        # lengths 0-300 present from the start (column order); longer reads are added (release 0.3.0: no limit at 300)
        read_length_dict = defaultdict(lambda: defaultdict(int, {i: 0 for i in range(0, 301)}))
        read_feature_distinctseqs_dict = defaultdict(lambda: defaultdict(int))
        read_biotypes_distinctseqs_dict = defaultdict(lambda: defaultdict(int))

        with tqdm.tqdm(desc=sample, unit=' reads') as pbar:

            current_read_name = None

            for read in bam.fetch(until_eof=True):

                pbar.update(1)

                read_name = read.query_name
                read_contig = read.reference_name
                read_length = len(read.query_sequence)

                if read.is_unmapped or read_length < 15:      # the 15-nt floor (the shortest insert trimming keeps)
                    count_value = 1
                    read_mappability_dict['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    read_length_dict['unmapped'][read_length] += count_value
                    output_bam.write(read)
                    continue

                if '#miND' in read_contig:
                    count_value = 1
                    read_mappability_dict['Mapped'] += count_value
                    feature = read_contig.split('.')[0].replace('#miND-', 'miND_')
                    biotype = 'miND'
                    read_feature_dict[feature] += count_value
                    read_biotype_dict[biotype] += count_value
                    read_length_dict[biotype][read_length] += count_value
                    read.set_tag('CS', feature)
                    output_bam.write(read)
                    continue

                ## check for spike in reads
                if re.search(ribospike_regex, read_contig):
                    count_value = 1
                    ribospike = re.search(ribospike_regex, read_contig).group(1)
                    sequence = read.query_sequence
                    read_length = len(sequence)

                    if read_length < (int(ribospike) - 2):
                        biotype = 'RiboSpike-Short'
                        read_biotype_dict[biotype] += count_value
                        read_length_dict[biotype][read_length] += count_value
                        read_mappability_dict['Mapped'] += count_value
                        read.set_tag('CS', 'RiboSpike-Fragment')   # 'RiboSpike-fragment' before 0.3.0
                        output_bam.write(read)
                        continue
                    else:
                        typespike = re.search(ribospike_regex, read_contig).group(2).split('.')[0]
                        typespike = 'T' + typespike
                        ribospike = f'RiboSpike-{ribospike}-{typespike}'
                        biotype = f'RiboSpike-{typespike}'
                        featuretypeD[ribospike] = biotype
                        read_feature_dict[ribospike] += count_value
                        read_biotype_dict[biotype] += count_value
                        read_length_dict[biotype][read_length] += count_value
                        read_mappability_dict['Mapped'] += count_value
                        read.set_tag('CS', ribospike)
                        output_bam.write(read)
                        continue

                if read_mappings_dict.get(read_name) is None:

                    if current_read_name is not None:
                        ## get the best alignment for the read
                        best_read, mapping, barcode, strandedness  = self.pick_bed_alignment(read_mappings_dict, current_read_name, genetype_hierarchy)
                        self._accumulate_full_coverage(read_mappings_dict[current_read_name], mapping, 1, full_cov)

                        ## check if the read is unique or multi-mapped
                        count_value = 1
                        # remove-at-source: skip counting + BAM write for multimap/antisense when requested
                        remove = (self.args.remove_multimap and mapping == 'Multi') or (self.args.remove_antisense and strandedness == 'antisense')
                        if not remove:
                            best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict = self.generate_metrics_bed_distinct(best_read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict)
                            output_bam.write(best_read)
                        # tally OUTSIDE the gate on purpose: removed multimap/antisense reads still
                        # count toward the mappingstats rate (design Z). Do NOT move these inside `if not remove`.
                        read_mappability_dict['Mapped'] += count_value
                        read_mappability_dict[mapping] += count_value
                        if strandedness == 'antisense':
                            read_mappability_dict['Antisense' if mapping == 'Unique' else 'AntisenseMulti'] += count_value
                        
                        read_mappings_dict.pop(current_read_name)

                    current_read_name = read_name

                chr = str(read.reference_name)
                read_start = read.reference_start
                read_end = read.reference_end
                read_length = len(read.query_sequence)
                read_strand = '-' if read.is_reverse else '+'
                ## reverse read strand if SMART is in the sample name
                if 'SMART' in sample:
                    read_strand = '+' if read_strand == '-' else '-'

                if chr in bed_dict:
                    d = bed_dict[chr]
                    feature, biotype, strand, barcode = self.find_feature(read_start, read_end, read_strand, d['start'], d['end'], d['biotype'], d['feature'], d['strand'], d['barcode'])
                else:      # a contig without annotation must not reach find_region
                    feature, biotype, strand, barcode = 'unannotated', 'unannotated', -2, 'unannotated'

                read_mappings_dict[read_name]['alignments'].append(read)
                ## Get the alignment score
                try:
                    alignment_score = read.get_tag('AS') ## Get the alignment score; bowtie2
                except KeyError:                       # Bowtie 1 records have no AS tag
                    alignment_score = read.get_tag('NM') * -1 ## Get the alignment score; bowtie

                read_mappings_dict[read_name]['alignmentscores'].append(alignment_score)
                read_mappings_dict[read_name]['annotations'].append(feature)
                read_mappings_dict[read_name]['biotypes'].append(biotype)
                read_mappings_dict[read_name]['strand'].append(strand)
                read_mappings_dict[read_name]['barcodes'].append(barcode)

        if current_read_name is not None:

            best_read, mapping, barcode, strandedness = self.pick_bed_alignment(read_mappings_dict, current_read_name, genetype_hierarchy)
            self._accumulate_full_coverage(read_mappings_dict[current_read_name], mapping, 1, full_cov)

            
            ## check if the read is unique or multi-mapped

            count_value = 1
            # remove-at-source: skip counting + BAM write for multimap/antisense when requested
            remove = (self.args.remove_multimap and mapping == 'Multi') or (self.args.remove_antisense and strandedness == 'antisense')
            if not remove:
                best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict = self.generate_metrics_bed_distinct(best_read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict)
                output_bam.write(best_read)
            # tally OUTSIDE the gate on purpose: removed multimap/antisense reads still
            # count toward the mappingstats rate (design Z). Do NOT move these inside `if not remove`.
            read_mappability_dict['Mapped'] += count_value
            read_mappability_dict[mapping] += count_value
            if strandedness == 'antisense':
                read_mappability_dict['Antisense' if mapping == 'Unique' else 'AntisenseMulti'] += count_value
            
            read_mappings_dict.pop(current_read_name)

        self._write_bedgraph_full(sample, full_cov)
        output_bam.close()
        bam.close()

        ## delete the defaultdict
        del read_mappings_dict

        ## from here on RNA classes instead of BED labels (class table; release 0.3.0)
        read_biotype_dict, read_length_dict, read_biotypes_distinctseqs_dict = to_classes(
            self.classes.class_of, read_biotype_dict, read_length_dict, read_biotypes_distinctseqs_dict)
        unannotated = {self.classes.class_of('unannotated'), self.classes.class_of('UniVec')}   # reads without a feature


        read_mappability_dict['Total'] = read_mappability_dict['Mapped'] + read_mappability_dict['Unmapped']

        sample_logdata['mapping_total'] = read_mappability_dict['Total']
        sample_logdata['mapping_mapped'] = read_mappability_dict['Mapped']
        sample_logdata['mapping_unmapped'] = read_mappability_dict['Unmapped']
        sample_logdata['mapping_unique'] = read_mappability_dict['Unique']
        sample_logdata['mapping_multi'] = read_mappability_dict['Multi']

        read_mappability_df = pd.DataFrame.from_dict(read_mappability_dict, orient='index').rename(columns={0:sample})

        read_length_df = pd.DataFrame.from_dict(read_length_dict, orient='index').T.fillna(0).astype(int)
        read_length_df['total'] = read_length_df.sum(axis=1)
        lengths = read_length_df[['total']].T
        read_length_df = read_length_df.drop(columns=['total']).T

        ## get the average read length
        total_lengths = []
        for length in lengths.columns:
            number_of_reads = lengths[length].values[0]
            ## add the length to the total lengths list number_of_reads times
            total_lengths.extend([length] * number_of_reads)

        median_read_length = np.median(total_lengths)
        sample_logdata['median_read_length'] = median_read_length
        del total_lengths


        read_feature_df = pd.DataFrame(read_feature_dict, index=[0]).T.rename(columns={0:sample})

        read_biotype_df = pd.DataFrame(read_biotype_dict, index=[0]).T.rename(columns={0:sample})

        sample_logdata['reads_annotated'] = read_biotype_df[~read_biotype_df.index.isin(unannotated)].sum().values[0]
        sample_logdata['reads_unannotated'] = read_biotype_df[read_biotype_df.index.isin(unannotated)].sum().values[0]
        sample_logdata['reads_total'] = sample_logdata['reads_annotated'] + sample_logdata['reads_unannotated']

        sample_logdata_df = pd.DataFrame.from_dict(sample_logdata, orient='index')

        ## make output directory if it does not exist (exist_ok: samples run in parallel,
        ## a separate exists-check followed by mkdir can race and raise FileExistsError)
        os.makedirs(f'{self.args.output}/metrics', exist_ok=True)

        sample_logdata_df.to_csv(f'{self.args.output}/metrics/{sample}.log', header=False)

        read_features = read_feature_df.index.tolist()
        read_feature_biotypes = [self.classes.class_of(featuretypeD.get(i)) for i in read_features]
        read_featuretypes_df = pd.DataFrame(read_feature_biotypes, index=read_features, columns=['rna_class']).reset_index().rename(columns={'index':'feature'})

        ## get a list of sorted biotypes by the count in biotype_df
        sorted_biotypes = read_biotype_df.sort_values(by=sample, ascending=False, kind='stable').index.tolist()

        read_biotype_df = read_biotype_df.loc[sorted_biotypes]

        read_length_df = read_length_df.loc[sorted_biotypes]
        read_length_df['sample'] = sample
        read_length_df = read_length_df.reset_index().rename(columns={'index':'rna_class'})
        longest = max([150] + [length for length in lengths.columns if lengths[length].values[0] > 0])
        read_length_df = read_length_df.reindex(columns=['sample', 'rna_class'] + list(range(0, longest + 1)), fill_value=0)

        ## Get the number of distinct sequences for each feature by finding the length of the dictionary
        ## loop through the dictionary, find the length of the dictionary to make a dataframe of feature,length_of_distinct_sequences
        read_feature_distinctseqs_df = pd.DataFrame({k: [len(v)] for k, v in read_feature_distinctseqs_dict.items()}).T.rename(columns={0:'distinct_sequences'}).reset_index().rename(columns={'index':'feature'})
        read_biotypes_distinctseqs_df = pd.DataFrame({k: [len(v)] for k, v in read_biotypes_distinctseqs_dict.items()}).T.rename(columns={0:'distinct_sequences'}).reset_index().rename(columns={'index':'rna_class'})
        read_biotypes_distinctseqs_df['sample'] = sample
        read_feature_distinctseqs_df['sample'] = sample
        read_feature_distinctseqs_df = read_feature_distinctseqs_df.pivot(index='feature', columns='sample', values='distinct_sequences').fillna(0).astype(int)
        read_biotypes_distinctseqs_df = read_biotypes_distinctseqs_df.pivot(index='rna_class', columns='sample', values='distinct_sequences').fillna(0).astype(int)

        del read_feature_dict
        del read_biotype_dict
        del read_length_dict
        del read_feature_distinctseqs_dict
        del read_biotypes_distinctseqs_dict

        return [read_feature_df, read_featuretypes_df, read_biotype_df, read_length_df, read_mappability_df, read_feature_distinctseqs_df, read_biotypes_distinctseqs_df]
    
    
    def _accumulate_full_coverage(self, group, mapping, count_value, full_cov):
        """Task 4 — full-coverage bedgraph. Every EQUALLY-BEST locus of the read
        contributes its full weight (no partitioning); a read flagged `Multi` also
        adds that weight to multimap_depth. Dedup per (read, feature, strand, pos)."""
        scores = group['alignmentscores']
        if not scores:
            return
        best = max(scores)
        is_multi = (mapping == 'Multi')
        covered = set()
        biotype_of = {}
        for i, sc in enumerate(scores):
            if sc != best:                      # only the equally-best loci (true footprint)
                continue
            feat = group['annotations'][i]
            if feat == 'unannotated':
                continue
            aln = group['alignments'][i]
            if aln.reference_start is None or aln.reference_end is None:
                continue
            strand_label = 'sense' if group['strand'][i] in (1, -2) else 'antisense'
            biotype_of[feat] = group['biotypes'][i]
            for p in range(aln.reference_start, aln.reference_end):
                covered.add((feat, strand_label, p))
        for feat, strand_label, p in covered:
            entry = full_cov[(feat, strand_label)]
            entry['depth'][p] += count_value
            if is_multi:
                entry['mm'][p] += count_value
            entry['biotype'] = biotype_of[feat]

    def _write_bedgraph_full(self, sample, full_cov):
        """Write bedgraph_full/{sample}.bedgraph:
        feature  start  end  depth  multimap_depth  strand  rna_class (class id; release 0.3.0)
        (run-length-encoded on constant (depth, multimap_depth) blocks)."""
        outdir = f"{self.args.output}/bedgraph_full"
        os.makedirs(outdir, exist_ok=True)
        with open(f"{outdir}/{sample}.bedgraph", 'w') as f:
            # sorted: full_cov is filled from a set, whose order changes between runs
            # (string hash randomisation); sorting makes the file byte-identical (fixed 2026-09-29)
            for (feat, strand_label), entry in sorted(full_cov.items(), key=lambda kv: kv[0]):
                depth, mm, biotype = entry['depth'], entry['mm'], self.classes.class_of(entry['biotype'])
                if not depth:
                    continue
                positions = sorted(depth)
                s = prev = positions[0]
                d0, m0 = depth[s], mm.get(s, 0)
                for p in positions[1:]:
                    d, m = depth[p], mm.get(p, 0)
                    if p == prev + 1 and d == d0 and m == m0:
                        prev = p
                    else:
                        f.write(f"{feat}\t{s}\t{prev+1}\t{d0}\t{m0}\t{strand_label}\t{biotype}\n")
                        s = prev = p; d0, m0 = d, m
                f.write(f"{feat}\t{s}\t{prev+1}\t{d0}\t{m0}\t{strand_label}\t{biotype}\n")

    def find_feature(self, read_start, read_end, read_strand, regions_starts, regions_ends, regions_biotypes, regions_features, regions_strands, regions_barcodes):

        ## check for overlap between read and annotation where the read overlaps with the annotation with 2 bases of clipping
        index = find_region(regions_starts, regions_ends, read_start + self.args.clipping, read_end - self.args.clipping)

        if index == -1:
            return 'unannotated', 'unannotated', -2, 'unannotated'
        else:
            biotype = regions_biotypes[index]
            strandedness = 1 if regions_strands[index] == read_strand else -1
            name = regions_features[index].strip()
            barcode = regions_barcodes[index]
            return name, biotype, strandedness, barcode

    def sort_and_index(self):

        bampath = self.args.bampath

        for sample in self.samplelist:

            if os.path.exists(f'{bampath}/{sample}.bam.bai'):
                run(f'rm -f {bampath}/{sample}.bam.bai')
            if os.path.exists(f'{bampath}/{sample}.bam.csi'):
                run(f'rm -f {bampath}/{sample}.bam.csi')

            cmd = f'samtools sort -@ {self.args.threads} -o {bampath}/{sample}.bam {bampath}/{sample}.unsorted.bam'
            run(cmd)

            try:
                bamfile = pysam.AlignmentFile(f'{bampath}/{sample}.bam', 'rb')
                need_csi = any([len(ref) > 512000000 for ref in bamfile.header.references])
                bamfile.close()
            except FileNotFoundError:
                print(f'{bampath}{sample}.bam does not exist. Check the path and rerun the script.')
                sys.exit(1)

            if need_csi:
                run(f'samtools index -c -@ {self.args.threads} {bampath}/{sample}.bam')
            else:
                run(f'samtools index -@ {self.args.threads} {bampath}/{sample}.bam')

            run(f'rm -f {bampath}/{sample}.unsorted.bam')
        
        return

    def generate_module_logfile(self):
        '''
        Generates a log file for the module.
        '''
        log_file = f'{self.args.output}/logs/parse.log'
        with open(log_file, 'w') as f:
            f.write(f'###\n')
            for key, value in self.logdata.items():
                if key == 'input_bams':
                    f.write(f'{key}:\n')
                    for sample, bam in value.items():
                        f.write(f'{sample}: {bam}\n')
                else:
                    f.write(f'{key}: {value}\n')
            f.write('###\n')
        f.close()
        return

def main():

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('-s','--samples', required=True, help='Sample metadata file. Format: sample,group,fastq')
    ap.add_argument('-@','--threads', required=False, help='Number of threads to process paraellel fastq files with. Default=4', default=4)
    ap.add_argument('-o', '--output', required=True, help='Output prefix.')
    ap.add_argument('-b','--bed', required=True, help='BED file to assign reads to')
    ap.add_argument('-y','--bampath', required=False, help='Path to the bam files', default='')
    ap.add_argument('-z','--clipping', required=False, help='Number of bases to clip from the start and end of the read. Default=2', default=2)    
    ap.add_argument('--remove-multimap', required=False, default=False, action='store_true', help='Drop multimapping reads from counts, featuretypes and the parsed BAM (still tallied in mappingstats).')
    ap.add_argument('--gene-symbols', required=False, default=None, help='gene-symbol table (docker/tables/<species>.gene_symbols.tsv); adds its columns to featureclasses')
    ap.add_argument('--contig-tiers', required=False, default=None, help='contig-to-tier map (contig<TAB>tier; Snakefile: <run>/tables/contig_tiers.tsv). Without it, reads on UniVec contigs stay genome_unannotated')
    ap.add_argument('--rna-classes', required=False, default=None, help='RNA class table (docker/scripts/refTables.py); default: docker/tables/default.rna_classes.tsv')
    ap.add_argument('--remove-antisense', required=False, default=False, action='store_true', help='Drop antisense reads from counts, featuretypes and the parsed BAM (still tallied in mappingstats).')
    args = ap.parse_args()
    
    startTime = datetime.now()

    pm = ParseMappings(args)

    print(f'-'*50)
    print(f'-'*50)
    print(f'Parsing mappings.')
    print(f'-'*50)
    print(f'-'*50)

    if args.bampath == '':
        args.bampath = f'{args.output}/alignments'
    
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(f'{args.output}/alignments'):
        os.mkdir(f'{args.output}/alignments')

    samplelist = pm.samplelist

    bampath = args.bampath.rstrip('/') + '/'

    for sample in samplelist:
        pm.logdata['input_bams'][sample] = f'{bampath}{sample}.merged.bam'


    # parse_mapping_bed_distinct is the only assignment path: the older parse_mapping_bed, unused because the
    # Snakefile always passed -q, was removed in 0.3.0 (-q itself in 0.4.0). BED labels are checked once here.
    if os.path.exists(args.bed) is False:
        print(f'{args.bed} does not exist. Check the path and rerun the script.')
        sys.exit(1)
    else:
        pm.classes.check_labels(pd.read_csv(args.bed, sep='\t', header=None, usecols=[4])[4].drop_duplicates().tolist())
        feature_data = []
        featuretype_data = []
        biotype_data = []
        length_data = []
        mapdata = []
        featuredistseqs_data = []
        biotypedistseqs_data = []
        with Pool(int(args.threads)) as p:

            sampledata = p.map(pm.parse_mapping_bed_distinct, samplelist)
            for data in sampledata:
                feature_data.append(data[0])
                featuretype_data.append(data[1])
                biotype_data.append(data[2])
                length_data.append(data[3])
                mapdata.append(data[4])
                featuredistseqs_data.append(data[5])
                biotypedistseqs_data.append(data[6])

    ## write the feature stats to a file
    feature_df = pd.concat(feature_data, axis=1).fillna(0).astype(int).sort_index(kind='stable')
    feature_df.index.name = 'feature'

    # The spike-in counts are in ribospike/<prefix>-ribospike.csv (countRiboSpike.py); the copy this step wrote
    # next to the read counts was removed in 0.4.0 (cleanup; the RiboSpike-* rows of *-readcounts.csv remain)

    biotype_df = pd.concat(biotype_data, axis=1).fillna(0).astype(int)
    biotype_df.index.name = 'rna_class'

    # columns 0 to the longest read of the run, at least 150 (release 0.3.0; 0-150 before)
    length_df = pd.concat(length_data, axis=0)
    length_columns = sorted(c for c in length_df.columns if c not in ('sample', 'rna_class'))
    length_df = length_df.reindex(columns=['sample', 'rna_class'] + length_columns)
    length_df[length_columns] = length_df[length_columns].fillna(0).astype(int)

    featuretypes = pd.concat(featuretype_data, axis=0).drop_duplicates(subset=['feature']).sort_values(by='feature', kind='stable')
    featuretypes = add_gene_length(featuretypes, args.bed)
    if args.gene_symbols:
        featuretypes = add_gene_symbols(featuretypes, args.bed, args.gene_symbols)

    mappability_df = pd.concat(mapdata, axis=1).fillna(0).astype(int).T   # integers (0.4.0)
    mappability_df['Total'] = mappability_df['Mapped'] + mappability_df['Unmapped']
    mappability_df.index.name = 'sample'


    feature_distinctseqs_df = pd.concat(featuredistseqs_data, axis=1).fillna(0).astype(int)
    biotype_distinctseqs_df = pd.concat(biotypedistseqs_data, axis=1).fillna(0).astype(int)
    
    print(f'-'*50)
    feature_df.to_csv(f'{args.output}/{args.output}-readcounts.csv')
    # file names of release 0.3.0 (before: typecounts, typelengths, featuretypes, biotypedistinctseqs)
    biotype_df.to_csv(f'{args.output}/{args.output}-classcounts.csv')
    length_df.to_csv(f'{args.output}/{args.output}-classlengths.csv', index=False)
    featuretypes.to_csv(f'{args.output}/{args.output}-featureclasses.csv', index=False)
    mappability_df.to_csv(f'{args.output}/{args.output}-mappingstats.csv')

    ## write the distinct sequences to a file
    feature_distinctseqs_df.to_csv(f'{args.output}/{args.output}-featuredistinctseqs.csv')
    biotype_distinctseqs_df.to_csv(f'{args.output}/{args.output}-classdistinctseqs.csv')

    print(f'Read counts written to {args.output}/{args.output}-readcounts.csv')
    print(f'Feature classes written to {args.output}/{args.output}-featureclasses.csv')
    print(f'Class counts written to {args.output}/{args.output}-classcounts.csv')
    print(f'Class read lengths written to {args.output}/{args.output}-classlengths.csv')
    print(f'Mapping stats written to {args.output}/{args.output}-mappingstats.csv')
    print(f'-'*50)

    pm.sort_and_index()
    pm.generate_module_logfile()    

    print(f'-'*50)
    print(f'-'*50)
    print('Time elasped for read parsing: ', datetime.now() - startTime)
    print(f'-'*50)
    print(f'-'*50)

if __name__ == "__main__":
    main()
