#!/usr/bin/env python3
#Aidan Manning 7.9.24
#Script to map small RNA reads to the genome, ncRNAs, and other RNAs using various mapping tools
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import warnings
from multiprocessing import Pool
from toolRun import run  # stops on a failed tool
import pysam
from collections import defaultdict, deque, Counter, defaultdict
import tqdm
import numpy as np
np.random.seed(42)
import re
import random
random.seed(42)
from numba import jit
warnings.filterwarnings("ignore")

def add_gene_length(featuretypes, bed_path):
    """Attach a gene_length column: the base-feature length from the BED.
    parseMappings suffixes feature names (_5p/_3p/_internal/_whole, _multimap,
    _antisense); strip them so each row matches its base feature in the BED.
    """
    bed = pd.read_csv(bed_path, sep='\t', header=None,
                      usecols=[1, 2, 3], names=['start', 'end', 'feature'])
    base_length = (bed['end'] - bed['start']).groupby(bed['feature']).max().to_dict()

    suffixes = ('_multimap', '_antisense', '_5p', '_3p', '_internal', '_whole')
    def base_name(name):
        name = name.split('|')[-1]
        changed = True
        while changed:
            changed = False
            for suf in suffixes:
                if name.endswith(suf):
                    name = name[:-len(suf)]
                    changed = True
        return name

    featuretypes = featuretypes.copy()
    featuretypes['gene_length'] = featuretypes['feature'].map(
        lambda f: base_length.get(base_name(f)))
    return featuretypes

@jit(nopython=True)
def find_region(starts, ends, read_start, read_end):
    for i in range(len(starts)):
        # Check for any overlap between the intervals
        if read_start <= ends[i] and read_end >= starts[i]:
            return i  # Return index of the first overlapping region
    return -1  # Indicate no overlap was found



class ParseMappings(object):

    def __init__(self, args):
        '''
        This class is used to parse the mappings from the various mapping tools.
        '''

        self.args = args

        self.sampleinfo = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq'])
        self.samplelist = self.sampleinfo['sample'].tolist()
        self.sampleDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['group']))
        self.fastqDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['fastq']))

        self.logdata = {}
        self.logdata['module'] = 'parseMappings'
        self.logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.logdata['samples'] = self.samplelist
        self.logdata['threads'] = self.args.threads
        self.logdata['output'] = self.args.output
        self.logdata['mode'] = 'bed' if self.args.bed else 'best'
        if self.args.bed:
            self.logdata['bed'] = self.args.bed
            self.logdata['mode'] = 'bed'
        else:
            self.logdata['mode'] = 'best'
        self.logdata['tools'] = []
        self.logdata['commands'] = []
        self.logdata['input_bams'] = {}

        return

    def partition_read_counts(self, reads):
        '''
        This function partitions the read counts based on the number of alignments.
        '''
        def partition_number(n, parts):
            # Divide the number by parts and calculate the integer part
            quotient = n // parts
            # Calculate the remainder
            remainder = n % parts
            # Create the partition
            partition = [quotient] * parts
            # Distribute the remainder across the first 'remainder' parts
            for i in range(remainder):
                partition[i] += 1
            return partition
        ## Get the read name
        read_name = reads[0].query_name
        ## get the read count
        read_count = int(read_name.split(';size=')[1])
        ## get the number of alignments
        num_alignments = len(reads)

        count_partitions = partition_number(read_count, num_alignments)

        return count_partitions

    def get_read_count(self, read_name):
        '''
        This function gets the read count from the read name if the reads are collapsed.
        '''
        if self.args.collapse == True:
            read_count = int(read_name.split(';size=')[1])
        elif self.args.collapse == False:
            read_count = int(1)

        # if ';size=' in read_name:
        #     read_count = int(read_name.split(';size=')[1])
        # else:
        #     read_count = int(1)

        return read_count

    def parse_mapping_bed(self, sample):
        '''
        Assigns reads to bed features.
        '''

        if self.args.bampath == '':
            bampath = f'{self.args.output}/alignments/'
        else:
            if self.args.bampath.endswith('/') is False:
                bampath = self.args.bampath + '/'

        genetype_hierarchy = {
            'unannotated':0,
            'protein_coding':1,
            'pre-mRNA':1,
            'sinRNA':1,
            'scdRNA':1,
            'sutRNA':1,
            'piRNA':2,
            'lncRNA':2,
            'SINE':3,
            'LINE':3,
            'pre-tRNA':4,
            'IG_V_pseudogene': 4,
            'IG_V_gene': 5,
            'IG_C_gene': 5,
            'IG_J_gene': 5,
            'TR_V_gene': 5,
            'TR_V_pseudogene': 4,
            'TR_D_gene': 5,
            'TR_J_gene': 5,
            'TR_C_gene': 5,
            'IG_C_pseudogene': 4,
            'TR_J_pseudogene': 4,
            'IG_J_pseudogene': 4,
            'IG_D_gene': 5,
            'Y_RNA':5,
            'vault_RNA':5,
            'scRNA': 5,
            'scaRNA': 5,
            'snoRNA':5,
            'snRNA':5,
            'pre-miRNA':4,
            'star-miRNA':4,
            'Mt_rRNA':6,
            'miRNA':6,
            'rRNA':7,
            'Mt_tRNA':8,
            'tRNA':9
        }

        try:
            bed = pd.read_csv(self.args.bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'feature', 'biotype', 'strand'])
        except pd.errors.EmptyDataError:
            print(f'{self.args.bed} is empty. Check the file and rerun the script.')
            sys.exit(1)
        except FileNotFoundError:
            print(f'{self.args.bed} does not exist. Check the path and rerun the script.')
            sys.exit(1)

        ## creates a unique identifier for each feature
        bed['barcode'] = bed['feature'].astype(str) + '::' + bed['chrom'].astype(str) + ':' + bed['start'].astype(str) + '-' + bed['end'].astype(str) + ':' + bed['strand']
        ## get the unique biotypes in the bed file
        biotypes = bed['biotype'].drop_duplicates().tolist()

        ## make a barcodes dictionary
        bed_barcode_dict = bed.set_index('barcode').T.to_dict('dict')

        uniquefeatures = bed.drop_duplicates(subset=['feature'], keep='first')
        featuretypeD = dict(zip(uniquefeatures['feature'], uniquefeatures['biotype']))
        del uniquefeatures

        ## get the contigs
        contigs = {contig:1 for contig in bed['chrom'].unique()}

        ## get the other biotypes in the bed file
        other_biotypes = [i for i in biotypes if i not in genetype_hierarchy.keys()]

        ## add the other biotypes to the genetype_hierarchy, assume they are important so rank them higher than the other biotypes
        if len(other_biotypes) > 0:
            for i in other_biotypes:
                genetype_hierarchy[i] = 10

        ## assign a sorting order to each feature based on the genetype_hierarchy
        bed['order'] = bed['biotype'].map(genetype_hierarchy).fillna(0)

        ## sort the bed file based on the sorting order
        bed = bed.sort_values(by=['order'], ascending=False).drop(columns=['order'])

        ## group the bed file by chromosome
        bed = bed.groupby('chrom')

        bed_strand_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {'start': list(), 'end': list(), 'feature': list(), 'biotype': list(), 'strand': list()}))))
        
        ## create a dictionary to store the annotations for each chromosome so that they can be accessed quickly
        for chrom, group_df in bed:
            for strand in ['+', '-']:
                bed_strand_dict[chrom][strand] = {
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
            # self.logdata['input_bams'].append([sample, f'{bampath}{sample}.merged.bam'])
            sample_logdata = {}
            sample_logdata['sample'] = sample
            sample_logdata['input_bam'] = f'{bampath}{sample}.merged.bam'
            sample_logdata['output_bam'] = f'{bampath}{sample}.bam'
            if self.args.unmapped:
                sample_logdata['unmapped_bam'] = f'{bampath}{sample}.unmapped.bam'
            sample_logdata['mode'] = 'bed'
            sample_logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            sample_logdata['bed'] = self.args.bed
        else:
            print(f'An error occurred while trying to read {bampath}{sample}.merged.bam. Check the path and rerun the script.')
            sys.exit(1)
 
        ## delete output bam file if it already exists
        if os.path.exists(f'{bampath}{sample}.unsorted.bam'):
            os.remove(f'{bampath}{sample}.unsorted.bam')

        if self.args.unmapped:
            if os.path.exists(f'{bampath}{sample}.unmapped.bam'):
                os.remove(f'{bampath}{sample}.unmapped.bam')
            ## make unmapped metrics dir
            if not os.path.exists(f'{self.args.output}/unmapped'):
                os.mkdir(f'{self.args.output}/unmapped')

        ##create new bam file to write out reads that pass filters
        try:
            output_bam = pysam.AlignmentFile(f'{bampath}{sample}.unsorted.bam', 'wb', template=bam)
        except FileNotFoundError:
            print(f'Could not write to {bampath}{sample}.unsorted.bam. Check the path and rerun the script.')
            sys.exit(1)

        if self.args.unmapped:
            try:
                unmapbam = pysam.AlignmentFile(f'{bampath}{sample}.unmapped.bam', 'wb', template=bam)
                unmap_sequences = defaultdict(lambda: 0)
            except FileNotFoundError:
                print(f'Could not write to {bampath}{sample}.unmapped.bam. Check the path and rerun the script.')
                sys.exit(1)
    
        ## check if the reads are collapsed, on a SEPARATE file handle: reading the first
        ## record of `bam` here made the main loop skip it, because fetch(until_eof=True)
        ## continues from the current file position (1 read lost per sample; fixed 2026-09-29)
        with pysam.AlignmentFile(f'{bampath}{sample}.merged.bam', 'rb') as peek:
            first_read = next(peek.fetch(until_eof=True), None)
        if first_read is not None:
            self.args.collapse = ';size=' in first_read.query_name

        ## define spike in regex
        ribospike_regex = r"SR_(\d+)_4N-T(\d+(\.\d+)?)"
        ribospike_dict = defaultdict(int)

        ## create a dictionary to store the read mappings
        read_mappability_dict = {'Mapped':0, 'Unmapped':0, 'Unique':0, 'Multi':0}
        read_mappings_dict = defaultdict(lambda: {'alignments': deque(), 'alignmentscores': deque(), 'annotations': deque(), 'biotypes': deque(), 'strand': deque(), 'barcodes': deque()})

        read_feature_dict = defaultdict(lambda: 0)
        read_biotype_dict = defaultdict(lambda: 0)
        read_length_dict = defaultdict(lambda: {i:0 for i in range(0, 301)})
        
        with tqdm.tqdm(desc=sample, unit=' reads') as pbar:

            current_read_name = None

            for read in bam.fetch(until_eof=True):

                pbar.update(1)

                read_name = read.query_name
                read_contig = read.reference_name

                if read.is_unmapped:
                    count_value = self.get_read_count(read_name)
                    read_mappability_dict['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    read_length_dict['unmapped'][read.query_length] += count_value
                    output_bam.write(read)
                    if self.args.unmapped:
                        unmapbam.write(read)
                        unmap_sequences[read.query_sequence] += count_value
                    continue

                if read.query_length < 15:
                    count_value = self.get_read_count(read_name)
                    read_mappability_dict['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    read_length_dict['unmapped'][read.query_length] += count_value
                    output_bam.write(read)
                    if self.args.unmapped:
                        unmapbam.write(read)
                        unmap_sequences[read.query_sequence] += count_value
                    continue

                ## check for spike in reads
                if re.search(ribospike_regex, read_contig):
                    count_value = self.get_read_count(read_name)
                    ribospike = re.search(ribospike_regex, read_contig).group(1)
                    sequence = read.query_sequence
                    read_length = len(sequence)

                    if read_length < (int(ribospike) - 2):
                        biotype = 'RiboSpike-Short'
                        read_biotype_dict[biotype] += count_value
                        read_length_dict[biotype][read_length] += count_value
                        read_mappability_dict['Mapped'] += count_value
                        read.set_tag('CS', 'RiboSpike-fragment')
                        output_bam.write(read)
                        continue
                    else:
                        typespike = re.search(ribospike_regex, read_contig).group(2).split('.')[0]
                        typespike = 'T' + typespike
                        ribospike = f'RiboSpike-{ribospike}-{typespike}'
                        biotype = f'RiboSpike-{typespike}'
                        featuretypeD[ribospike] = biotype
                        ribospike_dict[ribospike] += count_value
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
                        if mapping == 'Multi':
                            continue
                        ## check if the read is unique or multi-mapped
                        if mapping == 'Multi' and self.args.collapse == True:
                            count_values = self.partition_read_counts(best_read)
                            for cdx, count_value in enumerate(count_values):
                                if count_value == 0:
                                    continue
                                curr_best_read = best_read[cdx]
                                ## set the size= in the read name to the count value
                                oldname = curr_best_read.query_name
                                newname = oldname.split(';size=')[0] + f';size={count_value}'
                                curr_best_read.query_name = newname
                                curr_barcode = barcode[cdx]
                                curr_strandedness = strandedness[cdx]
                                curr_best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD = self.generate_metrics_bed(curr_best_read, curr_barcode, curr_strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value)
                                output_bam.write(curr_best_read)
                                read_mappability_dict['Mapped'] += count_value
                                read_mappability_dict[mapping] += count_value

                        else:
                            count_value = self.get_read_count(current_read_name)
                            best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD = self.generate_metrics_bed(best_read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value)
                            output_bam.write(best_read)
                            read_mappability_dict['Mapped'] += count_value
                            read_mappability_dict[mapping] += count_value
                            
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

                if contigs.get(chr) == 1:

                    pos_feature, pos_biotype, pos_strand, pos_barcode = self.find_feature(read_start, read_end, read_strand, bed_strand_dict[chr]['+']['start'], bed_strand_dict[chr]['+']['end'], bed_strand_dict[chr]['+']['biotype'], bed_strand_dict[chr]['+']['feature'], bed_strand_dict[chr]['+']['strand'], bed_strand_dict[chr]['+']['barcode'])
                    neg_feature, neg_biotype, neg_strand, neg_barcode = self.find_feature(read_start, read_end, read_strand, bed_strand_dict[chr]['-']['start'], bed_strand_dict[chr]['-']['end'], bed_strand_dict[chr]['-']['biotype'], bed_strand_dict[chr]['-']['feature'], bed_strand_dict[chr]['-']['strand'], bed_strand_dict[chr]['-']['barcode'])

                    if pos_feature == 'unannotated' and neg_feature == 'unannotated':
                        feature, biotype, strand, barcode = 'unannotated', 'unannotated', -2, 'unannotated'

                    elif pos_feature == 'unannotated' and neg_feature != 'unannotated':
                        feature, biotype, strand, barcode = neg_feature, neg_biotype, neg_strand, neg_barcode
                    
                    elif pos_feature != 'unannotated' and neg_feature == 'unannotated':
                        feature, biotype, strand, barcode = pos_feature, pos_biotype, pos_strand, pos_barcode

                    elif pos_feature != 'unannotated' and neg_feature != 'unannotated':
                        ## go with the highest strand
                        if pos_strand > neg_strand:
                            feature, biotype, strand, barcode = pos_feature, pos_biotype, pos_strand, pos_barcode
                        elif pos_strand < neg_strand:
                            feature, biotype, strand, barcode = neg_feature, neg_biotype, neg_strand, neg_barcode
                        elif pos_strand == neg_strand:
                            ## go with higher based on the genetype_hierarchy
                            if genetype_hierarchy[pos_biotype] > genetype_hierarchy[neg_biotype]:
                                feature, biotype, strand, barcode = pos_feature, pos_biotype, pos_strand, pos_barcode
                            elif genetype_hierarchy[pos_biotype] < genetype_hierarchy[neg_biotype]:
                                feature, biotype, strand, barcode = neg_feature, neg_biotype, neg_strand, neg_barcode
                            elif genetype_hierarchy[pos_biotype] == genetype_hierarchy[neg_biotype]:
                                # deterministic: pick the lexicographically smaller barcode
                                options = [(pos_feature, pos_biotype, pos_strand, pos_barcode),
                                        (neg_feature, neg_biotype, neg_strand, neg_barcode)]
                                feature, biotype, strand, barcode = min(options, key=lambda o: o[3])

                else:
                    # print(read_contig)
                    feature, biotype, strand, barcode = 'unannotated', 'unannotated', -2, 'unannotated'

                read_mappings_dict[read_name]['alignments'].append(read)
                ## Get the alignment score
                try:
                    alignment_score = read.get_tag('AS') ## Get the alignment score; bowtie2
                except:
                    alignment_score = read.get_tag('NM') * -1 ## Get the alignment score; bowtie

                read_mappings_dict[read_name]['alignmentscores'].append(alignment_score)
                read_mappings_dict[read_name]['annotations'].append(feature)
                read_mappings_dict[read_name]['biotypes'].append(biotype)
                read_mappings_dict[read_name]['strand'].append(strand)
                read_mappings_dict[read_name]['barcodes'].append(barcode)

        if current_read_name is not None:

            best_read, mapping, barcode, strandedness = self.pick_bed_alignment(read_mappings_dict, current_read_name, genetype_hierarchy)


            ## check if the read is unique or multi-mapped
            if mapping == 'Multi' and self.args.collapse == True:
                count_values = self.partition_read_counts(best_read)
                for cdx, count_value in enumerate(count_values):
                    if count_value == 0:
                        continue
                    curr_best_read = best_read[cdx]
                    ## set the size= in the read name to the count value
                    oldname = curr_best_read.query_name
                    newname = oldname.split(';size=')[0] + f';size={count_value}'
                    curr_best_read.query_name = newname
                    curr_barcode = barcode[cdx]
                    curr_strandedness = strandedness[cdx]
                    curr_best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD = self.generate_metrics_bed(curr_best_read, curr_barcode, curr_strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value)
                    output_bam.write(curr_best_read)
                    read_mappability_dict['Mapped'] += count_value
                    read_mappability_dict[mapping] += count_value
            else:
                count_value = self.get_read_count(current_read_name)
                best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD = self.generate_metrics_bed(best_read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value)
                output_bam.write(best_read)
                read_mappability_dict['Mapped'] += count_value
                read_mappability_dict[mapping] += count_value
                
            read_mappings_dict.pop(current_read_name)

        output_bam.close()
        bam.close()

        if self.args.unmapped:
            unmapbam.close()

        ## delete the defaultdict
        del read_mappings_dict


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

        # if self.args.unmapped:
        #     unmap_sequences_df = pd.DataFrame.from_dict(unmap_sequences, orient='index').T.fillna(0).astype(int).sort_values(by=0, ascending=False).head(100)
        #     unmap_sequences_df.to_csv(f'{bampath}{sample}-unmappedtop100.csv')
        #     del unmap_sequences

        read_feature_df = pd.DataFrame(read_feature_dict, index=[0]).T.rename(columns={0:sample})

        read_biotype_df = pd.DataFrame(read_biotype_dict, index=[0]).T.rename(columns={0:sample})

        sample_logdata['reads_annotated'] = read_biotype_df[read_biotype_df.index != 'unannotated'].sum().values[0]
        sample_logdata['reads_unannotated'] = read_biotype_df[read_biotype_df.index == 'unannotated'].sum().values[0]
        sample_logdata['reads_total'] = sample_logdata['reads_annotated'] + sample_logdata['reads_unannotated']

        sample_logdata_df = pd.DataFrame.from_dict(sample_logdata, orient='index')

        ## make output directory if it does not exist (exist_ok: samples run in parallel,
        ## a separate exists-check followed by mkdir can race and raise FileExistsError)
        os.makedirs(f'{self.args.output}/metrics', exist_ok=True)

        sample_logdata_df.to_csv(f'{self.args.output}/metrics/{sample}.log', header=False)

        read_features = read_feature_df.index.tolist()
        read_feature_biotypes = [featuretypeD.get(i) for i in read_features]
        read_featuretypes_df = pd.DataFrame(read_feature_biotypes, index=read_features, columns=['feature_type']).reset_index().rename(columns={'index':'feature'})

        ## get a list of sorted biotypes by the count in biotype_df
        sorted_biotypes = read_biotype_df.sort_values(by=sample, ascending=False).index.tolist()

        read_biotype_df = read_biotype_df.loc[sorted_biotypes]

        read_length_df = read_length_df.loc[sorted_biotypes]
        read_length_df['sample'] = sample
        read_length_df = read_length_df.reset_index().rename(columns={'index':'biotype'})
        read_length_df = read_length_df[['sample', 'biotype'] + [i for i in range(0, 151)]]

        # for k,v in ribospike_end_dict.items():
        #     df = pd.DataFrame.from_dict(v, orient='index').T.fillna(0).astype(int)
        #     df = df[['5p_1', '5p_2', '5p_3', '5p_4', '3p_4', '3p_3', '3p_2', '3p_1']]
        #     df = df.div(df.sum(axis=0), axis=1) * 100
        #     df = df.stack().reset_index().rename(columns={'level_0':'nucleotide', 'level_1':'position', 0:'count'})
        #     df['sample'] = sample
        #     df['spike'] = k

        return [read_feature_df, read_featuretypes_df, read_biotype_df, read_length_df, read_mappability_df]

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
            elif len(highest_alignscore_indices) > 1 and self.args.collapse == True:
                mappability = 'Multi'
                reads = [alignments[i] for i in highest_alignscore_indices]
                strands = [strands[i] for i in highest_alignscore_indices]
                barcodes = [barcodes[i] for i in highest_alignscore_indices]
                strandednesses = [assign_strandedness_tag(reads[i], strands[i], barcodes[i], mappability) for i in range(len(highest_alignscore_indices))]
                return reads, mappability, barcodes, strandednesses
            elif len(highest_alignscore_indices) > 1 and self.args.collapse == False:
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

    def generate_metrics_bed(self, read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value):

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
            
            if biotype in ['tRNA', 'Mt_tRNA']:

                # if biotype in ['Mt_tRNA']:
                #     biotype = 'tRNA'

                if strandedness == 'antisense':
                    biotype = 'antisense'
                    feature = feature + '_antisense'
                else:
                    feature = feature

                if read_start <= start + 10 and read_end < end - 10:
                    feature = feature + '_5p'
                elif read_start > start + 10 and read_end >= end - 10:
                    feature = feature + '_3p'
                elif read_start > start + 10 and read_end < end - 10:
                    feature = feature + '_internal'
                elif read_start <= start + 10 and read_end >= end - 10:
                    feature = feature + '_whole'

            elif biotype in ['protein_coding', 'pre-mRNA']:
                mappability = read.get_tag('XF').split('|')[-1]
                if mappability == 'Unique':
                    biotype = biotype
                    feature = feature
                    if strandedness == 'antisense':
                        feature = feature + '_antisense'
                        biotype = 'antisense'
                    else:
                        feature = feature
                else:
                    biotype = biotype
                    feature = feature
                    if strandedness == 'antisense':
                        feature = feature + '_antisense'
                        biotype = 'antisense'
                    else:
                        feature = feature
                    # biotype = 'unannotated'
                    # biotype = 'multimapped'
                    # feature = feature + '_multimap'

            else:
                if strandedness == 'antisense':
                    biotype = 'antisense'
                    feature = feature + '_antisense'
                else:
                    feature = feature
                    
            if featuretypeD.get(feature) == None:
                featuretypeD[feature] = biotype

            read_biotype_dict[biotype] += count_value
            read_length_dict[biotype][read_length] += count_value
            read_feature_dict[feature] += count_value
            read.set_tag('CS', f"{feature}")

            return read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD

        else:

            biotype = 'unannotated'
            read_length = len(read.query_sequence)

            read_biotype_dict[biotype] += count_value
            read_length_dict[biotype][read_length] += count_value
            read.set_tag('CS', f"unannotated")

            return read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD

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
            
            if biotype in ['tRNA', 'Mt_tRNA']:

                # if biotype in ['Mt_tRNA']:
                #     biotype = 'tRNA'

                if strandedness == 'antisense':
                    biotype = 'antisense'
                    feature = feature + '_antisense'
                else:
                    feature = feature

                if read_start <= start + 10 and read_end < end - 10:
                    feature = feature + '_5p'
                elif read_start > start + 10 and read_end >= end - 10:
                    feature = feature + '_3p'
                elif read_start > start + 10 and read_end < end - 10:
                    feature = feature + '_internal'
                elif read_start <= start + 10 and read_end >= end - 10:
                    feature = feature + '_whole'

            elif biotype in ['protein_coding', 'pre-mRNA']:
                mappability = read.get_tag('XF').split('|')[-1]
                if mappability == 'Unique':
                    biotype = biotype
                    feature = feature

            else:
                if strandedness == 'antisense':
                    biotype = 'antisense'
                    feature = feature + '_antisense'
                else:
                    feature = feature
                    
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

            if read_feature_distinctseqs_dict[feature][query_sequence] is None:
                read_feature_distinctseqs_dict[feature][query_sequence] = count_value
            else:
                read_feature_distinctseqs_dict[feature][query_sequence] += count_value

            if read_biotypes_distinctseqs_dict[biotype][query_sequence] is None:
                read_biotypes_distinctseqs_dict[biotype][query_sequence] = count_value
            else:
                read_biotypes_distinctseqs_dict[biotype][query_sequence] += count_value

            return read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict

        else:

            biotype = 'unannotated'
            read_length = len(read.query_sequence)
            query_sequence = read.query_sequence

            read_biotype_dict[biotype] += count_value
            read_length_dict[biotype][read_length] += count_value
            read.set_tag('CS', f"other")

            if read_biotypes_distinctseqs_dict[biotype][query_sequence] is None:
                read_biotypes_distinctseqs_dict[biotype][query_sequence] = count_value
            else:
                read_biotypes_distinctseqs_dict[biotype][query_sequence] += count_value

            return read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict

    def parse_mapping_bed_distinct(self, sample):
        '''
        Assigns reads to bed features.
        '''

        if self.args.bampath == '':
            bampath = f'{self.args.output}/alignments/'
        else:
            if self.args.bampath.endswith('/') is False:
                bampath = self.args.bampath + '/'

        genetype_hierarchy = {
            'unannotated':0,
            'protein_coding':1,
            'pre-mRNA':1,
            'sinRNA':1,
            'scdRNA':1,
            'sutRNA':1,
            'piRNA':2,
            'lncRNA':2,
            'SINE':3,
            'LINE':3,
            'pre-tRNA':4,
            'IG_V_pseudogene': 4,
            'IG_V_gene': 5,
            'IG_C_gene': 5,
            'IG_J_gene': 5,
            'TR_V_gene': 5,
            'TR_V_pseudogene': 4,
            'TR_D_gene': 5,
            'TR_J_gene': 5,
            'TR_C_gene': 5,
            'IG_C_pseudogene': 4,
            'TR_J_pseudogene': 4,
            'IG_J_pseudogene': 4,
            'IG_D_gene': 5,
            'Y_RNA':5,
            'vault_RNA':5,
            'scRNA': 5,
            'scaRNA': 5,
            'snoRNA':5,
            'snRNA':5,
            'pre-miRNA':4,
            'star-miRNA':4,
            'Mt_rRNA':6,
            'miRNA':6,
            'rRNA':7,
            'Mt_tRNA':8,
            'tRNA':9
        }

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

        ## get the contigs
        contigs = {contig:1 for contig in bed['chrom'].unique()}

        ## get the other biotypes in the bed file
        other_biotypes = [i for i in biotypes if i not in genetype_hierarchy.keys()]

        ## add the other biotypes to the genetype_hierarchy, assume they are important so rank them higher than the other biotypes
        if len(other_biotypes) > 0:
            for i in other_biotypes:
                genetype_hierarchy[i] = 10

        ## assign a sorting order to each feature based on the genetype_hierarchy
        bed['order'] = bed['biotype'].map(genetype_hierarchy).fillna(0)

        ## sort the bed file based on the sorting order
        bed = bed.sort_values(by=['order'], ascending=False).drop(columns=['order'])

        ## group the bed file by chromosome
        bed = bed.groupby('chrom')

        bed_strand_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {'start': list(), 'end': list(), 'feature': list(), 'biotype': list(), 'strand': list()}))))
        
        ## create a dictionary to store the annotations for each chromosome so that they can be accessed quickly
        for chrom, group_df in bed:
            for strand in ['+', '-']:
                bed_strand_dict[chrom][strand] = {
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
            # self.logdata['input_bams'].append([sample, f'{bampath}{sample}.merged.bam'])
            sample_logdata = {}
            sample_logdata['sample'] = sample
            sample_logdata['input_bam'] = f'{bampath}{sample}.merged.bam'
            sample_logdata['output_bam'] = f'{bampath}{sample}.bam'
            if self.args.unmapped:
                sample_logdata['unmapped_bam'] = f'{bampath}{sample}.unmapped.bam'
            sample_logdata['mode'] = 'bed'
            sample_logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            sample_logdata['bed'] = self.args.bed
        else:
            print(f'An error occurred while trying to read {bampath}{sample}.merged.bam. Check the path and rerun the script.')
            sys.exit(1)
 
        ## delete output bam file if it already exists
        if os.path.exists(f'{bampath}{sample}.unsorted.bam'):
            os.remove(f'{bampath}{sample}.unsorted.bam')

        if self.args.unmapped:
            if os.path.exists(f'{bampath}{sample}.unmapped.bam'):
                os.remove(f'{bampath}{sample}.unmapped.bam')

        ##create new bam file to write out reads that pass filters
        try:
            output_bam = pysam.AlignmentFile(f'{bampath}{sample}.unsorted.bam', 'wb', template=bam)
        except FileNotFoundError:
            print(f'Could not write to {bampath}{sample}.unsorted.bam. Check the path and rerun the script.')
            sys.exit(1)

        if self.args.unmapped:
            try:
                unmapbam = pysam.AlignmentFile(f'{bampath}{sample}.unmapped.bam', 'wb', template=bam)
                unmap_sequences = defaultdict(lambda: 0)
            except FileNotFoundError:
                print(f'Could not write to {bampath}{sample}.unmapped.bam. Check the path and rerun the script.')
                sys.exit(1)
    
        ## check if the reads are collapsed, on a SEPARATE file handle: reading the first
        ## record of `bam` here made the main loop skip it, because fetch(until_eof=True)
        ## continues from the current file position (1 read lost per sample; fixed 2026-09-29)
        with pysam.AlignmentFile(f'{bampath}{sample}.merged.bam', 'rb') as peek:
            first_read = next(peek.fetch(until_eof=True), None)
        if first_read is not None:
            self.args.collapse = ';size=' in first_read.query_name

        ## define spike in regex
        ribospike_regex = r"SR_(\d+)_4N-T(\d+(\.\d+)?)"
        ribospike_dict = defaultdict(int)

        ## create a dictionary to store the read mappings

        read_mappability_dict = {'Mapped':0, 'Unmapped':0, 'Unique':0, 'Multi':0}
        read_mappings_dict = defaultdict(lambda: {'alignments': deque(), 'alignmentscores': deque(), 'annotations': deque(), 'biotypes': deque(), 'strand': deque(), 'barcodes': deque()})

        full_cov = defaultdict(lambda: {'depth': Counter(), 'mm': Counter(), 'biotype': None})

        read_feature_dict = defaultdict(lambda: 0)
        read_biotype_dict = defaultdict(lambda: 0)
        read_length_dict = defaultdict(lambda: {i:0 for i in range(0, 301)})
        # read_mappability_dict = defaultdict(lambda: {'Unique':0, 'Multi':0})
        read_feature_distinctseqs_dict = defaultdict(lambda: defaultdict(int))
        read_biotypes_distinctseqs_dict = defaultdict(lambda: defaultdict(int))

        with tqdm.tqdm(desc=sample, unit=' reads') as pbar:

            current_read_name = None

            for read in bam.fetch(until_eof=True):

                pbar.update(1)

                read_name = read.query_name
                read_contig = read.reference_name
                read_length = len(read.query_sequence)

                if read.is_unmapped:
                    count_value = self.get_read_count(read_name)
                    read_mappability_dict['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    read_length_dict['unmapped'][read_length] += count_value
                    output_bam.write(read)
                    if self.args.unmapped:
                        unmapbam.write(read)
                        unmap_sequences[read.query_sequence] += count_value
                    continue

                if read_length < 15:
                    count_value = self.get_read_count(read_name)
                    read_mappability_dict['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    read_length_dict['unmapped'][read_length] += count_value
                    output_bam.write(read)
                    if self.args.unmapped:
                        unmapbam.write(read)
                        unmap_sequences[read.query_sequence] += count_value
                    continue

                if '#miND' in read_contig:
                    count_value = self.get_read_count(read_name)
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
                    count_value = self.get_read_count(read_name)
                    ribospike = re.search(ribospike_regex, read_contig).group(1)
                    sequence = read.query_sequence
                    read_length = len(sequence)

                    if read_length < (int(ribospike) - 2):
                        biotype = 'RiboSpike-Short'
                        read_biotype_dict[biotype] += count_value
                        read_length_dict[biotype][read_length] += count_value
                        read_mappability_dict['Mapped'] += count_value
                        read.set_tag('CS', 'RiboSpike-fragment')
                        output_bam.write(read)
                        continue
                    else:
                        typespike = re.search(ribospike_regex, read_contig).group(2).split('.')[0]
                        typespike = 'T' + typespike
                        ribospike = f'RiboSpike-{ribospike}-{typespike}'
                        biotype = f'RiboSpike-{typespike}'
                        featuretypeD[ribospike] = biotype
                        ribospike_dict[ribospike] += count_value
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
                        self._accumulate_full_coverage(read_mappings_dict[current_read_name], mapping, self.get_read_count(current_read_name), full_cov)

                        ## check if the read is unique or multi-mapped
                        if mapping == 'Multi' and self.args.collapse == True:

                            count_values = self.partition_read_counts(best_read)
                            for cdx, count_value in enumerate(count_values):
                                if count_value == 0:
                                    continue
                                curr_best_read = best_read[cdx]
                                ## set the size= in the read name to the count value
                                oldname = curr_best_read.query_name
                                newname = oldname.split(';size=')[0] + f';size={count_value}'
                                curr_best_read.query_name = newname
                                curr_barcode = barcode[cdx]
                                curr_strandedness = strandedness[cdx]
                                # remove-at-source: same gate as the else branch, applied per partition
                                remove = (self.args.remove_multimap and mapping == 'Multi') or (self.args.remove_antisense and curr_strandedness == 'antisense')
                                if not remove:
                                    curr_best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict = self.generate_metrics_bed_distinct(curr_best_read, curr_barcode, curr_strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict)
                                    output_bam.write(curr_best_read)
                                # tally OUTSIDE the gate on purpose: removed multimap/antisense reads still
                                # count toward the mappingstats rate (design Z). Do NOT move these inside `if not remove`.
                                read_mappability_dict['Mapped'] += count_value
                                read_mappability_dict[mapping] += count_value

                        else:
                            count_value = self.get_read_count(current_read_name)
                            # remove-at-source: skip counting + BAM write for multimap/antisense when requested
                            remove = (self.args.remove_multimap and mapping == 'Multi') or (self.args.remove_antisense and strandedness == 'antisense')
                            if not remove:
                                best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict = self.generate_metrics_bed_distinct(best_read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict)
                                output_bam.write(best_read)
                            read_mappability_dict['Mapped'] += count_value
                            read_mappability_dict[mapping] += count_value
                            
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

                if contigs.get(chr) == 1:

                    pos_feature, pos_biotype, pos_strand, pos_barcode = self.find_feature(read_start, read_end, read_strand, bed_strand_dict[chr]['+']['start'], bed_strand_dict[chr]['+']['end'], bed_strand_dict[chr]['+']['biotype'], bed_strand_dict[chr]['+']['feature'], bed_strand_dict[chr]['+']['strand'], bed_strand_dict[chr]['+']['barcode'])
                    neg_feature, neg_biotype, neg_strand, neg_barcode = self.find_feature(read_start, read_end, read_strand, bed_strand_dict[chr]['-']['start'], bed_strand_dict[chr]['-']['end'], bed_strand_dict[chr]['-']['biotype'], bed_strand_dict[chr]['-']['feature'], bed_strand_dict[chr]['-']['strand'], bed_strand_dict[chr]['-']['barcode'])

                    if pos_feature == 'unannotated' and neg_feature == 'unannotated':
                        feature, biotype, strand, barcode = 'unannotated', 'unannotated', -2, 'unannotated'

                    elif pos_feature == 'unannotated' and neg_feature != 'unannotated':
                        feature, biotype, strand, barcode = neg_feature, neg_biotype, neg_strand, neg_barcode
                    
                    elif pos_feature != 'unannotated' and neg_feature == 'unannotated':
                        feature, biotype, strand, barcode = pos_feature, pos_biotype, pos_strand, pos_barcode

                    elif pos_feature != 'unannotated' and neg_feature != 'unannotated':
                        ## go with the highest strand
                        if pos_strand > neg_strand:
                            feature, biotype, strand, barcode = pos_feature, pos_biotype, pos_strand, pos_barcode
                        elif pos_strand < neg_strand:
                            feature, biotype, strand, barcode = neg_feature, neg_biotype, neg_strand, neg_barcode
                        elif pos_strand == neg_strand:
                            ## go with higher based on the genetype_hierarchy
                            if genetype_hierarchy[pos_biotype] > genetype_hierarchy[neg_biotype]:
                                feature, biotype, strand, barcode = pos_feature, pos_biotype, pos_strand, pos_barcode
                            elif genetype_hierarchy[pos_biotype] < genetype_hierarchy[neg_biotype]:
                                feature, biotype, strand, barcode = neg_feature, neg_biotype, neg_strand, neg_barcode
                            elif genetype_hierarchy[pos_biotype] == genetype_hierarchy[neg_biotype]:
                                # deterministic: pick the lexicographically smaller barcode
                                options = [(pos_feature, pos_biotype, pos_strand, pos_barcode),
                                        (neg_feature, neg_biotype, neg_strand, neg_barcode)]
                                feature, biotype, strand, barcode = min(options, key=lambda o: o[3])
                            
                else:
                    # print(read_contig)
                    feature, biotype, strand, barcode = 'unannotated', 'unannotated', -2, 'unannotated'

                read_mappings_dict[read_name]['alignments'].append(read)
                ## Get the alignment score
                try:
                    alignment_score = read.get_tag('AS') ## Get the alignment score; bowtie2
                except:
                    alignment_score = read.get_tag('NM') * -1 ## Get the alignment score; bowtie

                read_mappings_dict[read_name]['alignmentscores'].append(alignment_score)
                read_mappings_dict[read_name]['annotations'].append(feature)
                read_mappings_dict[read_name]['biotypes'].append(biotype)
                read_mappings_dict[read_name]['strand'].append(strand)
                read_mappings_dict[read_name]['barcodes'].append(barcode)

        if current_read_name is not None:

            best_read, mapping, barcode, strandedness = self.pick_bed_alignment(read_mappings_dict, current_read_name, genetype_hierarchy)
            self._accumulate_full_coverage(read_mappings_dict[current_read_name], mapping, self.get_read_count(current_read_name), full_cov)

            
            ## check if the read is unique or multi-mapped

            if mapping == 'Multi' and self.args.collapse == True:
                count_values = self.partition_read_counts(best_read)
                for cdx, count_value in enumerate(count_values):
                    if count_value == 0:
                        continue
                    curr_best_read = best_read[cdx]
                    ## set the size= in the read name to the count value
                    oldname = curr_best_read.query_name
                    newname = oldname.split(';size=')[0] + f';size={count_value}'
                    curr_best_read.query_name = newname
                    curr_barcode = barcode[cdx]
                    curr_strandedness = strandedness[cdx]
                    # remove-at-source: same gate as the else branch, applied per partition
                    remove = (self.args.remove_multimap and mapping == 'Multi') or (self.args.remove_antisense and curr_strandedness == 'antisense')
                    if not remove:
                        curr_best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict = self.generate_metrics_bed_distinct(curr_best_read, curr_barcode, curr_strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict)
                        output_bam.write(curr_best_read)
                    # tally OUTSIDE the gate on purpose: removed multimap/antisense reads still
                    # count toward the mappingstats rate (design Z). Do NOT move these inside `if not remove`.
                    read_mappability_dict['Mapped'] += count_value
                    read_mappability_dict[mapping] += count_value
            else:
                count_value = self.get_read_count(current_read_name)
                # remove-at-source: skip counting + BAM write for multimap/antisense when requested
                remove = (self.args.remove_multimap and mapping == 'Multi') or (self.args.remove_antisense and strandedness == 'antisense')
                if not remove:
                    best_read, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict = self.generate_metrics_bed_distinct(best_read, barcode, strandedness, bed_barcode_dict, read_feature_dict, read_biotype_dict, read_length_dict, featuretypeD, count_value, read_feature_distinctseqs_dict, read_biotypes_distinctseqs_dict)
                    output_bam.write(best_read)
                read_mappability_dict['Mapped'] += count_value
                read_mappability_dict[mapping] += count_value
                
            read_mappings_dict.pop(current_read_name)

        self._write_bedgraph_full(sample, full_cov)
        output_bam.close()
        bam.close()

        if self.args.unmapped:
            unmapbam.close()

        ## delete the defaultdict
        del read_mappings_dict


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

        if self.args.unmapped:
            unmap_sequences_df = pd.DataFrame.from_dict(unmap_sequences, orient='index').reset_index().rename(columns={'index':'sequence', 0:'count'})
            ## get the % of the total
            unmap_sequences_df['%'] = unmap_sequences_df['count'] / unmap_sequences_df['count'].sum() * 100
            unmap_sequences_df = unmap_sequences_df.sort_values(by='count', ascending=False).head(100)
            unmap_sequences_df.to_csv(f'{bampath}{sample}-unmappedtop100.csv')
            del unmap_sequences

        read_feature_df = pd.DataFrame(read_feature_dict, index=[0]).T.rename(columns={0:sample})

        read_biotype_df = pd.DataFrame(read_biotype_dict, index=[0]).T.rename(columns={0:sample})

        sample_logdata['reads_annotated'] = read_biotype_df[read_biotype_df.index != 'unannotated'].sum().values[0]
        sample_logdata['reads_unannotated'] = read_biotype_df[read_biotype_df.index == 'unannotated'].sum().values[0]
        sample_logdata['reads_total'] = sample_logdata['reads_annotated'] + sample_logdata['reads_unannotated']

        sample_logdata_df = pd.DataFrame.from_dict(sample_logdata, orient='index')

        ## make output directory if it does not exist (exist_ok: samples run in parallel,
        ## a separate exists-check followed by mkdir can race and raise FileExistsError)
        os.makedirs(f'{self.args.output}/metrics', exist_ok=True)

        sample_logdata_df.to_csv(f'{self.args.output}/metrics/{sample}.log', header=False)

        read_features = read_feature_df.index.tolist()
        read_feature_biotypes = [featuretypeD.get(i) for i in read_features]
        read_featuretypes_df = pd.DataFrame(read_feature_biotypes, index=read_features, columns=['feature_type']).reset_index().rename(columns={'index':'feature'})

        ## get a list of sorted biotypes by the count in biotype_df
        sorted_biotypes = read_biotype_df.sort_values(by=sample, ascending=False).index.tolist()

        read_biotype_df = read_biotype_df.loc[sorted_biotypes]

        read_length_df = read_length_df.loc[sorted_biotypes]
        read_length_df['sample'] = sample
        read_length_df = read_length_df.reset_index().rename(columns={'index':'biotype'})
        read_length_df = read_length_df[['sample', 'biotype'] + [i for i in range(0, 151)]]

        ## Get the number of distinct sequences for each feature by finding the length of the dictionary
        # read_feature_distinctseqs_df = pd.DataFrame.from_dict(read_feature_distinctseqs_dict, orient='index').T.fillna(0).astype(int)
        # read_biotypes_distinctseqs_df = pd.DataFrame.from_dict(read_biotypes_distinctseqs_dict, orient='index').T.fillna(0).astype(int)
        ## loop through the dictionary, find the length of the dictionary to make a dataframe of feature,length_of_distinct_sequences
        read_feature_distinctseqs_df = pd.DataFrame({k: [len(v)] for k, v in read_feature_distinctseqs_dict.items()}).T.rename(columns={0:'distinct_sequences'}).reset_index().rename(columns={'index':'feature'})
        read_biotypes_distinctseqs_df = pd.DataFrame({k: [len(v)] for k, v in read_biotypes_distinctseqs_dict.items()}).T.rename(columns={0:'distinct_sequences'}).reset_index().rename(columns={'index':'biotype'})
        read_biotypes_distinctseqs_df['sample'] = sample
        read_feature_distinctseqs_df['sample'] = sample
        read_feature_distinctseqs_df = read_feature_distinctseqs_df.pivot(index='feature', columns='sample', values='distinct_sequences').fillna(0).astype(int)
        read_biotypes_distinctseqs_df = read_biotypes_distinctseqs_df.pivot(index='biotype', columns='sample', values='distinct_sequences').fillna(0).astype(int)

        del read_feature_dict
        del read_biotype_dict
        del read_length_dict
        del read_feature_distinctseqs_dict
        del read_biotypes_distinctseqs_dict

        # for k,v in ribospike_end_dict.items():
        #     df = pd.DataFrame.from_dict(v, orient='index').T.fillna(0).astype(int)
        #     df = df[['5p_1', '5p_2', '5p_3', '5p_4', '3p_4', '3p_3', '3p_2', '3p_1']]
        #     df = df.div(df.sum(axis=0), axis=1) * 100
        #     df = df.stack().reset_index().rename(columns={'level_0':'nucleotide', 'level_1':'position', 0:'count'})
        #     df['sample'] = sample
        #     df['spike'] = k

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
        feature  start  end  depth  multimap_depth  strand  type
        (run-length-encoded on constant (depth, multimap_depth) blocks)."""
        import os
        outdir = f"{self.args.output}/bedgraph_full"
        os.makedirs(outdir, exist_ok=True)
        with open(f"{outdir}/{sample}.bedgraph", 'w') as f:
            # sorted: full_cov is filled from a set, whose order changes between runs
            # (string hash randomisation); sorting makes the file byte-identical (fixed 2026-09-29)
            for (feat, strand_label), entry in sorted(full_cov.items(), key=lambda kv: kv[0]):
                depth, mm, biotype = entry['depth'], entry['mm'], entry['biotype']
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
            # os.system('rm ' + self.args.o + '/bams/' + sample + '.merged.bam')
        
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
    ap.add_argument('-u', '--unmapped', required=False, help='Whether or not to dump the unmapped reads to a bam file.', action='store_true')
    # ap.add_argument('-x','--best', required=False, help='Whether or not to assign the best alignment to the read. Default=False', action='store_true')
    ap.add_argument('-b','--bed', required=False, help='Bed file to assign reads to. Default=None', default=None)
    ap.add_argument('-y','--bampath', required=False, help='Path to the bam files', default='')
    ap.add_argument('-z','--clipping', required=False, help='Number of bases to clip from the start and end of the read. Default=2', default=2)    
    ap.add_argument('-c','--collapse', required=False, help='Whether or not the read input is collapsed from vsearch ', action='store_true')
    ap.add_argument('-q','--distinct', required=False, help='Whether or not to count the number of unique sequences for each feature & biotype', default=False, action='store_true')
    ap.add_argument('--remove-multimap', required=False, default=False, action='store_true', help='Drop multimapping reads from counts, featuretypes and the parsed BAM (still tallied in mappingstats).')
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

    #Collapse reads
    samplelist = pm.samplelist

    if args.bampath == '':
        bampath = f'{args.output}/alignments/'
    else:
        if args.bampath.endswith('/') is False:
            bampath = args.bampath + '/'

    for sample in samplelist:
        pm.logdata['input_bams'][sample] = f'{bampath}{sample}.merged.bam'


    if args.distinct == True:
        if os.path.exists(args.bed) is False:
            print(f'{args.bed} does not exist. Check the path and rerun the script.')
            sys.exit(1)
        else:
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
        feature_df = pd.concat(feature_data, axis=1).fillna(0).astype(int).sort_index()

        ## get the RiboSpike reads
        ribospike_df = feature_df.filter(like='RiboSpike', axis=0)
        if ribospike_df.empty is False:
            ribospike_df.to_csv(f'{args.output}/{args.output}-ribospike.csv')

        biotype_df = pd.concat(biotype_data, axis=1).fillna(0).astype(int)

        length_df = pd.concat(length_data, axis=0).fillna(0)

        featuretypes = pd.concat(featuretype_data, axis=0).drop_duplicates(subset=['feature']).sort_values(by='feature')
        featuretypes = add_gene_length(featuretypes, args.bed)

        mappability_df = pd.concat(mapdata, axis=1).fillna(0).T
        mappability_df['Total'] = mappability_df['Mapped'] + mappability_df['Unmapped']
        

        feature_distinctseqs_df = pd.concat(featuredistseqs_data, axis=1).fillna(0).astype(int)
        biotype_distinctseqs_df = pd.concat(biotypedistseqs_data, axis=1).fillna(0).astype(int)
        
        print(f'-'*50)
        feature_df.to_csv(f'{args.output}/{args.output}-readcounts.csv')
        biotype_df.to_csv(f'{args.output}/{args.output}-typecounts.csv')
        length_df.to_csv(f'{args.output}/{args.output}-typelengths.csv', index=False)
        featuretypes.to_csv(f'{args.output}/{args.output}-featuretypes.csv', index=False)
        mappability_df.to_csv(f'{args.output}/{args.output}-mappingstats.csv')

        ## write the distinct sequences to a file
        feature_distinctseqs_df.to_csv(f'{args.output}/{args.output}-featuredistinctseqs.csv')
        biotype_distinctseqs_df.to_csv(f'{args.output}/{args.output}-biotypedistinctseqs.csv')

        print(f'Read counts written to {args.output}/{args.output}-readcounts.csv')
        print(f'Feature types written to {args.output}/{args.output}-featuretypes.csv')
        print(f'Feature counts written to {args.output}/{args.output}-typecounts.csv')
        print(f'Feature lengths written to {args.output}/{args.output}-typelengths.csv')
        print(f'Mapping stats written to {args.output}/{args.output}-mappingstats.csv')
        print(f'-'*50)

    else:
        if os.path.exists(args.bed) is False:
            print(f'{args.bed} does not exist. Check the path and rerun the script.')
            sys.exit(1)
        else:
            feature_data = []
            featuretype_data = []
            biotype_data = []
            length_data = []
            mapdata = []
            with Pool(int(args.threads)) as p:

                sampledata = p.map(pm.parse_mapping_bed, samplelist)
                for data in sampledata:
                    feature_data.append(data[0])
                    featuretype_data.append(data[1])
                    biotype_data.append(data[2])
                    length_data.append(data[3])
                    mapdata.append(data[4])

        ## write the feature stats to a file
        feature_df = pd.concat(feature_data, axis=1).fillna(0).astype(int).sort_index()

        ## get the RiboSpike reads
        ribospike_df = feature_df.filter(like='RiboSpike', axis=0)
        if ribospike_df.empty is False:
            ribospike_df.to_csv(f'{args.output}/{args.output}-ribospike.csv')

        biotype_df = pd.concat(biotype_data, axis=1).fillna(0).astype(int)

        length_df = pd.concat(length_data, axis=0).fillna(0)

        featuretypes = pd.concat(featuretype_data, axis=0).drop_duplicates(subset=['feature']).sort_values(by='feature')
        featuretypes = add_gene_length(featuretypes, args.bed)

        mappability_df = pd.concat(mapdata, axis=1).fillna(0).T
        mappability_df['Total'] = mappability_df['Mapped'] + mappability_df['Unmapped']

        print(f'-'*50)
        feature_df.to_csv(f'{args.output}/{args.output}-readcounts.csv')
        biotype_df.to_csv(f'{args.output}/{args.output}-typecounts.csv')
        length_df.to_csv(f'{args.output}/{args.output}-typelengths.csv', index=False)
        featuretypes.to_csv(f'{args.output}/{args.output}-featuretypes.csv', index=False)
        mappability_df.to_csv(f'{args.output}/{args.output}-mappingstats.csv')

        print(f'Read counts written to {args.output}/{args.output}-readcounts.csv')
        print(f'Feature types written to {args.output}/{args.output}-featuretypes.csv')
        print(f'Feature counts written to {args.output}/{args.output}-typecounts.csv')
        print(f'Feature lengths written to {args.output}/{args.output}-typelengths.csv')
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
