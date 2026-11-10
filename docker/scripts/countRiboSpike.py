#!/usr/bin/env python3
#Aidan Manning 1.28.25
#Script to count the number of ribospike reads in a bam file
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
from multiprocessing import Pool
import pysam
from collections import defaultdict
import tqdm
import re

class CountRiboSpike(object):

    def __init__(self, args):
        '''
        This class is used to parse the mappings from the various mapping tools.
        '''

        self.args = args

        if args.paired == True:
            print('Paired-end reads detected.')
            self.sampleinfo = pd.read_csv(args.samples, sep=',', names=['sample', 'group', 'fastq', 'fastq2'])
        else:
            print('Single-end reads detected.')
            self.sampleinfo = pd.read_csv(args.samples, sep=',', names=['sample', 'group', 'fastq'])

        self.samplelist = self.sampleinfo['sample'].tolist()

        return
    
    def get_ends(self, ribospike_end_dict, sequence, count_value, spikelength):
        '''
        This function gets the 4 most 5' and 3' end of the spike in reads.
        '''
        ## get the 5' end
        five_prime = sequence[:4]
        ## get the 3' end
        three_prime = sequence[-4:]
        for idx, pos in enumerate(five_prime):
            pos = f'5p_{idx+1}'
            nucleotide = five_prime[idx]
            ribospike_end_dict[spikelength][pos][nucleotide] += count_value
        for idx, pos in enumerate(three_prime):
            nucleotide = three_prime[idx]
            ## reverse the index
            idx = len(three_prime) - idx - 1
            pos = f'3p_{idx+1}'
            ribospike_end_dict[spikelength][pos][nucleotide] += count_value

        return ribospike_end_dict

    def parse_mappings(self, sample):
        '''
        Parses the mappings from the bam files.
        '''

        ## Get the bam file
        if os.path.exists(f'{self.args.bampath}/{sample}.bam'):
            bam_file = f'{self.args.bampath}/{sample}.bam'
            ## Open the bam file
            bam = pysam.AlignmentFile(bam_file, 'rb')
        else:
            print(f'{self.args.bampath}/{sample}.bam does not exist. Skipping...')
            return
        
        ribospikes = ['RiboSpike-T1', 'RiboSpike-T2', 'RiboSpike-T3', 'RiboSpike-T4']
        
        ## define spike in regex
        ribospike_regex = r"SR_(\d+)_4N-T(\d+(\.\d+)?)"
        ribospike_dict = defaultdict(int)
        ribospike_end_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
        # the same, only for reads of exactly the designed length aligned over the whole spike (QC 2026-10-02)
        ribospike_end_exact_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))

        read_biotype_dict = defaultdict(int)
        # read_length_dict = defaultdict(lambda: {i:0 for i in range(0, 71)})
        read_length_dict = {rs: {i:0 for i in range(0, 71)} for rs in ribospikes}

        with tqdm.tqdm(desc=sample, unit=' reads') as pbar:

            for read in bam.fetch(until_eof=True):

                pbar.update(1)

                read_name = read.query_name
                read_contig = read.reference_name

                if read.is_unmapped or read.query_length < 15:     # the 15-nt floor (the shortest insert trimming keeps)
                    count_value = 1
                    read_biotype_dict['unmapped'] += count_value
                    continue

                ## check for spike in reads
                if re.search(ribospike_regex, read_contig):
                    count_value = 1
                    ribospike = re.search(ribospike_regex, read_contig).group(1)
                    designed_length = int(ribospike)
                    sequence = read.query_sequence
                    read_length = len(sequence)
                    
                    if read_length < (int(ribospike) - 2):
                        biotype = 'RiboSpike-Fragment'
                        typespike = re.search(ribospike_regex, read_contig).group(2).split('.')[0]
                        typespike = 'T' + typespike
                        ribotype = f'RiboSpike-{typespike}'
                        read_biotype_dict[biotype] += count_value
                        # read_length_dict[biotype][read_length] += count_value
                        continue
                    else:
                        typespike = re.search(ribospike_regex, read_contig).group(2).split('.')[0]
                        typespike = 'T' + typespike
                        ribospike = f'RiboSpike-{ribospike}-{typespike}'
                        biotype = f'RiboSpike-{typespike}'
                        ribospike_dict[ribospike] += count_value
                        read_biotype_dict[biotype] += count_value
                        read_length_dict[biotype][read_length] += count_value

                        ## get the 5' and 3' ends of the spike in reads
                        ribospike_end_dict = self.get_ends(ribospike_end_dict, sequence, count_value, ribospike)

                        ## exact-length reads aligned from the first to the last spike base without indels:
                        ## only for them are the first and last 4 bases certainly the random 4N positions
                        ## (reads 1-2 nt short or long put fixed spike bases into these positions)
                        if (read_length == designed_length and read.reference_start == 0
                                and read.cigarstring == f'{designed_length}M'):
                            ribospike_end_exact_dict = self.get_ends(ribospike_end_exact_dict, sequence, count_value, ribospike)

                        continue



                ## else, if the read is mapped, count it
                count_value = 1
                read_biotype_dict['Natural'] += count_value

        ## Close the bam files
        bam.close()

        read_biotype_df = pd.DataFrame(list(read_biotype_dict.items()), columns=['read_origin', sample]).set_index('read_origin')

        read_ribospike_df = pd.DataFrame(list(ribospike_dict.items()), columns=['spike', sample]).set_index('spike')

        ## flatten the end dictionaries
        end_df = self.ends_table(sample, ribospike_end_dict)
        end_exact_df = self.ends_table(sample, ribospike_end_exact_dict)

        read_length_df = pd.DataFrame.from_dict(read_length_dict, orient='index').T.fillna(0).astype(int)
        ## sort the columns
        read_length_df = read_length_df.reindex(sorted(read_length_df.columns), axis=1).T
        read_length_df['sample'] = sample
        read_length_df = read_length_df.reset_index().rename(columns={'index':'spike_type'})
        read_length_df = read_length_df[['sample', 'spike_type'] + [i for i in range(0, 71)]]

        return [read_biotype_df, read_ribospike_df, end_df, read_length_df, end_exact_df]

    @staticmethod
    def ends_table(sample, end_dict):
        '''One row per spike and end position with the base counts.'''
        enddata = []
        for spikelength, positions in end_dict.items():
            for end, counts in positions.items():
                enddata.append([sample, spikelength, end, counts['A'], counts['C'], counts['G'], counts['T'], counts['N']])
        return pd.DataFrame(enddata, columns=['sample', 'spike', 'end', 'A', 'C', 'G', 'T', 'N'])

def main():

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('-s','--samples', required=True, help='Sample metadata file. Format: sample,group,fastq')
    ap.add_argument('-@','--threads', required=False, help='Number of threads to process paraellel fastq files with. Default=4', default=4)
    ap.add_argument('-o', '--output', required=True, help='Output prefix.')
    ap.add_argument('-y','--bampath', required=False, help='Path to the bam files', default='')
    
    args = ap.parse_args()
    
    startTime = datetime.now()

    with open(args.samples) as f:                # 4 columns in the sample sheet = paired-end
        args.paired = len(f.readline().split(',')) == 4

    pm = CountRiboSpike(args)

    print(f'-'*50)
    print(f'-'*50)
    print(f'Counting ribospike reads...')
    print(f'-'*50)
    print(f'-'*50)

    if args.bampath == '':
        args.bampath = f'{args.output}/alignments'

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(f'{args.output}/ribospike'):
        os.mkdir(f'{args.output}/ribospike')

    #Collapse reads
    samplelist = pm.samplelist

    biotypedata = []
    ribospike_data = []
    end_data = []
    end_exact_data = []
    length_data = []

    with Pool(int(args.threads)) as p:
        outdata = p.map(pm.parse_mappings, samplelist)

    for data in outdata:
        biotypedata.append(data[0])
        ribospike_data.append(data[1])
        end_data.append(data[2])
        length_data.append(data[3])
        end_exact_data.append(data[4])

    biotypedata = pd.concat(biotypedata, axis=1).fillna(0)
    ribospike_data = pd.concat(ribospike_data, axis=1).fillna(0)
    end_data = pd.concat(end_data, axis=0).fillna(0)
    end_exact_data = pd.concat(end_exact_data, axis=0).fillna(0)
    length_data = pd.concat(length_data, axis=0).fillna(0)

    # file names and columns of release 0.3.0 (before: ribotype.csv with column biotype, ribolengths.csv with
    # column biotype, column ribospike in ribospike.csv and the ends tables; the value Unmapped is now unmapped)
    biotypedata.to_csv(f'{args.output}/ribospike/{args.output}-readorigin.csv')
    ribospike_data.to_csv(f'{args.output}/ribospike/{args.output}-ribospike.csv')
    end_data.to_csv(f'{args.output}/ribospike/{args.output}-ribospikeends.csv', index=False)
    end_exact_data.to_csv(f'{args.output}/ribospike/{args.output}-ribospikeends_exact.csv', index=False)
    length_data.to_csv(f'{args.output}/ribospike/{args.output}-spikelengths.csv', index=False)

    print(f'-'*50)
    print(f'-'*50)
    print('Time elasped for read parsing: ', datetime.now() - startTime)
    print(f'-'*50)
    print(f'-'*50)
    
if __name__ == "__main__":
    main()
