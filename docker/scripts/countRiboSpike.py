#!/usr/bin/env python3
#Aidan Manning 1.28.25
#Script to count the number of ribospike reads in a bam file
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import warnings
from multiprocessing import Pool
import pysam
from collections import defaultdict, deque
import tqdm
import numpy as np
import re
warnings.filterwarnings("ignore")

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
        self.sampleDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['group']))
        self.fastqDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['fastq']))

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

    def get_read_count(self, read_name):
        '''
        This function gets the read count from the read name if the reads are collapsed.
        '''
        if self.args.collapse == True:
            read_count = int(read_name.split(';size=')[1])
        elif self.args.collapse == False:
            read_count = int(1)

        return read_count

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
        
        ## Check whether or not to dump the unmapped reads to a bam file
        if self.args.unmapped:
            unmapped_bam = f'{self.args.output}/alignments/{sample}.unmapped.bam'
            unmapped_bam = pysam.AlignmentFile(unmapped_bam, 'wb', template=bam)

        ## define spike in regex
        ribospike_regex = r"SR_(\d+)_4N-T(\d+(\.\d+)?)"
        ribospike_dict = defaultdict(int)
        ribospike_end_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))
        # the same, only for reads of exactly the designed length aligned over the whole spike (QC 2026-10-02)
        ribospike_end_exact_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0)))

        read_mappability = defaultdict(int)
        read_biotype_dict = defaultdict(int)
        # read_length_dict = defaultdict(lambda: {i:0 for i in range(0, 71)})
        read_length_dict = {rs: {i:0 for i in range(0, 71)} for rs in ribospikes}
        unmap_sequences = defaultdict(int)

        ## check if the reads are collapsed, on a SEPARATE file handle: reading the first
        ## record of `bam` here made the main loop skip it, because fetch(until_eof=True)
        ## continues from the current file position (1 read lost per sample; fixed 2026-09-29)
        with pysam.AlignmentFile(bam_file, 'rb') as peek:
            first_read = next(peek.fetch(until_eof=True), None)
        if first_read is not None:
            self.args.collapse = ';size=' in first_read.query_name

        with tqdm.tqdm(desc=sample, unit=' reads') as pbar:

            for read in bam.fetch(until_eof=True):

                pbar.update(1)

                read_name = read.query_name
                read_contig = read.reference_name

                if read.is_unmapped:
                    count_value = self.get_read_count(read_name)
                    read_mappability['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    if self.args.unmapped:
                        unmapped_bam.write(read)
                        unmap_sequences[read.query_sequence] += count_value
                    continue

                if read.query_length < 15:
                    count_value = self.get_read_count(read_name)
                    read_mappability['Unmapped'] += count_value
                    read_biotype_dict['unmapped'] += count_value
                    if self.args.unmapped:
                        unmapped_bam.write(read)
                        unmap_sequences[read.query_sequence] += count_value
                    continue

                ## check for spike in reads
                if re.search(ribospike_regex, read_contig):
                    count_value = self.get_read_count(read_name)
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
                count_value = self.get_read_count(read_name)
                read_mappability['Mapped'] += count_value
                read_biotype_dict['Natural'] += count_value

        ## Close the bam files
        bam.close()

        read_mappability_df = pd.DataFrame(list(read_mappability.items()), columns=['', sample]).set_index('').T

        ## reindex, not select: a sample without unmapped (or without natural mapped) reads has no such key
        ## (found 2026-10-05 by the synthetic sample S3, whose reads all map)
        read_mappability_df = read_mappability_df.reindex(columns=['Mapped', 'Unmapped'], fill_value=0)

        read_biotype_df = pd.DataFrame(list(read_biotype_dict.items()), columns=['read_origin', sample]).set_index('read_origin')

        read_ribospike_df = pd.DataFrame(list(ribospike_dict.items()), columns=['spike', sample]).set_index('spike')

        if self.args.unmapped:
            unmapped_df = pd.DataFrame(list(unmap_sequences.items()), columns=['sequence', sample]).set_index('sequence').T
            unmapped_df = unmapped_df.T.sort_values(sample, ascending=False, kind='stable')

        ## flatten the end dictionaries
        end_df = self.ends_table(sample, ribospike_end_dict)
        end_exact_df = self.ends_table(sample, ribospike_end_exact_dict)

        read_length_df = pd.DataFrame.from_dict(read_length_dict, orient='index').T.fillna(0).astype(int)
        ## sort the columns
        read_length_df = read_length_df.reindex(sorted(read_length_df.columns), axis=1).T
        read_length_df['sample'] = sample
        read_length_df = read_length_df.reset_index().rename(columns={'index':'spike_type'})
        read_length_df = read_length_df[['sample', 'spike_type'] + [i for i in range(0, 71)]]

        return [read_mappability_df, read_biotype_df, read_ribospike_df, end_df, read_length_df, end_exact_df]

    @staticmethod
    def ends_table(sample, end_dict):
        '''One row per spike and end position with the base counts.'''
        enddata = []
        for spikelength, positions in end_dict.items():
            for end, counts in positions.items():
                enddata.append([sample, spikelength, end, counts['A'], counts['C'], counts['G'], counts['T'], counts['N']])
        return pd.DataFrame(enddata, columns=['sample', 'spike', 'end', 'A', 'C', 'G', 'T', 'N'])

    def plot_ribospike_ends(self):
        '''
        Plots a grouped stacked bar plot of the 5' and 3' ends of the spike in reads.
        '''

        ribospike_lengths = [20, 30, 40, 50, 60]
        ribospike_types = ['T1', 'T2', 'T3', 'T4']

        nucleotide_colors = {'A':'#7CAE00', 
                             'T':'#F8766D', 
                             'C':'#619CFF', 
                             'G':'#CD9600', 
                             'N':'grey', 
                             '-':'grey'}

        enddata = pd.read_csv(f'{self.args.output}/ribospike/{self.args.output}-ribospikeends.csv')
        enddata['spikelength'] = enddata['spike'].str.split('-').str[1].astype(int)
        enddata['spiketype'] = enddata['spike'].str.split('-').str[2]

        samplelist = enddata['sample'].drop_duplicates().tolist()

        for sample in samplelist:

            fig, ax = plt.subplots(len(ribospike_types), len(ribospike_lengths), figsize=(2 + 2*len(ribospike_lengths), 2 + 1.4*len(ribospike_types)), sharex=False, sharey=True)
            fig.subplots_adjust(hspace=1, wspace=0.25)

            for idx, spikelength in enumerate(ribospike_lengths):
                for jdx, spiketype in enumerate(ribospike_types):

                    data = enddata[(enddata['sample'] == sample) & (enddata['spikelength'] == spikelength) & (enddata['spiketype'] == spiketype)]

                    if data.shape[0] == 0:
                        continue

                    data = data[['end', 'A', 'C', 'G', 'T', 'N']]
                    data = data.set_index('end')

                    ## if the sum of all data is less than 20, skip
                    if data.sum().sum() < 40:
                        continue

                    ## convert to percentage
                    data = data.div(data.sum(axis=1), axis=0) * 100

                    data.plot(kind='bar', stacked=True, ax=ax[jdx, idx], color=['#7CAE00', '#619CFF', '#CD9600', '#F8766D', 'grey'], legend=False, width=0.8)
                    ax[jdx, idx].set_title(f'{spiketype} - {spikelength}nt', fontsize=10, fontweight='bold', pad=5)
                    ax[jdx, idx].set_ylabel('% of reads', fontsize=10)
                    ax[jdx, idx].set_xlabel('')

                    for axes in ['top', 'bottom', 'left', 'right']:
                        ax[jdx, idx].spines[axes].set_linewidth(0.2)

                ## add a legend in the middle right of the figure
                handles = [plt.Rectangle((0,0),1,1, color=nucleotide_colors[x]) for x in ['A', 'C', 'G', 'T', 'N']]
                labels = ['A', 'C', 'G', 'T', 'N']

                fig.legend(handles, labels, loc='center right', fontsize=8, title='Nucleotide', title_fontsize=8, bbox_to_anchor=(0.98, 0.5), frameon=False)

            plt.savefig(f'{self.args.output}/ribospike/{sample}-ribospikeends2.pdf', dpi=300, bbox_inches='tight')
            plt.close()

        return

    def plot_ribospike_ends_line(self):
        '''
        Plots a line plot of the 5' and 3' ends of the spike in reads.
        '''

        ribospike_lengths = [20, 30, 40, 50, 60]
        ribospike_types = ['T1', 'T2', 'T3', 'T4']

        nucleotide_colors = {'A':'#7CAE00', 
                             'T':'#F8766D', 
                             'C':'#619CFF', 
                             'G':'#CD9600', 
                             'N':'grey', 
                             '-':'grey'}

        enddata = pd.read_csv(f'{self.args.output}/ribospike/{self.args.output}-ribospikeends.csv')
        enddata['spikelength'] = enddata['spike'].str.split('-').str[1].astype(int)
        enddata['spiketype'] = enddata['spike'].str.split('-').str[2]
        enddata['group'] = enddata['sample'].map(self.sampleDict)

        grouplist = enddata['group'].drop_duplicates().tolist()

        for rtype in ribospike_types:

            fig, ax = plt.subplots(len(grouplist), len(ribospike_lengths), figsize=(2 + 2*len(ribospike_lengths), 2 + 1.8*len(grouplist)), sharex=False, sharey=True)
            fig.subplots_adjust(hspace=0.8, wspace=0.25)  

            for idx, spikelength in enumerate(ribospike_lengths):
                for jdx, group in enumerate(grouplist):

                    data = enddata[(enddata['group'] == group) & (enddata['spikelength'] == spikelength) & (enddata['spiketype'] == rtype)]

                    if data.shape[0] == 0:
                        continue

                    data = data[['end', 'A', 'C', 'G', 'T', 'N']]
                    data = data.set_index('end')

                    ## if the sum of all data is less than 20, skip
                    if data.sum().sum() < 40:
                        continue

                    ## convert to percentage
                    data = data.div(data.sum(axis=1), axis=0) * 100

                    data = data.stack().reset_index().rename(columns={0:'value', 'level_1':'nucleotide'})

                    sns.lineplot(data=data, x='end', y='value', hue='nucleotide', ax=ax[jdx, idx], palette=nucleotide_colors, linewidth=0.8, legend=False, err_style="bars", errorbar=("se", 1))
      

                    # data.plot(kind='line', ax=ax[jdx, idx], color=['#7CAE00', '#619CFF', '#CD9600', '#F8766D', 'grey'], legend=False, linewidth=0.8, )

                    ax[jdx, idx].set_title(f'{group} - {spikelength}nt', fontsize=7, pad=2)
                    ax[jdx, idx].set_ylabel('% of reads', fontsize=10)
                    ax[jdx, idx].set_xlabel('')
                    ## rotate the x-axis labels
                    ax[jdx, idx].set_xticklabels(ax[jdx, idx].get_xticklabels(), rotation=45, horizontalalignment='right')
                    ax[jdx, idx].set_yticks(np.arange(0, 101, 25))

                    ax[jdx, idx].tick_params(axis='both', which='major', labelsize=8)


                    for axes in ['top', 'bottom', 'left', 'right']:
                        ax[jdx, idx].spines[axes].set_linewidth(0.2)

                ## add a legend in the middle right of the figure
                handles = [plt.Rectangle((0,0),1,1, color=nucleotide_colors[x]) for x in ['A', 'C', 'G', 'T', 'N']]
                labels = ['A', 'C', 'G', 'T', 'N']

                fig.legend(handles, labels, loc='center right', fontsize=8, title='Nucleotide', title_fontsize=8, bbox_to_anchor=(0.98, 0.5), frameon=False)



            plt.savefig(f'{self.args.output}/ribospike/{rtype}-ribospikeends2.pdf', dpi=300, bbox_inches='tight')
            plt.close()

        return

def main():

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('-s','--samples', required=True, help='Sample metadata file. Format: sample,group,fastq')
    ap.add_argument('-@','--threads', required=False, help='Number of threads to process paraellel fastq files with. Default=4', default=4)
    ap.add_argument('-o', '--output', required=True, help='Output prefix.')
    ap.add_argument('-u', '--unmapped', required=False, help='Whether or not to dump the unmapped reads to a bam file.', action='store_true')
    ap.add_argument('-c','--collapse', required=False, help='Whether or not the read input is collapsed from vsearch ', action='store_true')
    ap.add_argument('-y','--bampath', required=False, help='Path to the bam files', default='')
    ap.add_argument('--paired', required=False, help='Whether or not the reads are paired end.', action='store_true')
    
    args = ap.parse_args()
    
    startTime = datetime.now()

        ## check if the samplefile has 3 or 4 columns so we can determine if the reads are paired-end or not
    if args.samples != 'default':
        with open(args.samples, 'r') as f:
            line = f.readline()
            if len(line.split(',')) == 4:
                args.paired = True
            else:
                args.paired = False
    else:
        print('No sample metadata file provided. Please provide a sample metadata file.')
        sys.exit(1)

    pm = CountRiboSpike(args)

    print(f'-'*50)
    print(f'-'*50)
    print(f'Counting ribospike reads...')
    print(f'-'*50)
    print(f'-'*50)

    # pm.plot_ribospike_ends()
    # pm.plot_ribospike_ends_line()
    # sys.exit()

    if args.bampath == '':
        args.bampath = f'{args.output}/alignments'

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if not os.path.exists(f'{args.output}/ribospike'):
        os.mkdir(f'{args.output}/ribospike')

    #Collapse reads
    samplelist = pm.samplelist

    mapdata = []
    biotypedata = []
    ribospike_data = []
    end_data = []
    end_exact_data = []
    length_data = []

    with Pool(int(args.threads)) as p:
        outdata = p.map(pm.parse_mappings, samplelist)

    for data in outdata:
        mapdata.append(data[0])
        biotypedata.append(data[1])
        ribospike_data.append(data[2])
        end_data.append(data[3])
        length_data.append(data[4])
        end_exact_data.append(data[5])

    mapdata = pd.concat(mapdata, axis=0).fillna(0)
    biotypedata = pd.concat(biotypedata, axis=1).fillna(0)
    ribospike_data = pd.concat(ribospike_data, axis=1).fillna(0)
    end_data = pd.concat(end_data, axis=0).fillna(0)
    end_exact_data = pd.concat(end_exact_data, axis=0).fillna(0)
    length_data = pd.concat(length_data, axis=0).fillna(0)

    # mapdata.to_csv(f'{args.output}/ribospike/{args.output}-mappability.csv')
    # file names and columns of release 0.3.0 (before: ribotype.csv with column biotype, ribolengths.csv with
    # column biotype, column ribospike in ribospike.csv and the ends tables; the value Unmapped is now unmapped)
    biotypedata.to_csv(f'{args.output}/ribospike/{args.output}-readorigin.csv')
    ribospike_data.to_csv(f'{args.output}/ribospike/{args.output}-ribospike.csv')
    end_data.to_csv(f'{args.output}/ribospike/{args.output}-ribospikeends.csv', index=False)
    end_exact_data.to_csv(f'{args.output}/ribospike/{args.output}-ribospikeends_exact.csv', index=False)
    length_data.to_csv(f'{args.output}/ribospike/{args.output}-spikelengths.csv', index=False)

    # pm.plot_ribospike_ends()
    # pm.plot_ribospike_ends_line()

    print(f'-'*50)
    print(f'-'*50)
    print('Time elasped for read parsing: ', datetime.now() - startTime)
    print(f'-'*50)
    print(f'-'*50)
    
if __name__ == "__main__":
    main()
