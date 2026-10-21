#!/usr/bin/env python3
# Aidan Manning 1.28.25
# Script to trim adapters from small RNA reads.

import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import warnings
import gzip
from collections import defaultdict
from toolRun import run  # stops on a failed tool
from Bio import SeqIO
import tqdm
import numpy as np
warnings.filterwarnings("ignore")
import json

class TrimReads():
    
    def __init__(self, args):
        '''
        This class is used to handle small RNA read processing and adapter trimming.
        '''

        self.args = args

        # Using a regular dictionary for logdata (avoiding lambda functions)
        self.logdata = {}
        self.logdata['module'] = 'trimReads'
        self.logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def generate_sampleinfo(self):
        '''
        Generates a list of samples from the sample metadata file if provided. Else, uses all fastq files in the current directory.
        '''

        if self.args.samples == 'default':

            if self.args.paired == True:
                print(f'No sample metadata file provided. Does not work with paired-end reads.')
                sys.exit(1)

            print(f'No sample metadata file provided. Using all fastq files in the current directory.')

            ## get all fastq files in the current directory
            fastqlist = [f for f in os.listdir() if f.endswith('.fastq') or f.endswith('.fq') or f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
            ## sorts the list
            fastqlist = sorted(fastqlist)
            samplelist = [f.split('.')[0] for f in fastqlist]

            self.sampleinfo = {}

            for idx, sample in enumerate(samplelist):

                self.sampleinfo[sample] = {'group': sample, 'fastq': fastqlist[idx]}

        else:

            if self.args.paired == False:
                print(f'Using sample metadata file: {self.args.samples}')
                sampledf = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq'])
                sampledf = sampledf.set_index('sample')

                ## check if the fastqdir argument is provided
                # os.path.join keeps an absolute name as it is (merged lanes)
                if self.args.fastqdir != 'default':
                    sampledf['fastq'] = sampledf['fastq'].apply(lambda x: os.path.join(self.args.fastqdir, x))

            elif self.args.paired == True:
                print(f'Using sample metadata file: {self.args.samples}')
                sampledf = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq', 'fastq_2'])
                sampledf = sampledf.set_index('sample')

                ## check if the fastqdir argument is provided
                if self.args.fastqdir != 'default':
                    sampledf['fastq'] = sampledf['fastq'].apply(lambda x: os.path.join(self.args.fastqdir, x))
                    sampledf['fastq_2'] = sampledf['fastq_2'].apply(lambda x: os.path.join(self.args.fastqdir, x))

            self.sampleinfo = sampledf.to_dict(orient='index')

        self.logdata['samples'] = sorted(self.sampleinfo.keys())
        self.logdata['threads'] = self.args.threads
        self.logdata['output'] = self.args.output
        self.logdata['tools'] = []
        self.logdata['commands'] = []

        print(f'Samples: {sorted(self.sampleinfo.keys())}')
        print(f'-'*50)

        return

    def feed_trimming(self):
        '''
        Feeds the samples to the appropriate trimming function based on the library type.
        '''

        samples = self.sampleinfo.keys()
        samples = sorted(samples)

        
        for sample in samples:

            fastq = self.sampleinfo[sample]['fastq']

            print(f'Trimming reads for {sample} | {fastq} | {self.args.library} | {self.args.tool}')

            if self.args.library == 'RSB':

                if self.args.paired == True:
                    fastq2 = self.sampleinfo[sample]['fastq_2']
                    self.trim_reads_rsb_paired(sample, fastq, fastq2, self.args.tool)

                elif self.args.paired == False:

                    self.trim_reads_rsb(sample, fastq, self.args.tool)
                    self.fastqc(sample, 'raw')
                    self.fastqc(sample, 'trimmed')
                    
            elif self.args.library == 'SMART':

                if self.args.paired == True:
                    self.trim_reads_smart_paired(sample, fastq, self.args.tool)

                elif self.args.paired == False:
                    print('Paired-end reads needs to be implemented.')
                    sys.exit(1)
                    self.trim_reads_smart(sample, fastq, self.args.tool)

            else:
                print(f'Library type not recognized. Please specify either RSB or SMART.')
                sys.exit(1)

            print('-' * 50)

        if self.args.library == 'RSB':
            if self.args.paired == False:
                self.parse_logfiles_rsb()
            elif self.args.paired == True:
                self.parse_logfiles_rsb_paired()
        # sys.exit()

        outdata = pd.DataFrame.from_dict(self.sampleinfo, orient='index')
        outdata.to_csv(f'{self.args.output}/{self.args.output}-trimstats.csv')
        print('-' * 50)
        print(f'Trimming results saved to {self.args.output}/{self.args.output}-trimstats.csv')
        print('-' * 50)
        return

    def trim_reads_rsb(self, sample, fastq, tool):
        '''
        Trims adapters from fastq files. RealSeq library.
        '''

        threads = self.args.threads
        trimmeddir = f'{self.args.output}/trimmed'
        trimlogdir = f'{self.args.output}/trimmed/logs'
        dimer_fq = f'{trimmeddir}/{sample}.dimer.fq.gz'
        output_fq = f'{trimmeddir}/{sample}.trimmed.fq.gz'

        if not os.path.exists(trimmeddir):
            os.mkdir(trimmeddir)
        if not os.path.exists(trimlogdir):
            os.mkdir(trimlogdir)

        # Skip trimming if output already exists
        if os.path.exists(output_fq):
            print(f'{output_fq} already exists. Skipping...')
            self.sampleinfo[sample]['trimmed_fastq'] = output_fq.replace(f'{self.args.output}/', '')
            if len(self.logdata['tools']) == 0:
                self.logdata['tools'].append(tool)
                self.logdata['commands'].append(f'cutadapt -O 1 -j {threads} -u 1 -m 15 -Z -a TGGAATTCTCGGGTGCCAAGG')
            return
        if 'NEB' in sample:

            if tool == 'cutadapt':
                cmd = f'cutadapt -O 1 -j {threads} -m 1 -Z -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {dimer_fq}  {fastq} > {trimlogdir}/{sample}.dimerlog'
                # cmd = f'cutadapt -O 1 -j {threads} -u 1 -m 1 -Z -a TGGAATTCTCGGGTGCCAAGG -o {dimer_fq} {fastq} > {trimlogdir}/{sample}.dimerlog'
                run(cmd)
                cmd = f'cutadapt -O 1 -j {threads} -m 15 -Z -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {output_fq} {fastq} > {trimlogdir}/{sample}.trimlog'
                # cmd = f'cutadapt -O 1 -j {threads} -u 1 -m 15 -Z -a TGGAATTCTCGGGTGCCAAGG -o {output_fq} {fastq} > {trimlogdir}/{sample}.trimlog'
                run(cmd)
        else:
            if tool == 'cutadapt':
                # cmd = f'cutadapt -O 1 -j {threads} -m 1 -Z -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {dimer_fq}  {fastq} > {self.args.output}/trimming/{sample}.dimerlog'
                cmd = f'cutadapt -O 1 -j {threads} -u 1 -m 1 -Z -a TGGAATTCTCGGGTGCCAAGG --json {trimlogdir}/{sample}.dimer.cutadapt.json -o {dimer_fq} {fastq} > {trimlogdir}/{sample}.dimerlog'
                run(cmd)
                # cmd = f'cutadapt -O 1 -j {threads} -m 15 -Z -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {output_fq} {fastq} > {self.args.output}/trimming/{sample}.trimlog'
                cmd = f'cutadapt -O 1 -j {threads} -u 1 -m 15 -Z -a TGGAATTCTCGGGTGCCAAGG --json {trimlogdir}/{sample}.trim.cutadapt.json -o {output_fq} {fastq} > {trimlogdir}/{sample}.trimlog'
                run(cmd)

            elif tool == 'fastp':

                cmd = f'fastp -i {fastq} \
                -o {dimer_fq} \
                -j {self.args.output}/trimming/{sample}.dimer.json \
                -w {threads} \
                --adapter_sequence=TGGAATTCTCGGGTGCCAAGG \
                -l 1 \
                --trim_front1 1 \
                --trim_poly_g \
                2> {self.args.output}/trimming/{sample}.dimerlog'

                run(cmd)

                cmd = f'fastp -i {fastq} \
                -o {output_fq} \
                -j {self.args.output}/trimming/{sample}.trim.json \
                -w {threads} \
                --adapter_sequence=TGGAATTCTCGGGTGCCAAGG \
                -l 15 \
                --trim_front1 1 \
                --trim_poly_g \
                2> {self.args.output}/trimming/{sample}.trimlog'

                run(cmd)

        ## Remove dimer file
        run(f'rm -f {dimer_fq}')

        self.sampleinfo[sample]['trimmed_fastq'] = output_fq.replace(f'{self.args.output}/', '')

        if len(self.logdata['tools']) == 0:
            self.logdata['tools'].append(tool)
            self.logdata['commands'].append(f'cutadapt -O 1 -j {threads} -u 1 -m 15 -Z -a TGGAATTCTCGGGTGCCAAGG')

        print(f'Finished trimming reads for {sample} | {output_fq}')
        print('-' * 50)

        return
    
    def trim_reads_rsb_paired(self, sample, fastq, fastq2, tool):
        '''
        Trims adapters from fastq files. RealSeq library using paired-end reads.
        '''

        print(f'Trimming paired end RSB reads...')

        threads = self.args.threads
        trimmeddir = f'{self.args.output}/trimmed'
        trimlogdir = f'{self.args.output}/trimmed/logs'
        dimer_fq_r1 = f'{trimmeddir}/{sample}.R1.dimer.fq.gz'
        output_fq_r1 = f'{trimmeddir}/{sample}.R1.trimmed.fq.gz'
        dimer_fq_r2 = f'{trimmeddir}/{sample}.R2.dimer.fq.gz'
        output_fq_r2 = f'{trimmeddir}/{sample}.R2.trimmed.fq.gz'

        if not os.path.exists(trimmeddir):
            os.mkdir(trimmeddir)
        if not os.path.exists(trimlogdir):
            os.mkdir(trimlogdir)

        # Skip trimming if output already exists
        if os.path.exists(f'{trimmeddir}/{sample}.merged.fq.gz'):
            print(f'{output_fq_r1} already exists. Skipping...')
            self.sampleinfo[sample]['trimmed_fastq'] = output_fq_r1.replace(f'{self.args.output}/', '')
            return
        
        # RiboMarker inserts carry ONE extra base at their 5' end (protocol, Appendix D), i.e. at the
        # start of R1 -> '-u 1'. R2 starts at the insert's 3' end, which has no extra base, so R2 is
        # not clipped (the former '-U 1' removed the last real base of every insert; fixed 2026-09-29).
        if tool == 'cutadapt':
            
            cmd = f'cutadapt -O 1 -j {threads} -u 1 -m 1 -Z -a TGGAATTCTCGGGTGCCAAGG -A TGATCGTCGGACT --json {trimlogdir}/{sample}.dimer.cutadapt.json -o {dimer_fq_r1} -p {dimer_fq_r2} {fastq} {fastq2} > {trimlogdir}/{sample}.dimerlog'
            run(cmd)

            cmd = f'cutadapt -O 1 -j {threads} -u 1 -m 15 -Z -a TGGAATTCTCGGGTGCCAAGG -A TGATCGTCGGACT --json {trimlogdir}/{sample}.trim.cutadapt.json -o {output_fq_r1} -p {output_fq_r2} {fastq} {fastq2} > {trimlogdir}/{sample}.trimlog'
            run(cmd)

            run(f'fastp -i {output_fq_r1} -I {output_fq_r2} --merged_out {trimmeddir}/{sample}.merged.fq.gz --merge --thread {threads} --overlap_len_require 10 --overlap_diff_limit 2 --disable_adapter_trimming -j {trimlogdir}/{sample}.fastp.json -h {trimlogdir}/{sample}.fastp.html')
            print(f'Merged {output_fq_r1} and {output_fq_r2} into {trimmeddir}/{sample}.merged.fq.gz')

        # elif tool == 'fastp':

        #     cmd = f'fastp -i {fastq} \
        #     -o {dimer_fq} \
        #     -j {self.args.output}/trimming/{sample}.dimer.json \
        #     -w {threads} \
        #     --adapter_sequence=TGGAATTCTCGGGTGCCAAGG \
        #     -l 1 \
        #     --trim_front1 1 \
        #     --trim_poly_g \
        #     2> {self.args.output}/trimming/{sample}.dimerlog'

        #     os.system(cmd)

        #     cmd = f'fastp -i {fastq} \
        #     -o {output_fq} \
        #     -j {self.args.output}/trimming/{sample}.trim.json \
        #     -w {threads} \
        #     --adapter_sequence=TGGAATTCTCGGGTGCCAAGG \
        #     -l 15 \
        #     --trim_front1 1 \
        #     --trim_poly_g \
        #     2> {self.args.output}/trimming/{sample}.trimlog'

        #     os.system(cmd)

        ## Remove dimer file
        run(f'rm -f {dimer_fq_r1}')
        run(f'rm -f {dimer_fq_r2}')

        self.sampleinfo[sample]['trimmed_fastq'] = output_fq_r1.replace(f'{self.args.output}/', '')
        self.sampleinfo[sample]['trimmed_fastq2'] = output_fq_r2.replace(f'{self.args.output}/', '')
        self.sampleinfo[sample]['merged_fastq'] = f'{trimmeddir}/{sample}.merged.fq.gz'

        print(f'Finished trimming reads for {sample} | {output_fq_r1} | {output_fq_r2}')
        print('-' * 50)

        return

    def trim_reads_smart_paired(self, sample, fastq, tool):
        '''
        Trims adapters from fastq files. SMART library using paired-end reads.
        '''
        
        threads = self.args.threads
        r1_trimmed_fq = f'{self.args.output}/{sample}.trimmed.R1.fq.gz'
        r2_trimmed_fq = f'{self.args.output}/{sample}.trimmed.R2.fq.gz'
        r1_tmp_fq = f'{self.args.output}/{sample}.tmp.R1.fq.gz'
        r2_tmp_fq = f'{self.args.output}/{sample}.tmp.R2.fq.gz'

        ## check if the trimmed fastq files already exist
        if os.path.isfile(r1_trimmed_fq) and os.path.isfile(r2_trimmed_fq):
            print(f'Trimmed fastq files for {sample} already exist. Skipping trimming...')
            return
        
        print(f'Removing UMI and dedup for: {sample}')
        run(f'fastp \
                --dedup \
                -i {fastq} \
                -o {r1_tmp_fq} \
                -I {fastq.replace("_R1", "_R2")} \
                -O {r2_tmp_fq} \
                -j {self.args.output}/logs/{sample}.umilog.json \
                -w {threads} \
                --umi --umi_loc read2 --umi_len 8 --trim_front2=6')
        
        print(f'Trimming sample: {sample}')
        run(f'fastp \
                --adapter_sequence AGATCGGAAGAGC \
                --adapter_sequence_r2 AGATCGGAAGAGC \
                -i {r1_tmp_fq} \
                -o {r1_trimmed_fq} \
                -I {r2_tmp_fq} \
                -O {r2_trimmed_fq} \
                -j {self.args.output}/logs/{sample}.trimlog.json \
                -w {threads}')
        
        print(f'Trimmed R1 file: {r1_trimmed_fq}')
        print(f'Trimmed R2 file: {r2_trimmed_fq}')

        return

    def parse_logfiles_rsb(self):
        '''
        Parses log files generated during adapter trimming.
        '''

        logpath = f'{self.args.output}/trimmed/logs'

        for sample in self.sampleinfo.keys():
            if os.path.isfile(f'{logpath}/{sample}.dimerlog'):
                for line in open(f'{logpath}/{sample}.dimerlog'):
                    if 'Total reads processed' in line:
                        total_reads = int(line.split(' ')[-1].rstrip().replace(',', ''))
                    elif 'Reads that were too short' in line:
                        dimer_reads = int(line.split(' ')[-2].rstrip().replace(',', ''))
                        dimer_percent = line.split(' ')[-1].rstrip().replace('(', '').replace(')', '')
                        self.sampleinfo[sample]['total_reads'] = total_reads
                        self.sampleinfo[sample]['dimer_reads'] = dimer_reads
                        self.sampleinfo[sample]['dimer_percent'] = dimer_percent
                        # dimer_data.append([sample, total_reads, dimer_reads, dimer_percent])
                        break
                    else:
                        pass
            else:
                print(f'No dimer log file found for {sample}. Skipping...')
                continue

        for sample in self.sampleinfo.keys():
            if os.path.isfile(f'{logpath}/{sample}.trimlog'):
                for line in open(f'{logpath}/{sample}.trimlog'):
                    if 'Reads that were too short' in line:
                        under15_reads = int(line.split(' ')[-2].rstrip().replace(',', ''))
                        under15_percent = line.split(' ')[-1].rstrip().replace('(', '').replace(')', '')
                        self.sampleinfo[sample]['<15nt_reads'] = under15_reads
                        self.sampleinfo[sample]['<15nt_percent'] = under15_percent
                    if 'Reads written' in line:
                        trimmed_reads = int(line.split(' ')[-2].rstrip().replace(',', ''))
                        trimmed_percent = line.split(' ')[-1].rstrip().replace('(', '').replace(')', '')
                        self.sampleinfo[sample]['trimmed_reads'] = trimmed_reads
                        self.sampleinfo[sample]['trimmed_percent'] = trimmed_percent
                        # trim_data.append([sample, trimmed_reads, trimmed_percent])
                        break
                    else:
                        pass
            else:
                print(f'No trim log file found for {sample}. Skipping...')
                continue

        return

    def parse_logfiles_rsb_paired(self):
        '''
        Parses log files generated during adapter trimming.
        '''

        logpath = f'{self.args.output}/trimmed/logs'

        for sample in self.sampleinfo.keys():
            if os.path.isfile(f'{logpath}/{sample}.trimlog'):
                for line in open(f'{logpath}/{sample}.trimlog'):
                    if 'Pairs that were too short' in line:
                        under15_reads = int(line.split(' ')[-2].rstrip().replace(',', ''))
                        under15_percent = line.split(' ')[-1].rstrip().replace('(', '').replace(')', '')
                        self.sampleinfo[sample]['dimer_reads'] = under15_reads
                        self.sampleinfo[sample]['dimer_percent'] = under15_percent
                    if 'Pairs written' in line:
                        trimmed_reads = int(line.split(' ')[-2].rstrip().replace(',', ''))
                        trimmed_percent = line.split(' ')[-1].rstrip().replace('(', '').replace(')', '')
                        self.sampleinfo[sample]['trimmed_reads'] = trimmed_reads
                        self.sampleinfo[sample]['trimmed_percent'] = trimmed_percent
                        # trim_data.append([sample, trimmed_reads, trimmed_percent])
                        break
                    else:
                        pass
            else:
                print(f'No trim log file found for {sample}. Skipping...')
                continue

        for sample in self.sampleinfo.keys():
            if os.path.isfile(f'{logpath}/{sample}.fastp.json'):
                ## load the json
                json_data = json.load(open(f'{logpath}/{sample}.fastp.json'))
                self.sampleinfo[sample]['merged_reads'] = int(json_data['summary']['after_filtering']['total_reads'])
            else:
                print(f'No dimer log file found for {sample}. Skipping...')
                continue

        return


    def fastqc(self, sample, outprefix):
        '''
        Runs FastQC on fastq files.
        '''

        if not os.path.exists(f'{self.args.output}/trimmed/logs/fastqc'):
            os.mkdir(f'{self.args.output}/trimmed/logs/fastqc')

        ## check if fastqc output already exists
        if os.path.exists(f'{self.args.output}/trimmed/logs/fastqc/{sample}_{outprefix}_fastqc.csv'):
            print(f'FastQC output already exists for {sample}. Skipping...')
            self.sampleinfo[sample][f'{outprefix}_fastqc'] = f'{self.args.output}/trimmed/logs/fastqc/{sample}_{outprefix}_fastqc.csv'
            return
        
        if outprefix == 'raw':
            fastq = self.sampleinfo[sample]['fastq']
        elif outprefix == 'trimmed':
            fastq = self.sampleinfo[sample]['trimmed_fastq']
            fastq = f'{self.args.output}/{fastq}'

        number_of_reads_to_sample = 1000000
        read_lengthD = defaultdict(int)
        read_topSeqD = defaultdict(int)
        reads_sampled = 0
        total_bases = 0
        total_gc_count = 0
        quality_scores = []
        quality_scores_by_position = defaultdict(list)

        print(f'Running FastQC on {fastq}...')

        if fastq.endswith('.gz'):
            file_to_open = gzip.open(fastq, 'rt')
        else:
            file_to_open = open(fastq, 'rt')

        read_length_max = 0

        with file_to_open as handle:
            with tqdm.tqdm(total=number_of_reads_to_sample, desc=sample, unit=' reads') as pbar:
                for i, record in enumerate(SeqIO.parse(handle, 'fastq')):
                    if i % 5 == 0 and i != number_of_reads_to_sample:  # Sample every 5th read
                        pbar.update(5)
                        seq = record.seq
                        total_bases += len(seq)
                        reads_sampled += 1
                        total_gc_count += seq.count('G') + seq.count('C')
                        quality_scores.extend(record.letter_annotations['phred_quality'])
                        read_topSeqD[str(seq)] += 1
                        for pos, score in enumerate(record.letter_annotations['phred_quality']):
                            quality_scores_by_position[pos + 1].append(score)

                        if len(seq) > read_length_max:
                            for q in range(read_length_max + 1, len(seq) + 1):
                                read_lengthD[q] = 0
                                read_length_max = len(seq)

                        read_lengthD[len(seq)] += 1

                    if i == number_of_reads_to_sample:
                        break

        # Calculating statistics
        avg_quality = np.mean(quality_scores) if quality_scores else 0
        gc_content = (total_gc_count / total_bases * 100) if total_bases else 0
        base_quality_stats = self.compute_quality_stats(quality_scores_by_position)

        sample_quality_df = pd.DataFrame(list(base_quality_stats.items()), columns=['position', 'stats'])
        sample_quality_df['min'] = sample_quality_df['stats'].apply(lambda x: x['min'])
        sample_quality_df['lower_quartile'] = sample_quality_df['stats'].apply(lambda x: x['lower_quartile'])
        sample_quality_df['mean'] = sample_quality_df['stats'].apply(lambda x: x['mean'])
        sample_quality_df['upper_quartile'] = sample_quality_df['stats'].apply(lambda x: x['upper_quartile'])
        sample_quality_df['max'] = sample_quality_df['stats'].apply(lambda x: x['max'])
        sample_quality_df = sample_quality_df.drop(columns=['stats'])

        sample_lengths_df = pd.DataFrame(list(read_lengthD.items()), columns=['length', 'count']).sort_values('length', kind='stable')
        sample_lengths_df['percent'] = (sample_lengths_df['count'] / sum(sample_lengths_df['count'])) * 100
        sample_lengths_df['cumulative_percent'] = sample_lengths_df['percent'].cumsum()

        sample_sequences_df = pd.DataFrame(list(read_topSeqD.items()), columns=['sequence', 'count']).sort_values('count', ascending=False, kind='stable')
        sample_sequences_df['percent'] = (sample_sequences_df['count'] / sum(sample_sequences_df['count'])) * 100
        sample_sequences_df['cumulative_percent'] = sample_sequences_df['percent'].cumsum()
        sample_sequences_df = sample_sequences_df.head(100)

        with open(f'{self.args.output}/trimmed/logs/fastqc/{sample}_{outprefix}_fastqc.csv', 'w') as f:
            f.write(f'Sample: {sample}\n')
            f.write(f'Reads sampled: {reads_sampled}\n')
            f.write(f'GC content: {gc_content}\n')
            f.write(f'Average quality: {avg_quality}\n')
            f.write('\n>>Quality\n')
            sample_quality_df.to_csv(f, sep=',', index=False)
            f.write('>>END\n')
            f.write('\n>>Lengths\n')
            sample_lengths_df.to_csv(f, sep=',', index=False)
            f.write('>>END\n')
            f.write('\n>>Overrep_Seqs\n')
            sample_sequences_df.to_csv(f, sep=',', index=False)
            f.write('>>END\n')

        file_to_open.close()
        f.close()



        print(f'FastQC output saved to {self.args.output}/trimmed/logs/fastqc/{sample}_{outprefix}_fastqc.csv')
        print('-' * 50)

        self.sampleinfo[sample][f'{outprefix}_fastqc'] = f'{self.args.output}/trimmed/logs/fastqc/{sample}_{outprefix}_fastqc.csv'

        return

    def compute_quality_stats(self, quality_scores_by_position):
        stats = defaultdict(dict)
        for pos, scores in quality_scores_by_position.items():
            scores_array = np.array(scores)
            stats[pos] = {
                'min': np.min(scores_array),
                'lower_quartile': np.percentile(scores_array, 25),
                'mean': np.mean(scores_array),
                'upper_quartile': np.percentile(scores_array, 75),
                'max': np.max(scores_array)
            }
        return stats

    def subsample_fastq(self):
        '''
        Subsamples fastq files to a specified number of reads.
        '''
        subsample_dir = f'{self.args.output}/subsample'
        if not os.path.exists(subsample_dir):
            os.mkdir(subsample_dir)

        subsample_amount = self.args.subsample
        self.logdata['subsample'] = f'Subsampled reads to {subsample_amount} reads.'
        self.logdata['tools'].append('seqtk')
        self.logdata['commands'].append(f'seqtk sample -s100 -{subsample_amount}')
        print(f'Subsampling reads to {subsample_amount} reads.')

        for sample in self.samplelist:
            output_file = f'{subsample_dir}/{sample}.subsample.fq.gz'

            # Skip if already subsampled
            if os.path.exists(output_file):
                print(f'{output_file} already exists. Skipping...')
                self.fastqDict[sample] = output_file
                continue

            fastq = self.fastqDict[sample]
            cmd = f'seqtk sample -s100 {fastq} {subsample_amount} | pigz > {output_file}'
            run(cmd)

            print(f'Finished subsampling {sample}')
            self.fastqDict[sample] = output_file

    def generate_module_logfile(self):
        '''
        Generates a log file for the module.
        '''

        logdir = f'{self.args.output}/logs'
        
        if not os.path.exists(logdir):
            os.mkdir(logdir)

        log_file = f'{logdir}/trim.log'
        with open(log_file, 'w') as f:
            f.write(f'###\n')
            for key, value in self.logdata.items():
                f.write(f'{key}: {value}\n')
            f.write('###\n')

    def generate_new_samplefile(self, sampleinfo, trimdir):
        '''
        Adds a new column to the sample metadata file with the new fastq file paths.
        '''
        
        if self.args.paired == True:
            data = pd.read_csv(sampleinfo, header=None, names=['sample', 'group', 'fastq', 'fastq_2'])
            data = data.drop(columns=['fastq_2'])
            data['fastq'] = data['sample'].apply(lambda x: f'{trimdir}/{x}.merged.fq.gz')
            data.to_csv(f'{self.args.output}/samples.trimmed.csv', index=False, header=False)

        elif self.args.paired == False:
            data = pd.read_csv(sampleinfo, header=None, names=['sample', 'group', 'fastq'])
            data['fastq'] = data['sample'].apply(lambda x: f'{trimdir}/{x}.trimmed.fq.gz')
            data.to_csv(f'{self.args.output}/samples.trimmed.csv', index=False, header=False)

        return

    def generate_new_samplefile_subsampleonly(self):
        '''
        Adds a new column to the sample metadata file with the new fastq file paths.
        '''

        if self.args.paired == True:
            self.sampleinfo['fastq'] = f'{self.args.output}/subsample/' + self.sampleinfo['sample'] + '.merged.subsample.fq.gz'
            self.sampleinfo = self.sampleinfo.drop(columns=['fastq_2'])
        elif self.args.paired == False:
            self.sampleinfo['fastq'] = f'{self.args.output}/subsample/' + self.sampleinfo['sample'] + '.subsample.fq.gz'
            
        self.sampleinfo.to_csv(f'{self.args.output}/samples.collapse.csv', index=False, header=False)

def main():
    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('-s', '--samples', required=False, help='Sample metadata file. Format: sample,group,fastq', default='default')
    ap.add_argument('-@', '--threads', required=False, help='Number of threads to process parallel fastq files with. Default=4', default=4)
    ap.add_argument('-o', '--output', required=True, help='Output prefix.')
    ap.add_argument('-x', '--subsample', required=False, help='Subsample reads to this number. Default=0; Use all reads.', default=0)
    ap.add_argument('--paired', required=False, help='Paired-end reads', action='store_true')
    ap.add_argument('--analyze', required=False, help='Whether or not to perform QC analyses.', action='store_true')
    ap.add_argument('--library', required=False, help='Library type. RSB (RealSeq); SMART (SMART-seq)', default='RSB')
    ap.add_argument('--tool', required=False, help='Tool to use for trimming.', default='cutadapt')
    ap.add_argument('--fastqdir', required=False, help='Directory containing fastq files. Appends this to the beginning of each file.', default='default')
    args = ap.parse_args()

    startTime = datetime.now()

    tr = TrimReads(args)
    print('-' * 50)
    print('-' * 50)
    print('Collapsing Reads...')
    print('-' * 50)
    print('-' * 50)

    ## check if the samplefile has 3 or 4 columns so we can determine if the reads are paired-end or not
    #if args.samples != 'default':
    #    with open(args.samples, 'r') as f:
    #        line = f.readline()
    #        if len(line.split(',')) == 4:
    #            args.paired = True
    #        else:
    #            args.paired = False
    #else:
    #    print('No sample metadata file provided. Please provide a sample metadata file.')
    #    sys.exit()
        
    if args.samples == 'default':
        print('No sample metadata file provided. Please provide a sample metadata file.')
        sys.exit(1)

    if not args.paired:          # only auto-detect when flag was not explicitly passed
        with open(args.samples) as f:
            args.paired = len(f.readline().strip().split(',')) == 4

    # Ensure output directories exist
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    trimdir = f'{args.output}/trimmed'

    ## Generate sample metadata
    tr.generate_sampleinfo()

    ## Feed the samples to the appropriate trimming function based on the library type
    tr.feed_trimming()

    ## Parse log files
    tr.parse_logfiles_rsb()

    tr.generate_module_logfile()

    tr.generate_new_samplefile(args.samples, trimdir)

    print('-' * 50)
    print('-' * 50)
    print(f'Time elapsed for read trimming: {datetime.now() - startTime}')
    print('-' * 50)
    print('-' * 50)

if __name__ == "__main__":

    main()

