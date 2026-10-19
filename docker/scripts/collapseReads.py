#!/usr/bin/env python3
# Aidan Manning 7.9.24
# Script to map small RNA reads to the genome, ncRNAs, and other RNAs using various mapping tools

import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import warnings
from multiprocessing import Pool
import subprocess
from toolRun import run, ToolError  # stops on a failed tool

warnings.filterwarnings("ignore")

class CollapseReads:
    
    def __init__(self, args):
        '''
        This class is used to handle small RNA read processing and collapsing using vsearch.
        '''
        self.args = args

        self.sampleinfo = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq'])
        self.samplelist = self.sampleinfo['sample'].tolist()
        self.sampleDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['group']))
        self.fastqDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['fastq']))

        # Using a regular dictionary for logdata (avoiding lambda functions)
        self.logdata = {}
        self.logdata['module'] = 'collapseReads'
        self.logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.logdata['samples'] = self.samplelist
        self.logdata['threads'] = self.args.threads
        self.logdata['output'] = self.args.output
        self.logdata['tools'] = []
        self.logdata['commands'] = []

    def collapse_reads(self, sample):
        '''
        Collapses reads in a fastq file to unique sequences and counts the number of times each sequence is observed.
        '''

        fastq = self.fastqDict[sample]

        if self.args.fastqdir != 'default':

            if self.args.subsample > 0:
                fastq = fastq
            else:
                if self.args.fastqdir[-1] != '/':
                    fastq = f'{self.args.fastqdir}/{fastq}'
                else:
                    fastq = f'{self.args.fastqdir}{fastq}'

        output_file = f'{self.args.output}/collapse/{sample}.collapse.fa'

        # Skip collapsing if output already exists
        if os.path.exists(output_file):
            print(f'{output_file} already exists. Skipping...')
            return

        cmd = [
            'vsearch',
            '--fastx_uniques', fastq,
            '--fastaout', output_file,
            '--sizeout',
            '--quiet'
        ]

        # Running the command using subprocess
        result = subprocess.run(cmd, capture_output=True)
        if result.returncode == 0:
            print(f'Finished collapsing reads for {sample}')
        else:
            raise ToolError(f'vsearch failed for {sample} (exit status {result.returncode}):\n{result.stderr.decode()}')

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
                if self.args.fastqdir != 'default':
                    if self.args.fastqdir[-1] != '/':
                        output_file = output_file.replace(f'{self.args.fastqdir}/', '')
                    else:
                        output_file = output_file.replace(f'{self.args.fastqdir}', '')

                self.fastqDict[sample] = output_file
                continue

            fastq = self.fastqDict[sample]

            if self.args.fastqdir != 'default':
                if self.args.fastqdir[-1] != '/':
                    fastq = f'{self.args.fastqdir}/{fastq}'
                else:
                    fastq = f'{self.args.fastqdir}{fastq}'

            cmd = f'seqtk sample -s100 {fastq} {subsample_amount} | pigz > {output_file}'
            # print(cmd)
            run(cmd)

            if self.args.fastqdir != 'default':
                if self.args.fastqdir[-1] != '/':
                    output_file = output_file.replace(f'{self.args.fastqdir}/', '')
                else:
                    output_file = output_file.replace(f'{self.args.fastqdir}', '')
            
            print(f'Finished subsampling {sample}')
            self.fastqDict[sample] = output_file

    def generate_module_logfile(self):
        '''
        Generates a log file for the module.
        '''
        log_file = f'{self.args.output}/logs/collapse.log'
        with open(log_file, 'w') as f:
            f.write(f'###\n')
            for key, value in self.logdata.items():
                f.write(f'{key}: {value}\n')
            f.write('###\n')

    def generate_new_samplefile(self):
        '''
        Adds a new column to the sample metadata file with the new fastq file paths.
        '''
        self.sampleinfo['fastq'] = f'{self.args.output}/collapse/' + self.sampleinfo['sample'] + '.collapse.fa'
        self.sampleinfo.to_csv(f'{self.args.output}/samples.collapse.csv', index=False, header=False)

    def generate_new_samplefile_subsampleonly(self):
        '''
        Adds a new column to the sample metadata file with the new fastq file paths.
        '''
        self.sampleinfo['fastq'] = f'{self.args.output}/subsample/' + self.sampleinfo['sample'] + '.subsample.fq.gz'
        self.sampleinfo.to_csv(f'{self.args.output}/samples.collapse.csv', index=False, header=False)

def main():

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('-s', '--samples', required=True, help='Sample metadata file. Format: sample,group,fastq')
    ap.add_argument('-@', '--threads', required=False, help='Number of threads to process parallel fastq files with. Default=4', default=4)
    ap.add_argument('-o', '--output', required=True, help='Output prefix.')
    ap.add_argument('-x', '--subsample', required=False, help='Subsample reads to this number. Default=0; Use all reads.', default=0, type=int)
    ap.add_argument('--subsampleonly', required=False, help='Only subsample reads.', action='store_true')
    ap.add_argument('--fastqdir', required=False, help='Directory containing fastq files. Appends this to the beginning of each file.', default='default')
    args = ap.parse_args()

    startTime = datetime.now()

    cr = CollapseReads(args)
    print('-' * 50)
    print('-' * 50)
    print('Collapsing Reads...')
    print('-' * 50)
    print('-' * 50)

    # Ensure output directories exist
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    collapse_dir = f'{args.output}/collapse'
    if not os.path.exists(collapse_dir):
        os.mkdir(collapse_dir)

    logfile_dir = f'{args.output}/logs'
    if not os.path.exists(logfile_dir):
        os.mkdir(logfile_dir)

    # Subsample reads if requested
    if int(args.subsample) > 0:
        cr.subsample_fastq()
        print('-' * 50)

        if args.subsampleonly:
            cr.generate_new_samplefile_subsampleonly()
            cr.logdata['tools'].append('seqtk')
            cr.logdata['commands'].append(f'seqtk sample -s100 {args.subsample}')
            cr.generate_module_logfile()
            print(f'Time elapsed for read subsampling: {datetime.now() - startTime}')
            print('-' * 50)
            print('-' * 50)
            sys.exit()

    # Collapse reads in parallel
    with Pool(int(args.threads)) as p:
        p.map(cr.collapse_reads, cr.samplelist)

    cr.logdata['tools'].append('vsearch')
    cr.logdata['commands'].append('vsearch --fastx_uniques --sizeout')

    # Generate log file
    cr.generate_module_logfile()

    # Generate new sample metadata file
    cr.generate_new_samplefile()

    print('-' * 50)
    print('-' * 50)
    print(f'Time elapsed for read collapse: {datetime.now() - startTime}')
    print('-' * 50)
    print('-' * 50)


if __name__ == "__main__":
    main()
