#!/usr/bin/env python3
#Aidan Manning 1.28.25
#Script to map small RNA reads to the genome, ncRNAs, and other RNAs using various mapping tools.
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import warnings
from collections import defaultdict
from dbPaths import read_dbpaths
from toolRun import run  # stops on a failed tool
warnings.filterwarnings("ignore")

class MapReads(object):

    def __init__(self, args):
        '''
        This class is used to run FastQC on a set of fastq files. It will run FastQC on each fastq file in parallel using the number of threads specified by the user.
        '''

        self.args = args

        self.sampleinfo = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq'])
        self.samplelist = self.sampleinfo['sample'].tolist()
        self.sampleDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['group']))
        self.fastqDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['fastq']))

        # relative index paths are made absolute: against --dbs-base (the references folder, for a replaced
        # tier list) or the folder of the dbpaths file (dbPaths.py)
        self.dbs = pd.DataFrame(read_dbpaths(self.args.dbs, getattr(self.args, 'dbs_base', None)), columns=['map_cmd', 'db_path', 'retain'])

        # Using a regular dictionary for logdata (avoiding lambda functions)
        self.logdata = {}
        self.logdata['module'] = 'mapReads'
        self.logdata['date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.logdata['sample_input'] = self.fastqDict
        self.logdata['databases'] = self.dbs
        self.logdata['threads'] = self.args.threads
        self.logdata['output'] = f'{self.args.output}/alignments'
        self.logdata['tools'] = []
        self.logdata['commands'] = []

        ## get bowtie and bowtie2 versions
        bowtie_version = os.popen('bowtie --version').read().split('\n')[0].split(' ')[2]
        self.logdata['tools'].append(f'bowtie={bowtie_version}')
        bowtie2_version = os.popen('bowtie2 --version').read().split('\n')[0].split(' ')[2]
        self.logdata['tools'].append(f'bowtie2={bowtie2_version}')
        samtools_version = os.popen('samtools --version').read().split('\n')[0].split(' ')[1]
        self.logdata['tools'].append(f'samtools={samtools_version}')

        return
    
    #def sam_to_bam(self, sam_file):
    #    '''
    #    Converts a sam file to a bam file.
    #    '''
    #    print('-------------------SAM TO BAM----------------------')
    #    threads = self.args.threads#

    #    bam_file = sam_file.replace('.sam', '.bam')
    #    cmd = f'samtools view --threads {threads} -bS {sam_file} > {bam_file}'
    #    os.system(cmd)

    #    ## remove the sam file
    #    os.remove(sam_file)

    #    return bam_file
    
    def sam_to_bam(self, sam_file, filter_unaligned=True):
        '''
        Converts a sam file to a bam file.
    
        filter_unaligned: when True, strips unmapped reads (-F 4) from the BAM.
        This is necessary because bowtie1 does not support --no-unal, so 
        unaligned reads appear in its SAM output. Without filtering, these
        duplicate unmapped records inflate the Unmapped/Total counts in 
        mappingstats after BAM merging. The last tier should pass 
        filter_unaligned=False so that truly unmapped reads are preserved 
        for parseMappings to count.
        '''
        print('-------------------SAM TO BAM----------------------')
        threads = self.args.threads

        bam_file = sam_file.replace('.sam', '.bam')
        flag_filter = '-F 4 ' if filter_unaligned else ''
        cmd = f'samtools view --threads {threads} {flag_filter}-bS {sam_file} > {bam_file}'
        run(cmd)

        ## remove the sam file
        os.remove(sam_file)

        return bam_file

    def merge_bams(self, sample, bams):
        '''
        Merges bam files.
        '''
        threads = self.args.threads

        if len(bams) == 1:
            ## Move the file to the correct location
            cmd = f'mv {bams[0]} {self.args.output}/alignments/{sample}.tmp.bam'
            run(cmd)
            ## Sort the bam file by name
            cmd = f'samtools sort -n -@ {threads} {self.args.output}/alignments/{sample}.tmp.bam -o {self.args.output}/alignments/{sample}.merged.bam'
            run(cmd)

            ## remove the tmp bam file
            os.remove(f'{self.args.output}/alignments/{sample}.tmp.bam')

        else:
            cmd = f'samtools merge -@ {threads} {self.args.output}/alignments/{sample}.tmp.bam {" ".join(bams)}'
            run(cmd)

            ## Sort the bam file by name
            cmd = f'samtools sort -n -@ {threads} {self.args.output}/alignments/{sample}.tmp.bam -o {self.args.output}/alignments/{sample}.merged.bam'
            run(cmd)

            ## remove the tmp bam file
            os.remove(f'{self.args.output}/alignments/{sample}.tmp.bam')

        return f'{self.args.output}/alignments/{sample}.merged.bam'

    def verify_input_reads(self):
        '''
        Verifies that the input reads are in the correct format.
        '''
        num_samples = len(self.samplelist)

        for sample in self.samplelist:
            fastq = self.fastqDict[sample]
            if not os.path.exists(fastq):
                print(f'{fastq} does not exist. Exiting...')
                sys.exit(1)

        print(f'{num_samples}/{num_samples} input read files verified.')

        return
    
    def verify_bowtie_database(self):
        '''
        Verifies that the bowtie database exists.
        '''

        num_dbs = len(self.dbs)

        for db in self.dbs['db_path']:
            if not os.path.exists(f'{db}.1.bt2'):
                ## check for a normal bowtie database
                if not os.path.exists(f'{db}.1.ebwt'):
                    ## check for bt2l bowtie database
                    if not os.path.exists(f'{db}.1.bt2l'):
                        print(f'{db} does not exist. Exiting...')
                        sys.exit(1)
        print(f'{num_dbs}/{num_dbs} bowtie databases verified.')

        return

    def map_reads(self, sample):
        '''
        Main function to map reads to the bowtie databases.
        '''
        multimap = self.args.multimap
        threads = self.args.threads


                  
        print("FASTQ inputs:")
        for s, fq in self.fastqDict.items():
            print(s, fq, os.path.exists(fq))

        # Define mapping commands
        #mapping_cmds = {
        #    'mind': 'bowtie --fullref --no-unal -k 1 -v 0 --best --norc -x {genome} -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq -S {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'bowtie': 'bowtie --no-unal -k 1 -n 1 -q -l 18 -e 80 --chunkmbs 256 -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq -x {genome} -S {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'bowtie1': 'bowtie --fullref --no-unal --best --strata -n 1 -q -l 18 -e 80 --chunkmbs 256 -m {multimap} -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq -x {genome} -S {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'ribospike': 'bowtie --fullref --no-unal -k 1 -n 2 --best --norc -x {genome} -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq -S {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'ribospikebt2': 'bowtie2 --sam-no-qname-trunc --no-unal --no-1mm-upfront --mp 6 --rdg 6,3 --rfg 6,3 --norc --n-ceil L,5,0.15 -k 1 -D 20 -R 2 -N 0 -L 10 -i S,1,0.40 --score-min L,0,-0.6 --np 1 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'bowtie2': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -k {multimap} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'bowtie2best': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #    'bowtie2sense': 'bowtie2 --sam-no-qname-trunc --norc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
        #}
        
        
        mapping_cmds = {

            'bowtie': (
                'zcat -f {fastq_file} | '
                'bowtie -S -k 1 -n 1 -l 18 -e 80 --chunkmbs 256 '
                '-p {threads} '
                '--un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq '
                '{genome} - {outprefix}/alignments/{sample}.{dbprefix}.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log'
            ),
            

            'bowtie1': (
                'bowtie --fullref --no-unal --best --strata -n 1 -l 18 -e 80 --chunkmbs 256 '
                '-m {multimap} -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log'
            ),
            'mind': (
                'bowtie --fullref --no-unal -k 1 -v 0 --best --norc '
                '-p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log'
            ),
            'ribospike': (
                'bowtie --fullref --no-unal -k 1 -n 2 --best --norc '
                '-p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log'
            ),        
            'ribospikebt2': 'bowtie2 --sam-no-qname-trunc --no-unal --no-1mm-upfront --mp 6 --rdg 6,3 --rfg 6,3 --norc --n-ceil L,5,0.15 -k 1 -D 20 -R 2 -N 0 -L 10 -i S,1,0.40 --score-min L,0,-0.6 --np 1 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
            'bowtie2': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -k {multimap} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
            'bowtie2best': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
            'bowtie2sense': 'bowtie2 --sam-no-qname-trunc --norc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.log',
            }
        
        
       
        
        dblist = self.dbs['db_path'].tolist()
        dblist = [db.split('/')[-1] for db in dblist]

        # Log commands if not already logged
        for idx, row in self.dbs.iterrows():
            cmd = row['map_cmd']
            if f'{cmd} : {mapping_cmds[cmd]}' not in self.logdata['commands']:
                self.logdata['commands'].append(f'{cmd} : {mapping_cmds[cmd]}')

        # Check if the final BAM file already exists
        final_bam = f'{self.args.output}/alignments/{sample}.merged.bam'
        if os.path.exists(final_bam):
            print(f'{sample} already has all alignments. Skipping...')
            return

        sample_start = defaultdict(lambda: 0)
        sample_bams = defaultdict(list)
        sample_bam_rm = defaultdict(list)
        sample_fq_rm = defaultdict(list)

        fastq = self.fastqDict[sample]

        # Check if the input file is a fasta or fastq file and set the appropriate flag
        if fastq.endswith('fa.gz') or fastq.endswith('fasta.gz') or fastq.endswith('fna.gz') or fastq.endswith('fa') or fastq.endswith('fasta') or fastq.endswith('fna'):
            toggle_fasta = True
        elif fastq.endswith('fq.gz') or fastq.endswith('fastq.gz') or fastq.endswith('fq') or fastq.endswith('fastq'):
            toggle_fasta = False
        else:
            print(f'{fastq} is not a valid fastq file. Exiting...')
            sys.exit(1)

        # Check if any BAM files already exist
        for idx, row in self.dbs.iterrows():

            db = row['db_path']
            dbprefix = db.split('/')[-1]
            retain = row['retain']
            
            db = row['db_path']
            retain = row['retain']
            map_cmd = row['map_cmd']

            # Check for valid mapping command
            if map_cmd not in mapping_cmds:
                print(f'{map_cmd} is not a valid mapping command. Exiting...')
                sys.exit(1)

            print(f'Mapping reads for {sample} to {db}')
            dbprefix = dblist[idx].split('/')[-1]

            if not os.path.exists(fastq):
                print(f'{fastq} does not exist. Exiting...')
                sys.exit(1)

            # Prepare and run the mapping command
            cmd = mapping_cmds[map_cmd].format(threads=threads, genome=db, multimap=multimap, fastq_file=fastq, outprefix=self.args.output, sample=sample, dbprefix=dbprefix)
            
            if idx == len(dblist) - 1:
                cmd = cmd.replace('--no-unal', '')  # Keep unaligned reads for the last mapping

            if toggle_fasta:
                cmd = cmd.replace('-S', '-f -S') # Add -f flag for fasta files
           
            
            run(cmd)





            if map_cmd in ['bowtie2', 'bowtie2best', 'bowtie2sense', 'ribospikebt2']:
                fastq = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.fq.gz'
            elif map_cmd in ['bowtie', 'ribospike', 'bowtie1', 'mind']:
                fastq = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.fq'

            # Convert SAM to BAM
            sam_file = f'{self.args.output}/alignments/{sample}.{dbprefix}.sam'
            #bam_file = self.sam_to_bam(sam_file)
            is_last_tier = (idx == len(dblist) - 1)
            bam_file = self.sam_to_bam(sam_file, filter_unaligned=not is_last_tier)

            # Track BAM files for merging
            if retain == 'T':
                sample_bams[sample].append(bam_file)
            sample_bam_rm[sample].append(bam_file)

        # Merge the BAM files and clean up
        self.merge_bams(sample, sample_bams[sample])

        ## delete all unaligned bam files
        for bam in sample_bam_rm[sample]:
            try:
                os.remove(bam)
            except:
                pass

        ## delete all unaligned fastq files
        run(f'rm -f {self.args.output}/alignments/unal/{sample}*.fq*')

        return

    @staticmethod
    def _is_fastq_empty(fq_path):
        if not os.path.exists(fq_path) or os.path.getsize(fq_path) == 0:
            return True
        import gzip
        try:
            open_fn = gzip.open if fq_path.endswith('.gz') else open
            with open_fn(fq_path, 'rt') as f:
                line = f.readline()
                return len(line.strip()) == 0
        except Exception:
            return True

    @staticmethod
    def _filter_fastq_by_length(in_fq, out_fq, min_len=20):
        import gzip
        open_in = gzip.open if in_fq.endswith('.gz') else open
        open_out = gzip.open if out_fq.endswith('.gz') else open
        count = 0
        with open_in(in_fq, 'rt') as fin, open_out(out_fq, 'wt') as fout:
            while True:
                h = fin.readline()
                if not h:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                if not qual:
                    break
                if len(seq.strip()) >= min_len:
                    fout.write(h)
                    fout.write(seq)
                    fout.write(plus)
                    fout.write(qual)
                    count += 1
        return count

    def map_reads_iterative(self, sample):
        '''
        Hierarchical Exact-First Fallback Alignment for Degraded cfRNA:
        Pass 1: Strict exact matching across all tiers (m=0).
        Pass 2: Biological variation fallback on residual unmapped reads:
                - Compact tiers (tRNA, miRNA, piRNA, etc.): allow 1 mismatch down to 15 nt.
                - Massive tiers (ncRNA, cDNA, genome): allow 1 mismatch only for reads >= 20 nt.
        '''
        multimap = self.args.multimap
        threads = self.args.threads

        # Check if the final BAM file already exists
        final_bam = f'{self.args.output}/alignments/{sample}.merged.bam'
        if os.path.exists(final_bam):
            print(f'{sample} already has all alignments. Skipping...')
            return

        fastq = self.fastqDict[sample]
        if fastq.endswith('fa.gz') or fastq.endswith('fasta.gz') or fastq.endswith('fa'):
            toggle_fasta = True
        else:
            toggle_fasta = False

        dblist = [db.split('/')[-1] for db in self.dbs['db_path'].tolist()]
        massive_keywords = ['ncrna', 'cdna', 'genome']

        # Pass 1: Strict Exact Matching
        # Bowtie 1: -v 0
        # Bowtie 2: --score-min L,0,0
        mapping_cmds_pass1 = {
            'bowtie': (
                'zcat -f {fastq_file} | '
                'bowtie -S -k 1 -v 0 --chunkmbs 256 '
                '-p {threads} '
                '--un {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq '
                '{genome} - {outprefix}/alignments/{sample}.{dbprefix}.p1.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log'
            ),
            'bowtie1': (
                'bowtie --fullref --no-unal --best --strata -v 0 --chunkmbs 256 '
                '-m {multimap} -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.p1.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log'
            ),
            'mind': (
                'bowtie --fullref --no-unal -k 1 -v 0 --best --norc '
                '-p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.p1.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log'
            ),
            'ribospike': (
                'bowtie --fullref --no-unal -k 1 -v 0 --best --norc '
                '-p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.p1.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log'
            ),
            'ribospikebt2': 'bowtie2 --sam-no-qname-trunc --no-unal --no-1mm-upfront --mp 6 --rdg 6,3 --rfg 6,3 --norc --n-ceil L,5,0.15 -k 1 -D 20 -R 2 -N 0 -L 10 -i S,1,0.40 --score-min L,0,0 --np 1 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p1.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log',
            'bowtie2': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --score-min L,0,0 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq.gz -x {genome} -k {multimap} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p1.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log',
            'bowtie2best': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --score-min L,0,0 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p1.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log',
            'bowtie2sense': 'bowtie2 --sam-no-qname-trunc --norc --no-unal --np 5 --very-sensitive --score-min L,0,0 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p1.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p1.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p1.log',
        }

        # Pass 2: Fallback with 1 Mismatch
        mapping_cmds_pass2 = {
            'bowtie': (
                'zcat -f {fastq_file} | '
                'bowtie -S -k 1 -n 1 -l 18 -e 80 --chunkmbs 256 '
                '-p {threads} '
                '--un {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq '
                '{genome} - {outprefix}/alignments/{sample}.{dbprefix}.p2.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log'
            ),
            'bowtie1': (
                'bowtie --fullref --no-unal --best --strata -n 1 -l 18 -e 80 --chunkmbs 256 '
                '-m {multimap} -p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.p2.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log'
            ),
            'mind': (
                'bowtie --fullref --no-unal -k 1 -v 1 --best --norc '
                '-p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.p2.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log'
            ),
            'ribospike': (
                'bowtie --fullref --no-unal -k 1 -n 2 --best --norc '
                '-p {threads} --un {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq '
                '-S {genome} {fastq_file} {outprefix}/alignments/{sample}.{dbprefix}.p2.sam '
                '2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log'
            ),
            'ribospikebt2': 'bowtie2 --sam-no-qname-trunc --no-unal --no-1mm-upfront --mp 6 --rdg 6,3 --rfg 6,3 --norc --n-ceil L,5,0.15 -k 1 -D 20 -R 2 -N 0 -L 10 -i S,1,0.40 --score-min L,0,-0.6 --np 1 --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p2.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log',
            'bowtie2': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq.gz -x {genome} -k {multimap} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p2.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log',
            'bowtie2best': 'bowtie2 --sam-no-qname-trunc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p2.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log',
            'bowtie2sense': 'bowtie2 --sam-no-qname-trunc --norc --no-unal --np 5 --very-sensitive --ignore-quals --threads {threads} --un-gz {outprefix}/alignments/unal/{sample}.{dbprefix}.p2.fq.gz -x {genome} -U {fastq_file} -S {outprefix}/alignments/{sample}.{dbprefix}.p2.sam 2> {outprefix}/alignments/logs/{sample}.{dbprefix}.p2.log',
        }

        sample_bams_to_merge = []
        files_to_clean = []

        # ==========================================
        # PASS 1: Exact Match Cascade (All Tiers)
        # ==========================================
        print(f"[Iterative Alignment] Pass 1 (Exact Match, m=0) for {sample}")
        curr_fq = fastq

        for idx, row in self.dbs.iterrows():
            db = row['db_path']
            dbprefix = dblist[idx]
            retain = row['retain']
            map_cmd = row['map_cmd']

            if not os.path.exists(curr_fq) or self._is_fastq_empty(curr_fq):
                break

            cmd = mapping_cmds_pass1[map_cmd].format(
                threads=threads, genome=db, multimap=multimap,
                fastq_file=curr_fq, outprefix=self.args.output,
                sample=sample, dbprefix=dbprefix
            )
            if toggle_fasta:
                cmd = cmd.replace('-S', '-f -S')

            run(cmd)

            if map_cmd in ['bowtie2', 'bowtie2best', 'bowtie2sense', 'ribospikebt2']:
                next_fq = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.p1.fq.gz'
            else:
                next_fq = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.p1.fq'

            curr_fq = next_fq
            files_to_clean.append(curr_fq)

            sam_file = f'{self.args.output}/alignments/{sample}.{dbprefix}.p1.sam'
            if os.path.exists(sam_file):
                bam_file = self.sam_to_bam(sam_file, filter_unaligned=True)
                if retain == 'T':
                    sample_bams_to_merge.append(bam_file)
                files_to_clean.append(bam_file)

        unaligned_pass1_fq = curr_fq

        # ==========================================
        # PASS 2: Biological Variation Fallback
        # ==========================================
        print(f"[Iterative Alignment] Pass 2 (Biological Fallback) for {sample}")
        curr_fq_p2 = unaligned_pass1_fq

        for idx, row in self.dbs.iterrows():
            db = row['db_path']
            dbprefix = dblist[idx]
            retain = row['retain']
            map_cmd = row['map_cmd']
            is_last_tier = (idx == len(dblist) - 1)

            if not os.path.exists(curr_fq_p2) or self._is_fastq_empty(curr_fq_p2):
                break

            is_massive = any(kw in dbprefix.lower() for kw in massive_keywords)

            if is_massive:
                filtered_fq = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.ge20.p2.fq'
                kept_count = self._filter_fastq_by_length(curr_fq_p2, filtered_fq, min_len=20)
                files_to_clean.append(filtered_fq)
                if kept_count == 0:
                    continue
                in_fq_for_tier = filtered_fq
            else:
                in_fq_for_tier = curr_fq_p2

            cmd = mapping_cmds_pass2[map_cmd].format(
                threads=threads, genome=db, multimap=multimap,
                fastq_file=in_fq_for_tier, outprefix=self.args.output,
                sample=sample, dbprefix=dbprefix
            )
            if is_last_tier:
                cmd = cmd.replace('--no-unal', '')
            if toggle_fasta:
                cmd = cmd.replace('-S', '-f -S')

            run(cmd)

            if map_cmd in ['bowtie2', 'bowtie2best', 'bowtie2sense', 'ribospikebt2']:
                next_fq_p2 = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.p2.fq.gz'
            else:
                next_fq_p2 = f'{self.args.output}/alignments/unal/{sample}.{dbprefix}.p2.fq'

            curr_fq_p2 = next_fq_p2
            files_to_clean.append(curr_fq_p2)

            sam_file = f'{self.args.output}/alignments/{sample}.{dbprefix}.p2.sam'
            if os.path.exists(sam_file):
                bam_file = self.sam_to_bam(sam_file, filter_unaligned=not is_last_tier)
                if retain == 'T' or is_last_tier:
                    sample_bams_to_merge.append(bam_file)
                files_to_clean.append(bam_file)

        # Merge all BAMs
        self.merge_bams(sample, sample_bams_to_merge)

        # Clean intermediate files
        for f in files_to_clean:
            try:
                os.remove(f)
            except:
                pass
        run(f'rm -f {self.args.output}/alignments/unal/{sample}*.fq*')
        return

    def generate_module_logfile(self):
        '''
        Generates a log file for the module.
        '''
        log_file = f'{self.args.output}/logs/map.log'
        with open(log_file, 'w') as f:
            f.write(f'###\n')
            for key, value in self.logdata.items():
                if key == 'sample_input':
                    f.write(f'{key}:\n')
                    for sample, fastq in value.items():
                        f.write(f'{sample}: {fastq}\n')
                elif key == 'commands':
                    f.write(f'{key}:\n')
                    for cmd in value:
                        f.write(f'{cmd}\n')
                elif key == 'databases':
                    f.write(f'{key}:\n')
                    for idx, row in value.iterrows():
                        f.write(f'{row["map_cmd"]}: {row["db_path"]} {row["retain"]}\n')
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
    ap.add_argument('-d', '--dbs', required=True, help='Bowtie databases to map reads to. Format: map_cmd,db_path,T/F retain reads that map to db.')
    ap.add_argument('--dbs-base', required=False, default=None, help='Folder that relative index paths in --dbs are resolved against. Default: the folder of the --dbs file.')
    ap.add_argument('-m', '--multimap', required=False, help='Number of multimapping loci a read can map to. Default=100', default=100)
    ap.add_argument('--align-mode', dest='align_mode', choices=['standard', 'iterative'], default='standard', help='Alignment mode: standard (single-pass) or iterative (exact-first fallback for degraded cfRNA)')
    args = ap.parse_args()
    
    startTime = datetime.now()

    mr =  MapReads(args)

    print(f'-'*50)
    print(f'-'*50)
    print(f'Mapping reads.')
    print(f'-'*50)
    print(f'-'*50)

    ## Create all the necessary directories
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    logfile_dir = f'{args.output}/logs'
    if not os.path.exists(logfile_dir):
        os.mkdir(logfile_dir)

    if not os.path.exists(f'{args.output}/alignments'):
        os.mkdir(f'{args.output}/alignments')

    if not os.path.exists(f'{args.output}/alignments/logs'):
        os.mkdir(f'{args.output}/alignments/logs')

    if not os.path.exists(f'{args.output}/alignments/unal'):
        os.mkdir(f'{args.output}/alignments/unal')

    ## Verify that the input reads exist
    mr.verify_input_reads()

    ## Verify that the bowtie databases exist
    mr.verify_bowtie_database()

    print(f'-'*50)
    print(f'-'*50)

    #Map reads
    samplelist = mr.samplelist

    for idx, sample in enumerate(samplelist):

        print(f'Mapping reads for {sample} ({idx+1}/{len(samplelist)})')

                # Check if sample's BAM file exists and skip if already processed
        if os.path.exists(f'{args.output}/alignments/{sample}.merged.bam'):
            print(f'{sample} already has all alignments. Skipping...')
            continue

        if getattr(args, 'align_mode', 'standard') == 'iterative':
            mr.map_reads_iterative(sample)
        else:
            mr.map_reads(sample)
        print(f'-'*50)

    mr.generate_module_logfile()
    
    print(f'-'*50)
    print(f'-'*50)
    print('Time elasped for read mapping: ', datetime.now() - startTime)
    print(f'-'*50)
    print(f'-'*50)

if __name__ == "__main__":
    main()