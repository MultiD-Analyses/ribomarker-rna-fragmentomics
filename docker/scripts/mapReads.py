#!/usr/bin/env python3
#Aidan Manning 1.28.25
#Script to map small RNA reads to the genome, ncRNAs, and other RNAs using various mapping tools.
import argparse
import gzip
import shutil
import sys
import pandas as pd
import os
from datetime import datetime
from dbPaths import read_dbpaths
import mapSchedule
from toolRun import run  # stops on a failed tool


# ----------------------------------------------------------------------------- FASTQ helpers (release 0.4.0)
def open_fastq(path, mode='rt'):
    """A FASTQ file, gzipped if its name ends in .gz."""
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def fastq_records(path):
    """(header line without '@', sequence, quality) of each record."""
    with open_fastq(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            sequence, _, quality = f.readline().rstrip('\n'), f.readline(), f.readline().rstrip('\n')
            yield header.rstrip('\n')[1:], sequence, quality


def has_reads(path):
    with open_fastq(path) as f:
        return bool(f.readline())


def split_by_length(path, inside, outside, min_length, max_length):
    """Writes the reads of path whose length is within [min_length, max_length] ('-' = no bound) to inside and the
    others to outside, unchanged and in order; returns (reads inside, reads outside)."""
    low = int(min_length) if min_length != '-' else 0
    high = int(max_length) if max_length != '-' else float('inf')
    counts = [0, 0]
    with open_fastq(inside, 'wt') as a, open_fastq(outside, 'wt') as b:
        for name, sequence, quality in fastq_records(path):
            keep = low <= len(sequence) <= high
            (a if keep else b).write(f'@{name}\n{sequence}\n+\n{quality}\n')
            counts[0 if keep else 1] += 1
    return counts[0], counts[1]


def append_fastq(source, target):
    """Appends the reads of source to target (each gzipped or not by its name)."""
    with open_fastq(source) as src, open_fastq(target, 'at') as dst:
        shutil.copyfileobj(src, dst)


def unmapped_record(header, name, sequence, quality, aligner):
    """An unaligned SAM record in the form the aligner writes it: Bowtie 2 (run with
    --sam-no-qname-trunc) keeps the whole name line and adds YT:Z:UU; Bowtie 1 cuts the name at the first whitespace and
    adds XM:i:0 [background knowledge; compared with the aligners' own records in docker/tests/tools]."""
    import pysam
    record = pysam.AlignedSegment(header)
    record.query_name = name if aligner == 'bowtie2' else name.split()[0]
    record.flag = 4
    record.query_sequence = sequence
    record.query_qualities = pysam.qualitystring_to_array(quality)
    record.set_tag('YT' if aligner == 'bowtie2' else 'XM', 'UU' if aligner == 'bowtie2' else 0)
    return record

class MapReads(object):

    def __init__(self, args):
        '''
        Maps each sample's reads to the references (tiers) of the tier list, one after another.
        '''

        self.args = args

        self.sampleinfo = pd.read_csv(self.args.samples, sep=',', names=['sample', 'group', 'fastq'])
        self.samplelist = self.sampleinfo['sample'].tolist()
        self.fastqDict = dict(zip(self.sampleinfo['sample'], self.sampleinfo['fastq']))

        # the cells to run, in order (tables/schedule.tsv, mapSchedule.py); the index paths come from the references
        # folder's dbpaths file, relative paths made absolute against the folder of that file (dbPaths.py)
        self.cells = mapSchedule.read_schedule(self.args.schedule)
        self.index_of = {os.path.basename(index): index for _, index, _ in read_dbpaths(self.args.dbs)}
        # each reference once, in the order in effect (tables/order.tsv beside the schedule; also when pass 1 does not
        # cover every reference)
        order_file = os.path.join(os.path.dirname(os.path.abspath(self.args.schedule)), 'order.tsv')
        names = mapSchedule.read_order(order_file) if os.path.exists(order_file) else None
        cell_of = {c['reference']: c for c in self.cells}
        references = [(cell_of[n]['map_cmd'], cell_of[n]['index_name']) for n in mapSchedule.in_order(self.cells, names)]
        self.dbs = pd.DataFrame([(cmd, self.index_of[name], 'T') for cmd, name in references],
                                columns=['map_cmd', 'db_path', 'retain'])     # map.log "databases:", as before 0.4.0

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
        # one command per map_cmd, with the placeholders of the templates before 0.4.0 (compare_reads.py fills them in)
        for cell in self.cells:
            prefix = '{outprefix}/alignments'
            unal = f'{prefix}/unal/{{sample}}.{{dbprefix}}.fq' + ('.gz' if cell['aligner'] == 'bowtie2' else '')
            line = f"{cell['map_cmd']} : " + mapSchedule.command(
                cell, '{genome}', '{fastq_file}', unal, f'{prefix}/{{sample}}.{{dbprefix}}.sam',
                f'{prefix}/logs/{{sample}}.{{dbprefix}}.log', '{threads}')
            if not any(c.startswith(f"{cell['map_cmd']} : ") for c in self.logdata['commands']):
                self.logdata['commands'].append(line)

        ## get bowtie and bowtie2 versions
        bowtie_version = os.popen('bowtie --version').read().split('\n')[0].split(' ')[2]
        self.logdata['tools'].append(f'bowtie={bowtie_version}')
        bowtie2_version = os.popen('bowtie2 --version').read().split('\n')[0].split(' ')[2]
        self.logdata['tools'].append(f'bowtie2={bowtie2_version}')
        samtools_version = os.popen('samtools --version').read().split('\n')[0].split(' ')[1]
        self.logdata['tools'].append(f'samtools={samtools_version}')

        return
    
    
    def sam_to_bam(self, sam_file, filter_unaligned=True, pass_tag=None):
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
        if pass_tag is not None:             # with --passes: XP:i:<pass> on every mapped record
            import pysam
            with pysam.AlignmentFile(sam_file, 'r') as sam, pysam.AlignmentFile(bam_file, 'wb', template=sam) as bam:
                for record in sam:
                    if record.is_unmapped and filter_unaligned:
                        continue
                    if not record.is_unmapped:
                        record.set_tag('XP', pass_tag, 'i')
                    bam.write(record)
            os.remove(sam_file)
            return bam_file
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
        Maps one sample's reads through the cells of the schedule (tables/schedule.tsv, mapSchedule.py), in order;
        the reads that do not align in a cell go on to the next one.
        '''
        threads = self.args.threads
        out = self.args.output

        # Check if the final BAM file already exists
        final_bam = f'{out}/alignments/{sample}.merged.bam'
        if os.path.exists(final_bam):
            print(f'{sample} already has all alignments. Skipping...')
            return

        fastq = self.fastqDict[sample]

        # FASTQ input only (the collapsed FASTA input was removed in release 0.4.0)
        if not fastq.endswith(('fq.gz', 'fastq.gz', 'fq', 'fastq')):
            print(f'{fastq} is not a valid fastq file. Exiting...')
            sys.exit(1)

        bams, done = [], []
        for i, cell in enumerate(self.cells):
            index = self.index_of[cell['index_name']]
            name = f'{sample}.p{cell["pass"]}.{cell["index_name"]}'        # release 0.4.0: the pass in every name
            prefix = f'{out}/alignments/{name}'
            unal = f'{out}/alignments/unal/{name}.fq' + ('.gz' if cell['aligner'] == 'bowtie2' else '')
            is_last = i == len(self.cells) - 1        # the last cell keeps the unaligned reads in its BAM
            if not os.path.exists(fastq):
                print(f'{fastq} does not exist. Exiting...')
                sys.exit(1)
            skipped, outside = 0, None
            if cell['min_length'] != '-' or cell['max_length'] != '-':     # a window: split the input first
                ending = '.fq.gz' if fastq.endswith('.gz') else '.fq'
                inside, outside = f'{prefix}.window{ending}', f'{prefix}.skipped{ending}'
                _, skipped = split_by_length(fastq, inside, outside, cell['min_length'], cell['max_length'])
                fastq = inside
            ran = has_reads(fastq)
            if ran:
                print(f'Mapping reads for {sample} to {index} (pass {cell["pass"]})')
                run(mapSchedule.command(cell, index, fastq, unal, f'{prefix}.sam',
                                        f'{out}/alignments/logs/{name}.log', threads))
                bams.append(self.sam_to_bam(f'{prefix}.sam', filter_unaligned=not is_last,
                                            pass_tag=int(cell['pass']) if self.args.xp else None))
                if not os.path.exists(unal):         # every read aligned: the next cell gets an empty input
                    open_fastq(unal, 'wt').close()
            else:                                    # no read to align: the cell is not run
                print(f'No reads left for {sample} at {index} (pass {cell["pass"]}): not run')
                open_fastq(unal, 'wt').close()
            if outside is not None and skipped:
                if is_last:                          # the leftovers of the last cell: unmapped records in the BAM
                    bams.append(self.write_unmapped(outside, f'{prefix}.skipped.bam', cell['aligner'], bams))
                else:
                    append_fastq(outside, unal)
            done.append((cell['pass'], cell['index_name'], int(ran), skipped))
            fastq = unal
        with open(f'{out}/alignments/logs/{sample}.cells.tsv', 'w') as f:       # read by qcCollect.py
            f.write('pass\tindex_name\tran\tskipped_window\n')
            f.write(''.join(f'{p}\t{n}\t{r}\t{s}\n' for p, n, r, s in done))

        if not bams:                                 # no read reached any cell: an empty BAM (header only)
            bams.append(self.write_unmapped(os.devnull, f'{out}/alignments/{sample}.empty.bam', 'bowtie2', []))
        # Merge the BAM files and clean up
        self.merge_bams(sample, bams)
        for bam in bams:
            try:
                os.remove(bam)
            except FileNotFoundError:          # with one cell, merge_bams moved it
                pass

        ## delete all unaligned fastq files and the window files
        run(f'rm -f {out}/alignments/unal/{sample}.*.fq* {out}/alignments/{sample}.*.window.fq* '
            f'{out}/alignments/{sample}.*.skipped.fq*')

        return

    def write_unmapped(self, fastq, bam_file, aligner, bams):
        '''The reads of fastq as unmapped records in bam_file (the header of an earlier BAM of the sample, if any).'''
        import pysam
        header = pysam.AlignmentFile(bams[-1], 'rb').header.to_dict() if bams else {'HD': {'VN': '1.0'}}
        with pysam.AlignmentFile(bam_file, 'wb', header=header) as bam:
            for name, sequence, quality in fastq_records(fastq):
                bam.write(unmapped_record(bam.header, name, sequence, quality, aligner))
        return bam_file

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
    ap.add_argument('-d', '--dbs', required=True,
                    help="the references folder's dbpaths file (map_cmd,index,T): where the indexes are")
    ap.add_argument('--schedule', required=True,
                    help='tables/schedule.tsv of the run (mapSchedule.py): the cells to run, in order')
    ap.add_argument('--xp', action='store_true', help='tag every mapped record with XP:i:<pass> (runs with --passes)')
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