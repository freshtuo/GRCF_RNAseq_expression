#!/usr/bin/env python
# star_align.py
# align the trimmed reads to the reference genome using STAR
# 

import sys
import subprocess
import os.path
from contextlib import redirect_stdout

# sequencing type (default: single-end)
seqtype = 'SE'

# input fastq files
fq_files = snakemake.input['read1']
if len(snakemake.input) > 1:
    seqtype = 'PE'
    fq_files = '{} {}'.format(snakemake.input['read1'], snakemake.input['read2'])

# wildcards.sample
sample = snakemake.wildcards[1]

# output bam file
bamfile = snakemake.output['bam']

# output bam prefix
bamprefix= bamfile.replace('Aligned.sortedByCoord.out.bam','')

# output bam folder
bamdir = os.path.dirname(bamfile)

# library type
libtype = snakemake.config['libtype']

# indexed reference genome folder
staridx = snakemake.config["staridx"]

# parameters for STAR
seedLmax = snakemake.config["seedLmax"]
multimapNmax = snakemake.config["multimapNmax"]
overhangMin = snakemake.config["overhangMin"]
DBoverhangMin = snakemake.config["DBoverhangMin"]

# threads
threads = snakemake.threads

# log file
star_log_file = snakemake.log['bam']
samtools_log_file = snakemake.log['bai']
err_log_file = snakemake.log['err']

# compatibility for Cufflinks
flag = ''
if libtype == 'fr-unstranded':
    flag = '--outSAMstrandField intronMotif'

# run STAR to generate alignment
star_command = 'STAR --runThreadN {} --genomeDir {} --readFilesIn {} --readFilesCommand zcat --outFileNamePrefix {} {} --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax {} --alignSJoverhangMin {} --alignSJDBoverhangMin {} --seedSearchStartLmax {} --twopassMode Basic >{} 2>&1'.format(threads,staridx,fq_files,bamprefix,flag,multimapNmax,overhangMin,DBoverhangMin,seedLmax,star_log_file)

gzip_command = 'gzip {}/{}.Unmapped.out.mate*'.format(bamdir,sample)

samtools_command = 'samtools index -@ {} {} >{} 2>&1'.format(threads,bamfile,samtools_log_file)

# functions
def run_command(command):
    print(command)
    res = subprocess.run(args=command, shell=True, capture_output=True, text=True)
    print(res.stdout)
    sys.stdout.flush()
    print(res.stderr)
    sys.stdout.flush()

# run
with open(err_log_file, 'w') as flog:
    with redirect_stdout(flog):
        print('This is a {} run.'.format(seqtype))
        # run STAR
        run_command(star_command)
        # zip unmapped reads
        run_command(gzip_command)
        # index bam file
        run_command(samtools_command)
        print('Complete!')

