#!/usr/bin/env python
# htseq_count_reads.py
# count reads per gene using HTSeq-count
# 

import sys
import subprocess
import gzip
from contextlib import redirect_stdout

# input bam file
bamfile = snakemake.input['bam']

# output counts file
countfile = snakemake.output['counts']

# library type
libtype = snakemake.params['libtype']

# minimal quality
minaqual = snakemake.params['minaqual']

# reference gtf file
refgtf = snakemake.params['refgtf']

# threads
threads = snakemake.threads

# memory
mem = snakemake.params['mem']

# log file
htseq_log_file = snakemake.log['htseq']

# strandness information
strand = ''
if libtype == 'fr-firststrand':
    strand = 'reverse'
elif libtype == 'fr-unstranded':
    strand = 'no'
elif libtype == 'fr-secondstrand':
    strand = 'yes'
else:
    print('Unknown strandness: {}'.format(libtype))
    sys.exit(1)

# run HTSeq-count to generate per-gene read counts
command = 'samtools sort -n -@ {} -m {} -T results/htseq_count_reads/sort.tmp.{} -O sam {} |htseq-count --format=sam --stranded={} --minaqual={} --type=exon --idattr=gene_id --mode=union - {} >{}'.format(threads,mem,snakemake.wildcards,bamfile,strand,minaqual,refgtf,countfile)

with gzip.open(htseq_log_file, 'wt') as flog:
    with redirect_stdout(flog):
        print(command)
        sys.stdout.flush()
        res = subprocess.run(args=command, shell=True, capture_output=True, text=True)
        print(res.stdout)
        sys.stdout.flush()
        print(res.stderr)
        print('Complete!')

