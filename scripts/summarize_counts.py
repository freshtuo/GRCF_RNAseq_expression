#!/usr/bin/env python
# summarize_counts.py
# summarize the HTSeq-count counts result
# 

import sys
import pandas as pd
from contextlib import redirect_stdout

# sample ids
samples = snakemake.config['read1'].keys()

# input counts tables
counts_files = snakemake.input['counts']

# reference annotation file
refann_file = snakemake.config['refann']

# output merged counts table
merged_counts_file = snakemake.output['counts']

# output merged other info table (e.g. '__no_feature')
merged_other_file = snakemake.output['other']

# log file
log_file = snakemake.log[0]

# functions
def merge_counts(counts_files):
    """merge counts tables from multiple samples"""
    # read in counts tables
    counts_list = []
    for counts_file in counts_files:
        # extract sample id
        ##results/htseq_count_reads/DMSO_1.genes.HTSeq.count
        sid = counts_file.split('/')[2].replace('.genes.HTSeq.count','')
        print('Load {}: '.format(counts_file), end='')
        # read into dataframe
        counts = pd.read_table(counts_file, header=None, sep='\t', names=['gene_id',sid], low_memory=False)
        # add to list
        counts_list.append(counts)
        print(counts_list[-1].shape)
    # merge counts tables
    merged_table = counts_list[0]
    for counts in counts_list[1:]:
        merged_table = pd.merge(merged_table, counts, how='outer', on='gene_id')
    print('merged table: {}'.format(merged_table.shape))
    print(merged_table.head())
    return merged_table

def annotate_genes(merged_table, refann_file):
    """add annotations to genes"""
    # read in reference annotations
    refann = pd.read_table(refann_file, header=0, sep='\t', low_memory=False)
    # rename columns
    refann.rename(columns={'Gene stable ID':'gene_id','Chromosome/scaffold name':'chromosome','Gene type':'gene_type','Gene description':'description'}, inplace=True)
    # fix strand info
    refann['strand'] = refann['Strand'].apply(lambda x: '+' if x > 0 else '-')
    # dedup (use the first entry in case multiple annotations are found for a give ensembl id)
    refann.drop_duplicates(subset=['gene_id'], keep='first', inplace=True)
    # merge annotations into the FPKM table
    return pd.merge(merged_table, refann[['gene_id','chromosome','strand','gene_type','description']], how='left', on='gene_id')

def write_counts(counts_table, out_counts_file, out_other_file):
    """write counts table to a text file"""
    # entry representing a gene?
    mark = ~counts_table['gene_id'].str.contains('^__')
    # gene counts
    counts_table[mark].to_csv(out_counts_file, sep='\t', index=False)
    # other information
    counts_table[~mark].drop(['chromosome','strand','gene_type','description'], axis=1).to_csv(out_other_file, sep='\t', index=False)

with open(log_file, 'w') as flog:
    with redirect_stdout(flog):
        # merge counts tables
        merged_counts = merge_counts(counts_files)
        # annotate counts table
        annotated_merged_counts = annotate_genes(merged_counts, refann_file)
        # save to file
        write_counts(annotated_merged_counts, merged_counts_file, merged_other_file)
        print('Complete!')

