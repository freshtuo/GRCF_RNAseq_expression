#!/usr/bin/env python
# summarize_fpkm.py
# summarize the cufflink FPKM gene/isoform tables
# 

import sys
import math
import gzip
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from contextlib import redirect_stdout

# sample ids
samples = snakemake.config['read1'].keys()

# input FPKM gene tables
fpkm_gene_files = snakemake.input['gene']

# input FPKM isoform tables
fpkm_isoform_files = snakemake.input['isoform']

# reference annotation file
refann_file = snakemake.config['refann']

# sample grouping information
grpdic = snakemake.config["group"]

# output merged FPKM gene table
merged_fpkm_gene_file = str(snakemake.output['gene'])

# output merged FPKM isoform table
merged_fpkm_isoform_file = str(snakemake.output['isoform'])

# output box plot
out_box_file = str(snakemake.output['box'])

# log file
log_file = snakemake.log[0]

# functions
def assign_group(grpdic):
    """create a dictionary to assign each sample to a group"""
    group = {}
    for grp in grpdic:
        for sample in grpdic[grp]:
            group[sample] = grp
    return group

def merge_fpkm(fpkm_files):
    """merge FPKM tables from multiple samples"""
    # columns to keep
    mycols = ['tracking_id','gene_id','gene_short_name','locus','FPKM']
    # read in fpkm tables
    fpkm_list = []
    for fpkm_file in fpkm_files:
        # extract sample id
        ##results/project/cufflinks_fpkm/DMSO_1/genes.fpkm_tracking
        sid = fpkm_file.split('/')[-2]
        print('Load {}: '.format(fpkm_file), end='')
        # read into dataframe
        fpkm = pd.read_table(fpkm_file, header=0, sep='\t', low_memory=False)
        # rename columns and add to list
        fpkm_list.append(fpkm[mycols].rename(columns={'FPKM':sid}))
        print(fpkm_list[-1].shape)
    # merge fpkm tables
    merged_table = fpkm_list[0]
    for fpkm in fpkm_list[1:]:
        merged_table = pd.merge(merged_table, fpkm, how='outer', on=['tracking_id','gene_id','gene_short_name','locus'])
    print('merged table: {}'.format(merged_table.shape))
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

def write_fpkm(fpkm_table, out_file, sheet_name='gene'):
    """write FPKM table to an excel file"""
    with pd.ExcelWriter(out_file, engine='xlsxwriter') as writer:
        fpkm_table.to_excel(writer, sheet_name=sheet_name, index=False)
        # apply format
        nrows, ncols = fpkm_table.shape
        worksheet = writer.sheets[sheet_name]
        worksheet.set_column(0, 0, 16)
        worksheet.set_column(1, 1, 16)
        worksheet.set_column(2, 2, 24)
        worksheet.set_column(3, 3, 17)
        for k in range(4, ncols-4):
            worksheet.set_column(k, k, 15)
        worksheet.set_column(ncols-4, ncols-4, 16.5)
        worksheet.set_column(ncols-3, ncols-3, 10)
        worksheet.set_column(ncols-2, ncols-2, 18.5)
        worksheet.set_column(ncols-1, ncols-1, 50)
        worksheet.autofilter(0, 0, nrows-1, ncols-1)
        worksheet.freeze_panes(1, 4)
        ##writer.save()

def fpkm_boxplot(fpkm_table, grpdic, out_file):
    # convert to long format for plotting
    fpkm_long = pd.melt(fpkm_table, id_vars=['tracking_id','gene_id','gene_short_name','locus','chromosome','strand','gene_type','description'], var_name='sample', value_name='fpkm')
    # log transform on FPKM
    fpkm_long['logfpkm'] = fpkm_long['fpkm'].apply(lambda x: math.log2(x+1))
    # add a column indicating grouping
    group = assign_group(grpdic)
    fpkm_long['group'] = fpkm_long['sample'].apply(lambda x: group[x])
    # make box plot
    sns.set(rc = {'figure.figsize':(10.5,8.5)})
    sns.set_style('white')
    ##sns_plot = sns.violinplot(data=fpkm_long, x='sample', y='logfpkm', hue='group', scale='width', palette="Pastel1")
    sns_plot = sns.boxplot(x='sample', y='logfpkm', hue='group', data=fpkm_long, palette="Pastel1", linewidth=1.5)
    sns_plot.set(xlabel='Samples', ylabel='Log2(FPKM+1)')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    ##plt.show()
    fig = sns_plot.get_figure()
    fig.savefig(out_file, dpi=300)

with open(log_file, 'w') as flog:
    with redirect_stdout(flog):
        # merge FPKM gene tables
        merged_fpkm_gene = merge_fpkm(fpkm_gene_files)
        # annotate FPKM gene tables
        annotated_merged_fpkm_gene = annotate_genes(merged_fpkm_gene, refann_file)
        # save to file
        write_fpkm(annotated_merged_fpkm_gene, merged_fpkm_gene_file, 'gene')
        # make box plot on per-gene FPKM
        fpkm_boxplot(annotated_merged_fpkm_gene, grpdic, out_box_file)
        # merge FPKM isoform tables
        merged_fpkm_isoform = merge_fpkm(fpkm_isoform_files)
        # annotate FPKM isoform tables
        annotated_merged_fpkm_isoform = annotate_genes(merged_fpkm_isoform, refann_file)
        # save to file
        write_fpkm(annotated_merged_fpkm_isoform, merged_fpkm_isoform_file, 'isoform')
        print('Complete!')

