import pandas as pd
import numpy as np
import os
import sys
import re
import warnings
warnings.filterwarnings('ignore')

# Read cosmic/TCGA data 
cosmic_GRCh37 = pd.read_table("../raw_data/CosmicCompleteGeneExpression_GRCh37.tsv")

# Read GTF file 
ens69_gtf = pd.read_table("../raw_data/Homo_sapiens.GRCh37.69.gtf", header = None, sep = '\t')
ens69_gtf.columns = ['chr', 'rna_type', 'transcript_region', 'start', 'end', '5', 'strand', '7', 'transcript_info']

# SUMMARY STATISTICS 
## Non-normal distribution stats
### Median Z 
cosmic_GRCh37_median = cosmic_GRCh37.groupby('GENE_NAME')[['GENE_NAME', 'Z_SCORE']].median()
cosmic_GRCh37_median['GENE_NAME'] = cosmic_GRCh37_median.index
### IQR Z
cosmic_GRCh37_iqr = cosmic_GRCh37.groupby('GENE_NAME')[['GENE_NAME', 'Z_SCORE']].quantile([.25, .5, .75]).unstack()
cosmic_GRCh37_iqr['GENE_NAME'] = cosmic_GRCh37_iqr.index
cols = ['Q1', 'Q2', 'Q3', 'gene_name']
cosmic_GRCh37_iqr.columns = cols
cosmic_GRCh37_iqr = cosmic_GRCh37_iqr.reset_index(drop = True)
cosmic_GRCh37_iqr['IQR'] = cosmic_GRCh37_iqr['Q3'] - cosmic_GRCh37_iqr['Q1']
cosmic_GRCh37_iqr['lower_1.5IQR'] = cosmic_GRCh37_iqr['Q1'] - 1.5*cosmic_GRCh37_iqr['IQR']
cosmic_GRCh37_iqr['upper_1.5IQR'] = cosmic_GRCh37_iqr['Q3'] + 1.5*cosmic_GRCh37_iqr['IQR']

## Normal stats - Mean, sd, etc
cosmic_GRCh37_summary = cosmic_GRCh37.groupby('GENE_NAME')[['GENE_NAME', 'Z_SCORE']].describe()
cosmic_GRCh37_summary['gene_name'] = cosmic_GRCh37_summary.index
cosmic_GRCh37_summary = cosmic_GRCh37_summary.reset_index(drop = True)
cols= ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max', 'gene_name']
cosmic_GRCh37_summary.columns = cols

## Join non-normal and normal stats 
cosmic_GRCh37_summary_full = pd.merge(cosmic_GRCh37_summary,cosmic_GRCh37_iqr, how = 'inner') 

# cosmic_GRCh37_summary_full.to_csv("../tables_output/cosmic_GRCh37_summary_full.tsv", sep = '\t', index = False)

# TRANSFORM Z_SCORE DISTRIBUTION VALUES  
## Find minimum (most negative) value
min_median_Z_SCORE = min(cosmic_GRCh37_summary_full['Q2'])
## Make all median (Q2) Z_SCORE values positive 
cosmic_GRCh37_summary_full['Q2_Z_positive'] = cosmic_GRCh37_summary_full['Q2'] - (min_median_Z_SCORE)
## Find total Q2_Z_positive (to normalize)
total_median_Z_SCORE = sum(cosmic_GRCh37_summary_full['Q2_Z_positive'])

# FUNCTIONS 
# Define function tidy_split()
def tidy_split(df, column, sep=';', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.
    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df.reset_index(drop=True)

info_data = ens69_gtf['transcript_info'].str.split(';', expand=True)
info_data = info_data.iloc[:, :7]
info_data.columns = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 
                     'gene_biotype', 'transcript_name', 'exon_id']
info_data = info_data.fillna('exon id "none"')
info_data['exon_id'] = info_data['exon_id'].replace('', 'exon id "none"')

# Define update_columns()
def update_columns(data, column):
    data[column] = data[column].map(lambda x: x.split('"')[1])
    return data
    
for i in info_data.columns:
    info_data = update_columns(info_data, i)
    
    
# PROCESS GTF FILE 
## Ens69 GTF data processing 
ens69_gtf_expand = pd.concat([ens69_gtf, info_data], axis=1)
ens69_gtf_expand['chr'] = ens69_gtf_expand['chr'].astype(str)

chromosomes_lst = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 
                   '10', '11', '12', '13', '14', '15', '16', '17', 
                   '18', '19', '20', '21', 'X', 'Y']

ens69_gtf_chr = ens69_gtf_expand[ens69_gtf_expand['chr'].isin(chromosomes_lst )]


# Get real start and end of transcripts (not exon level)
## START
transcripts_min = ens69_gtf_chr.groupby(['transcript_id'])[['start']].first().drop_duplicates()
transcripts_min = transcripts_min.dropna()
transcripts_min['transcript_id'] = transcripts_min.index
transcripts_min = transcripts_min.dropna()
transcripts_min = transcripts_min.reset_index(drop = True)
## END 
transcripts_max = ens69_gtf_chr.groupby(['transcript_id'])[['end']].max().drop_duplicates()
transcripts_max = transcripts_max.dropna()
transcripts_max['transcript_id'] = transcripts_max.index 
transcripts_max = transcripts_max.dropna()
transcripts_max = transcripts_max.reset_index(drop  = True)

## JOIN START-END
transcripts_start_end = pd.merge(transcripts_min, transcripts_max, how= 'left')
transcripts_start_end = transcripts_start_end.rename(columns={'start': 'tx_start', 'end': 'tx_end'})


## SUBSET GTF 
ens69_gtf_chr_subset = ens69_gtf_chr[['chr', 'gene_name', 'strand', 'transcript_id', 'gene_id', 'gene_biotype']]

## EXPAND GTF START-END DF
transcripts_start_end_genes = pd.merge(transcripts_start_end, ens69_gtf_chr_subset, how = 'left')
transcripts_start_end_genes = transcripts_start_end_genes.dropna()


## Ensure transcript start and end are integers
transcripts_start_end_genes['tx_end'] = transcripts_start_end_genes['tx_end'].astype(int)
transcripts_start_end_genes['tx_start'] = transcripts_start_end_genes['tx_start'].astype(int)
transcripts_start_end_genes['tx_len'] = transcripts_start_end_genes['tx_end'] - transcripts_start_end_genes['tx_start']

## Sort data, reset index 
transcripts_start_end_genes = transcripts_start_end_genes.drop_duplicates()
transcripts_start_end_genes = transcripts_start_end_genes.sort_values(['chr','tx_start',  'gene_id' , 'gene_name'])
transcripts_start_end_genes.reset_index(inplace = True, drop = True)


# Compute proportion of the full gene length for each transcript 
transcripts_start_end_genes['sum_tx_len'] = transcripts_start_end_genes.groupby('gene_name')['tx_len'].transform('sum')
transcripts_start_end_genes['tx_proportion'] = transcripts_start_end_genes['tx_len'] / transcripts_start_end_genes['sum_tx_len']

# Rename column expression 
cosmic_GRCh37_summary_full = cosmic_GRCh37_summary_full.rename(columns = {'GENE_NAME':'gene_name'})

# Create expressed fraction and expressed number columns 
## Count the number of transcripts and assume equal contribution of RNA from each transcript not accounting for length
cosmic_GRCh37_summary_full['count_transcripts'] = cosmic_GRCh37_summary_full.groupby('gene_name')['gene_name'].transform('count')

## * PRO Expressed Fraction col *
cosmic_GRCh37_summary_full['Expressed Fraction'] = cosmic_GRCh37_summary_full['Q2_Z_positive']/(total_median_Z_SCORE) 
# Adjust for number of transcripts 
cosmic_GRCh37_summary_full['Expressed Fraction - Transcript'] = cosmic_GRCh37_summary_full['Expressed Fraction']/cosmic_GRCh37_summary_full['count_transcripts']

## * PRO Expressed Number col *
cosmic_GRCh37_summary_full['Expressed Number'] = cosmic_GRCh37_summary_full['Expressed Fraction'] * 100000
# Adjust for number of transcripts 
cosmic_GRCh37_summary_full['Expressed Number - Transcript'] = cosmic_GRCh37_summary_full['Expressed Number']/cosmic_GRCh37_summary_full['count_transcripts']

# Round to integer for Expressed Number 
cosmic_GRCh37_summary_full = cosmic_GRCh37_summary_full.round({'Expressed Number - Transcript':0})

## MERGE EXPRESSION + GTF INFORMATION 
exp_gtf = pd.merge(cosmic_GRCh37_summary_full, transcripts_start_end_genes, how = 'inner')

## Count the number of transcripts and assume equal contribution of RNA from each transcript not accounting for length
exp_gtf['count_transcripts'] = exp_gtf.groupby('gene_name')['gene_name'].transform('count')

## Strand system switch 
exp_gtf['strand_WC'] = pd.Series(exp_gtf['strand']).str.replace('+', 'W')
exp_gtf['strand_WC'] = pd.Series(exp_gtf['strand']).str.replace('-', 'C')

## * PRO Locus col *
exp_gtf['Locus'] = exp_gtf['chr'].astype(str) + ":" + exp_gtf['tx_start'].astype(str) + "-" + exp_gtf['tx_end'].astype(str) + exp_gtf['strand_WC'].astype(str)

## Function to determine if transcript CDS (coding) or NC (noncoding)
def determine_coding(x):
    if x=='protein_coding':
        return('CDS')
    else:
        return('NC')

## * PRO Coding col *     
exp_gtf['Coding'] = exp_gtf['gene_biotype'].map(determine_coding)    

## * PRO Length col * 
exp_gtf['Length'] = exp_gtf['tx_end'] - exp_gtf['tx_start']

# FINAL PRO 
PRO_final = exp_gtf[['Locus', 'transcript_id', 'Coding', 'Length', 'Expressed Fraction - Transcript', 'Expressed Number - Transcript']]

## Remove all transcripts < 10 nt 
PRO_final_clean = PRO_final[PRO_final['Length'] >= 10]

PRO_final_clean['Expressed Number - Transcript'] = (PRO_final_clean['Expressed Number - Transcript']).astype(int)
# Compute total number of molecules 
total_N_MOLECULES = sum(PRO_final_clean['Expressed Number - Transcript']).astype(int)

## SAVE TO FILE 
pro_file_name = 'TCGA_dataEXP_median_N_MOLECULES_' + total_N_MOLECULES.astype(str) + ".PRO"

PRO_final_clean.to_csv(('../tables_output/' + pro_file_name), sep = '\t', index = False, header = None)

#####################################

## SUBSET CHR12 and CHR21
exp_gtf_chr12_chr21 = exp_gtf[exp_gtf['chr'].isin(['12', '21'])] 
PRO_final_clean_chr12_chr21 = exp_gtf_chr12_chr21[['Locus', 'transcript_id', 'Coding', 'Length', 'Expressed Fraction - Transcript', 'Expressed Number - Transcript']]


## THEN RE-COMPUTE TOTAL READS 
total_N_MOLECULES = sum(PRO_final_clean_chr12_chr21['Expressed Number - Transcript']).astype(int)

## PROPORTION
PRO_final_clean_chr12_chr21['Expressed Fraction - Transcript'] = PRO_final_clean_chr12_chr21['Expressed Number - Transcript']/total_N_MOLECULES

## INTEGER
PRO_final_clean_chr12_chr21['Expressed Number - Transcript'] = PRO_final_clean_chr12_chr21['Expressed Number - Transcript'].astype(int)
pro_file_name = 'TCGA_dataEXP_median_N_MOLECULES_' + total_N_MOLECULES.astype(str) + ".PRO"

# ISSUE. The TCGA-COSMIC data does not have information on all the genes included in the Ensembl69 GTF
    # Therefore, there are 5 thousands transcripts missing in the TCGA-COSMIC PRO expression file

# Solution: Merge both files and populate the transcripts/genes not present in TCGA-COSMIC with zero

## Read "default PRO" file that Flux Simulator puts together based on GTF
pro_man = pd.read_table("../tables_output/12_21_test.pro", sep = '\t', header = None)

pro_man.column = ['Locus', 'transcript_id', 'Coding', 'Length', 'Expressed Fraction - Transcript', 'Expressed Number - Transcript']
pro_man_v2 = pro_man[[0, 1, 2, 3, 4, 5]]
pro_man_v2.columns = ['Locus', 'transcript_id', 'Coding', 'Length', 'Expressed Fraction - Transcript', 'Expressed Number - Transcript']

PRO_MLT = PRO_final_clean_chr12_chr21[['transcript_id', 'Expressed Fraction - Transcript', 'Expressed Number - Transcript']]

only_man = [x for x in pro_man_v2['transcript_id'].unique() if x not in PRO_MLT['transcript_id'].unique()]

# Filter for only in the "default PRO" flux simulator 
ens69_gtf_only_man = pro_man_v2[pro_man_v2['transcript_id'].isin(only_man)]

# MERGE the "default" PRO file and the PRO file generated from TCGA/COSMIC data 
MERGE_PRO_TCGA = pd.concat([ens69_gtf_only_man, PRO_MLT])
final_PRO_TCGA = final_PRO_TCGA.drop_duplicates()

total_N_MOLECULES = sum(final_PRO_TCGA['Expressed Number - Transcript'])

final_PRO_TCGA['Expressed Fraction - Transcript'] = final_PRO_TCGA['Expressed Number - Transcript']/total_N_MOLECULES

final_PRO_TCGA.to_csv(('../tables_output/final_PRO_TCGA_expanded_N_MOLECULES_' + total_N_MOLECULES.astype(str) + '.PRO'), sep = '\t', index = False)

# Make a smaller PRO file (not as many reads)
small_final_PRO_TCGA = final_PRO_TCGA
small_final_PRO_TCGA['Expressed Number - Transcript'] = small_final_PRO_TCGA['Expressed Number - Transcript'] / 10
small_final_PRO_TCGA['Expressed Number - Transcript'] = small_final_PRO_TCGA['Expressed Number - Transcript'].astype(int)

total_N_MOLECULES = sum(small_final_PRO_TCGA['Expressed Number - Transcript'])

small_final_PRO_TCGA.to_csv(('../tables_output/final_SMALL_PRO_TCGA_expanded_N_MOLECULES_' + total_N_MOLECULES.astype(str) + '.PRO'), sep = '\t', index = False)

## This file 