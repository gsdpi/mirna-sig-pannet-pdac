import os
import gzip
import shutil

import numpy as np
import pandas as pd
import urllib.request

##################################################################
# 1 ) Source data: TCGA datasets
# Download of gene expression, miRNA and clinical datasets.
# 
# Note: This require gnu utils installed in your system (gunzip)
# this is by default in Linux or MacOS, but not in Windows.
# You can also manually uncompress the data.
#
# 1. a ) Gene expression and miRNA datasets
# Download of gene expression and miRNA datasets from Xena Browser.
# at https://xenabrowser.net/datapages/?hub=https://tcga.xenahubs.net:443
###################################################################


# Create folders for data

try:  
    os.mkdir('downloaded_data')
except OSError as error:  
    print(error)   

try:  
    os.mkdir('datasets')
except OSError as error:  
    print(error)   


# Start the download    
print('downloading Pancreatic Cancer data from xenabrowser ...')
    
urllib.request.urlretrieve('https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/miRNA_HiSeq_gene.gz ','downloaded_data/miRNA_HiSeq_gene.gz')
with gzip.open('downloaded_data/miRNA_HiSeq_gene.gz', 'rb') as f_in:
    with open('downloaded_data/miRNA_HiSeq_gene.csv', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
   
print('download completed')

################################################################
# 2 ) Preprocessing of source data
# Preprocessing of gene expression, miRNA and clinical datasets.
# 
# 2. a ) Gene expression and miRNA datasets
# Gene expression datasets: we remove rows containing the same value in all columns.
# miRNA datasets: we remove rows containing missing values.
# Finally, we merge all datasets into one called 'df_data'.
################################################################

print ('Preparing miRNA data ...')

genes_and_mirna = []

# load csv, delete file and drop rows which contain missing values
df_data = pd.read_csv('downloaded_data/miRNA_HiSeq_gene.csv',delim_whitespace=True)
os.unlink('downloaded_data/miRNA_HiSeq_gene.csv')
df_data.set_index('sample', inplace=True)
df_data.dropna(axis='index',inplace=True)
    
# Fill some additional fields
df_data = df_data.T
df_data['cancer']  = 'PAAD'
df_data['tumor'] = df_data.index

# homogenize type of tumor in a new variable called 'type'
df_data['type'] = 'other'
df_data.loc[df_data.index.str[-2:] == '01', 'type'] = 'primary'
df_data.loc[df_data.index.str[-2:] == '02', 'type'] = 'recurrent'
df_data.loc[(df_data.index.str[-2:] == '06') |
            (df_data.index.str[-2:] == '07'), 'type'] = 'metastatic'
df_data.loc[(df_data.index.str[-2:] == '10') |
            (df_data.index.str[-2:] == '11'), 'type'] = 'sane'



#########################################################
# 3 ) Final dataset
#########################################################

print ('Creating final dataset ...')

# change index to submitter_id
df_data.index = ['-'.join(i.split('-')[:3]) for i in df_data.index]
df_data.index.name = 'submitter_id'

df_data.to_hdf('datasets/df_TCGA.hdf', key='raiz')


