
import GEOparse
import pandas as pd

import urllib.request
import os
import gzip
import shutil

# Create folders for data
try:  
    os.mkdir('downloaded_data')
except OSError as error:  
    print(error)   

try:  
    os.mkdir('datasets')
except OSError as error:  
    print(error)   
####################################################
# GSE163031: Formalin-fixed and paraffin-embedded tissue (FFPE) from 25 PDAC surgical specimens and 13 non cancerous pancreatic tissue samples 
####################################################

# Automatic download with get_GEO seems not to work for this set, so we use a direct download link
print('downloading GSE163031 data ...')
urllib.request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163031/soft/GSE163031_family.soft.gz','./downloaded_data/GSE163031_family.soft.gz')

# get_GEO will work with the already downloaeded data.
gse = GEOparse.get_GEO(geo= "GSE163031", destdir="./downloaded_data")

# Arrange data, set type and miRNA
lcancer = [i for i in gse.phenotype_data["title"] if "Cancer" in i]
df2 = gse.phenotype_data[["title"]]
df2=df2[df2['title'].isin(lcancer)]
lcancer = df2.index.to_list()
df = gse.pivot_samples('VALUE')
df = df.T
df['Type']='Normal'
df.loc[lcancer,'Type']='Pancre'
df.columns=df.columns.str.strip('_st')

# Save dataframe
df.to_hdf('datasets/df_GSE163031.hdf', key='raiz')


####################################################
# GSE73367: A cross species and multi-omics (including metabolomics) analysis in pancreatic neuroendocrine tumours (miRNA)
####################################################

gse = GEOparse.get_GEO(geo= "GSE73367", destdir="./downloaded_data")

# For this dataset, the GEOparse library seems not to download the actual
# miRNA expression values, which are in another file...
# We'll do this by hand again.

print('downloading miRNA expression data ...')
urllib.request.urlretrieve('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE73367&format=file&file=GSE73367%5Fprocessed%5Fdata%2Etxt%2Egz','downloaded_data/GSE73367_data.csv.gz')
with gzip.open('downloaded_data/GSE73367_data.csv.gz', 'rb') as f_in:
    with open('downloaded_data/GSE73367_data.csv', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


df = pd.read_csv('downloaded_data/GSE73367_data.csv', sep='\t')
df = df.iloc[: , 1:]
df = df.set_index('miRNA')
df = df.T
df = df.iloc[:,df.columns.str.contains('hsa')]

# Save dataframe
df.to_hdf('datasets/df_GSE73367.hdf', key='raiz')

####################################################
# GSE43796: Solid-pseudopapillary neoplasm of pancreas (SPN), ductal adenocarcinoma (PCA), neuroendocrine tumor (NET) and non-neoplastic pancreas.
# comparison with gene expression of tumors and non-tumors
####################################################

gse = GEOparse.get_GEO(geo= "GSE43796" , destdir="./downloaded_data")

# Get miRNA data from samples of interest (pancreatic)
pivoted_control_samples = gse.pivot_samples('VALUE') 
pivoted_control_samples.head()
df = pivoted_control_samples.T

lpca = [i for i in gse.phenotype_data["title"] if "PCA" in i]

df2 = gse.phenotype_data[["title"]]
df2=df2[df2['title'].isin(lpca)]
lpca = df2.index.to_list()

lnet = [i for i in gse.phenotype_data["title"] if "NET" in i]
df2 = gse.phenotype_data[["title"]]
df2=df2[df2['title'].isin(lnet)]
lnet = df2.index.to_list()

lnormal = [i for i in gse.phenotype_data["title"] if "Normal" in i]
df2 = gse.phenotype_data[["title"]]
df2=df2[df2['title'].isin(lnormal)]
lnormal = df2.index.to_list()

df['Type']='SPN'
df.loc[lpca,'Type']='PDAC'
df.loc[lnet,'Type']='PNET'
df.loc[lnormal,'Type']='Normal'

# Save dataframe
df.to_hdf('datasets/df_GSE43796.hdf', key='raiz')
