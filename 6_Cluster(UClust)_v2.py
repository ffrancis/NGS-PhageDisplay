# v2 on 10/31/2016: include cluster stats.csv. For each cluster, it will summarize the 1. # of seq in each cluster, 2. alignment %ID, 3 ranking score
# Andrew Chang on 9/27/2016

# Note: This requires installing Usearch package (use Terminal for installation) 
# Usearch webpage: http://drive5.com/usearch/
# how to use python to call terminal command? http://stackoverflow.com/questions/89228/calling-an-external-command-in-python
# Download USEARCH file, http://www.drive5.com/usearch/download.html
# open the terminal, navigate to the folder that has the usearch file 
# rename it at the terminal using: mv usearch6.0.98_i86linux32 usearch 
# move the renamed file into the $PATH: /usr/local/bin by typing: cp usearch /usr/local/bin
# allow permission by typing: chmod +x /usr/bin/usearch
# if you want to clean 'Alignment' files using terminal: find . -name '*Alignment*' -delete

# Run this script using Terminal: Go to Tool => Canopy Terminal: 
# when in Terminal, navigate to the folder that contains this script => type: python 5_Cluster_peptides....py

# Purpose: from a .csv or .xlsx file that contain 'Peptide' column, cluster them into groups (UCLUST uses BLOSUM62 scoring matrix)

# Procedure:
# 1. select a .csv or .xlsx file that has the column name 'Peptide' (like the cleaned or analyzed file)
# 2. output the result as a .fasta file
# 3. perform Usearch clust_fast to cluster the .fasta file
# 4. from clustered result (.clstr file), find each cluster and mark corresponding rows of the input excel file

import Tkinter
import tkFileDialog
import os
import pandas as pd
import string
import time
#from numpy import mean

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

##############################################
## Convert selected file into .fasta format ##
##############################################
root3=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root3.withdraw()
root3.update()
file_path = tkFileDialog.askopenfilename(title='Please select the file that need to be clustered')
root3.destroy()

file_name = os.path.basename(file_path).split('.')[0]

## User Inputs ##
threshold = raw_input("Input the sequence alignment threshold (0.85): ")
PerOrRank = raw_input("Generate Cluster Heatmap based on either (1)Percentage or (2)RankScore: ")

## import the test file that need to be compared ##
if file_path.endswith('.csv'):
    test = pd.read_csv(file_path, header=0)
if file_path.endswith('.xlsx') or file_path.endswith('.xls'):
    test = pd.read_excel(file_path, header=0)

my_seqs=[]
# Identifying which data seq is the parent seq and mark it templt 'N'
for i, peptide in enumerate(test["Peptide"]):
    my_seqs.append(SeqRecord(Seq(str(peptide),IUPAC.protein), id = str(i)))
    #print(i)

fasta_file = file_path.split('.')[0] +'.fasta'

with open(fasta_file, "w") as handle:
  SeqIO.write(my_seqs, handle, "fasta")   

#########################################
## Cluster the fasta file using CD-HIT ##
#########################################
cluster_file = os.path.dirname(fasta_file) + "/Alignment_" + file_name 

# CD-HIT clustering
infile = fasta_file
outfile1 = cluster_file + '_centroids.fasta'
outfile2 = cluster_file + '_clusters.fasta'
outfile3 = cluster_file + '_clusters.uc'
BLOSUM62 = '/Volumes/PET2/Andrew/NGS/NGS_codes/BLOSUM62.txt' # the default already uses BLOSUM62 scoring, but just to make sure it does that, force it to the scoring matrix.
BLOSUM80 = '/Volumes/PET2/Andrew/NGS/NGS_codes/BLOSUM80.txt' # the default already uses BLOSUM80 scoring, but just to make sure it does that, force it to the scoring matrix.

#cmd = "usearch -cluster_fast " + infile + " -sort length " + " -id " + str(threshold) + " -centroids " + outfile1 + " -uc " + outfile3
cmd = "usearch -cluster_fast " + infile + " -sort length " + " -id " + str(threshold) + " -centroids " + outfile1 + " -msaout " + outfile2 + " -uc " + outfile3
#cmd = "usearch -cluster_fast " + infile + " -sort length " + " -id " + str(threshold) + " -centroids " + outfile1 + " -matrix " + BLOSUM80 + " -msaout " + outfile2 + " -uc " + outfile3

# cmd = "usearch -cluster_fast " + infile + " -sort length " + " -id " + str(threshold) + " -centroids " + outfile1 + " -idsuffix 30 -maxdiffs 6 -query_cov 0.9 " + " -matrix " + BLOSUM62 + " -msaout " + outfile2 + " -uc " + outfile3 # for scf28 (last 30 charc are the same, and only allow 4 positions of substitution)
# cmd = "usearch -cluster_fast " + infile + " -sort length " + " -id " + str(threshold) + " -centroids " + outfile1 + " -idsuffix 30 -idprefix 13 " + " -matrix " + BLOSUM62 + " -msaout " + outfile2 + " -uc " + outfile3 # for scf28 (last 30 charc are the same, and only allow 4 positions of substitution)

print(cmd)
os.system(cmd)

#######################################################################
## Cluster the input excel file based on the info from CD-HIT output ##
#######################################################################
print('Clustering the input file...')
time.sleep(3) # pause for 3 seconds, make sure the CD_HIT has genrated .clstr file

dir_name = os.path.dirname(file_path)
# cluster_file_path = dir_name + "/Alignment_" + file_name +".txt.clstr"

#filename = '/Volumes/Andrew/NGS/Example_Data/Clustering_Test/Peptide_LIb_test/COUNT/ANALYSIS/usearch/clusters.uc'
data = pd.read_table(outfile3, 
    sep='\t',
    names=['Record Type','Cluster_ID','Seq Length','H_%id','H_strand','na','na','note','Label','Cluster rep'],
    na_values=['*'])

data['Label'] = data['Label'].str.split(' ').str.get(0).astype(int)
data['Cluster rep'] = data['Cluster rep'].str.split(' ').str.get(0) # this is still string not integer. need to convert it during accessing.
cluster_num = data['Cluster_ID'].unique()
test['Cluster'] = 'NA' # initialize a new column for Cluster

# pull up the # of seq for each cluster from the .uc file and save it to a new dataframe, cluster_stat
cluster_stat = data.loc[data['Record Type'] == 'C',['Cluster_ID','Seq Length']] 
cluster_stat['%ID'] = 'NA'

# Cluster members
for value in cluster_num:
    idx = data[(data["Cluster_ID"] == value) & (data['Record Type']!='C')]['Label'].tolist()
    test.loc[idx,['Cluster']] = value # add the cluster id to the 'Cluster' column
    cluster_stat.loc[cluster_stat['Cluster_ID']==value,'%ID'] = data.loc[(data['Record Type'] == 'H')&(data['Cluster_ID']==value),'H_%id'].mean() # add average %ID to cluster_stat
    

##############################################
## Sum over the Rank Score for each cluster ##
##############################################
print('Generating the Heatmap...')

Cluster_Stats_dict = {}
if int(PerOrRank) == 1:
    idx_column = test.columns.str.contains('No.') - test.columns.str.contains('Score')
    total_peptide = test.loc[:,idx_column].sum()
if int(PerOrRank) == 2:
    idx_column = test.columns.str.contains('Score')

for cluster in cluster_num:
    idx = test['Cluster'] == cluster # select only a particular cluster row
    test_temp = test.loc[idx,idx_column] # select only a particular cluster row
    if int(PerOrRank) == 1:
        fraction_peptide = test_temp.sum()/total_peptide
    if int(PerOrRank) == 2:
        fraction_peptide = test_temp.sum()
    Cluster_Stats_dict[cluster] = fraction_peptide

df = pd.DataFrame(Cluster_Stats_dict).T
cluster_stat = cluster_stat.set_index('Cluster_ID') # re-define the index based on the cluster #
cluster_stat = pd.concat([cluster_stat,df], axis=1)

#####################################################
## Prepare to color code columns and save the file ##
#####################################################
# Create a Pandas Excel writer using XlsxWriter as the engine.
print('Saving File...')
excel_file1 = dir_name + "/Clustered_" + file_name + ".xlsx"
excel_file2 = dir_name + "/Clustered_Stats_" + file_name + ".xlsx"
sheet_name = 'Sheet1'

writer1 = pd.ExcelWriter(excel_file1, engine='xlsxwriter')
writer2 = pd.ExcelWriter(excel_file2, engine='xlsxwriter')
test.to_excel(writer1, sheet_name=sheet_name, index=False, header=True)
cluster_stat.to_excel(writer2, sheet_name=sheet_name, index=True, header=True)
workbook1 = writer1.book
worksheet1 = writer1.sheets[sheet_name]
workbook2 = writer2.book
worksheet2 = writer2.sheets[sheet_name]
######################################################
## convert the column number to excel column letter ##
######################################################
def colnum_string(num):
    title = ''
    alist = string.uppercase
    while num:
        mod = (num-1) % 26
        num = int((num - mod) / 26)  
        title += alist[mod]
    return title[::-1]

## For TEST file ##
# color code (conditional formatting) those columns with numbers and percentates
cols = test.columns.tolist()
cond_cols_search = ['No.','Percentage','SN','COUNT','%','S/N','COUNT'] # columns with these words will be conditionally formatted
cond_cols=[]
for cond in cond_cols_search:
    s = [s for s in cols if cond in s]
    cond_cols.extend(s)

# convert the columns that need to be formatting into excel column letter
cols_excel = []    
for i,name in enumerate(cols):
    if name in cond_cols:
        cols_excel.append(colnum_string(i+1))

# Apply a conditional format to the cell range.
num_useq = len(test) # number of unique sequences
for col in cols_excel:
    colrange = col + '2:' + col +str(num_useq+1) 
    worksheet1.conditional_format(colrange, {'type': '3_color_scale'}) # conditional_format('A2:A602', {'type': '3_color_scale'})

writer1.save()

## For Stats file ##
# color code (conditional formatting) those columns with numbers and percentates
cols2 = cluster_stat.columns.tolist()
#cond_cols_search = ['No.','Percentage','SN','COUNT','%','S/N','COUNT'] # columns with these words will be conditionally formatted
cond_cols2=[]
for cond in cond_cols_search:
    s = [s for s in cols2 if cond in s]
    cond_cols2.extend(s)

# convert the columns that need to be formatting into excel column letter
cols_excel2 = []    
for i,name in enumerate(cols2):
    if name in cond_cols2:
        cols_excel2.append(colnum_string(i+2)) # i+2 is because cluster_stat's index becomes the first column, so the total column # +1

# Apply a conditional format to the cell range.
num_useq = len(cluster_stat) # number of unique sequences
for col in cols_excel2:
    colrange = col + '2:' + col +str(num_useq+1) 
    worksheet2.conditional_format(colrange, {'type': '3_color_scale'}) # conditional_format('A2:A602', {'type': '3_color_scale'})

writer2.save()
print('DONE!')