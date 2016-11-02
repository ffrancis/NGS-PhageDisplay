# v2-1 Andrew Chang on 4/30/2016: modified the saved filename

# Purpose: compare all the input files from different rounds of sorting and find overplapping sequences and color code all the numbers

# Procedure:
# 1. User select analyzed files from multiple rounds that need to be compared
# 2. Script will color-code all the numbers and save output as .xlss file

import Tkinter
import tkFileDialog
import os
import pandas as pd
import string

root2=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root2.withdraw()
root2.update()
file_paths = tkFileDialog.askopenfilenames(title='Please select the files in orders')
root2.destroy()

# if the file order is wrong, please manually set it up like the following:
# file_paths = ['/folder/file_1.xlsx','folder/file_2.xlsx','folder/file_3.xlsx',...etc]

file_names = [os.path.basename(j).strip('.xlsx') for j in file_paths] # get filenames without the path (ex of list comprehension) 
# extract all the interesting columns (Peptide name, lib, S/N, normalized S-N) of each selected file into 'files'
files=[]
for i,file_path in enumerate(file_paths):
    file_temp = pd.read_excel(file_path, header=0)
    file_temp.rename(columns={'No. peptide': 'No. peptide_'+str(i), 'Percentage': 'Percentage_'+str(i)}, inplace=True) # change the column name to differentiate from other files that may have the same name
    files.append(file_temp) # put all the cleaned files into 'files'

# compare sequences across multiple rounds (across all the selected files) 
df_final = reduce(lambda left,right: pd.merge(left,right,on='Peptide',how='outer'), files)

# Prepare to color code columns and save the file
# convert the column number to excel column letter
def colnum_string(num):
    title = ''
    alist = string.uppercase
    while num:
        mod = (num-1) % 26
        num = int((num - mod) / 26)  
        title += alist[mod]
    return title[::-1]

# color code (conditional formatting) those columns with numbers and percentates
cols = df_final.columns.tolist()
cond_cols_search = ['No.','Percentage','SN','COUNT','%','S/N'] # columns with these words will be conditionally formatted
cond_cols=[]
for cond in cond_cols_search:
    s = [s for s in cols if cond in s]
    cond_cols.extend(s)

# convert the columns that need to be formatting into excel column letter
cols_excel = []    
for i,name in enumerate(cols):
    if name in cond_cols:
        cols_excel.append(colnum_string(i+1))

excel_file = os.path.dirname(file_paths[0]) + "/Enrichment_" + file_names[0] +".xlsx"
writer = pd.ExcelWriter(excel_file)
df_final.to_excel(writer,'Sheet1',index=False, header=True)
workbook = writer.book
worksheet = writer.sheets['Sheet1']

num_useq = len(df_final)
# Apply a conditional format to the cell range.
for col in cols_excel:
    colrange = col + '2:' + col +str(num_useq+1) 
    worksheet.conditional_format(colrange, {'type': '3_color_scale'}) # conditional_format('A2:A602', {'type': '3_color_scale'})

writer.save()