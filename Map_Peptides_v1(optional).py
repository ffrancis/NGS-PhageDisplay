# v2 Andrew Chang on 4/30/2016

# Purpose: compare sequences between two files

import Tkinter
import tkFileDialog
import os
import pandas as pd
import string
"""
    template sequence file: has to be .csv with column names: Peptide and frame work
"""
root3=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root3.withdraw()
root3.update()
file_path_1 = tkFileDialog.askopenfilename(title='Please select file 1',filetypes=[('Excel','*.csv *.xlsx *.xls')])
file_path_2 = tkFileDialog.askopenfilename(title='Please select file 2',filetypes=[('Excel','*.csv *.xlsx *.xls')])
root3.destroy()

file_name = os.path.basename(file_path_1).split('.')[0]

# import the test file that need to be compared
if file_path_1.endswith('.csv'):
    test = pd.read_csv(file_path_1, header=0)
if file_path_1.endswith('.xlsx') or file_path_1.endswith('.xls'):
    test = pd.read_excel(file_path_1, header=0)

if file_path_2.endswith('.csv'):
    par = pd.read_csv(file_path_2, header=0)
if file_path_2.endswith('.xlsx') or file_path_2.endswith('.xls'):
    par = pd.read_excel(file_path_2, header=0)

overlap=pd.merge(test,par, how='outer', on='Peptide') 
    
#csv_file = os.path.dirname(file_path_1) + "/Seq_Compared_" + file_name + ".csv"
#overlap.to_csv(csv_file, index=False)

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
print("prepare to save the file ...")

cols = overlap.columns.tolist()
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

excel_file = os.path.dirname(file_path_1) + "/Seq_Compared_" + file_name + ".xlsx"
writer = pd.ExcelWriter(excel_file)
overlap.to_excel(writer,'Sheet1',index=False, header=True)
workbook = writer.book
worksheet = writer.sheets['Sheet1']

num_useq = len(overlap)
# Apply a conditional format to the cell range.
for col in cols_excel:
    colrange = col + '2:' + col +str(num_useq+1) 
    worksheet.conditional_format(colrange, {'type': '3_color_scale'}) # conditional_format('A2:A602', {'type': '3_color_scale'})

writer.save()