# This script computes the rank score 100 - (Ri/N) * 100, based on Anders Christiansen et al., Scientific Reports, Nature, 2015

import Tkinter
import tkFileDialog
import os
import pandas as pd
import string

root3=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root3.withdraw()
root3.update()
file_path = tkFileDialog.askopenfilename(title='Please select the file that need to be clustered')
root3.destroy()

file_name = os.path.basename(file_path).split('.')[0]
dir_name = os.path.dirname(file_path)

## import the test file that need to be compared ##
if file_path.endswith('.csv'):
    test = pd.read_csv(file_path, header=0)
if file_path.endswith('.xlsx') or file_path.endswith('.xls'):
    test = pd.read_excel(file_path, header=0)

###################################################################################################
## Compute the Rank Score, based on Anders Christiansen et al., Scientific Reports, Nature, 2015 ##
###################################################################################################
print('Calculating the Rank Score...')

idx_column = test.columns.str.contains('No.') # select those columns that contain 'No. ' in the column name
idx_column_name = test.columns[idx_column]

for name in idx_column_name:
    # new_name = 'Score_' + name
    temp_rank = test.loc[:,name].sort_values(ascending=False).dropna().rank(method='dense', ascending=False) - 1 # most frequent peptide is ranked #0, second most is #1 ...
    temp_score = 100 - temp_rank/(len(temp_rank.unique())-1) * 100 # Compute the rank score here. Number of unique rank - 1, -1 is to make the rank score for the last rank to be 0
    test.loc[temp_score.index,'Score_'+name] = temp_score

#####################################################
## Prepare to color code columns and save the file ##
#####################################################
# Create a Pandas Excel writer using XlsxWriter as the engine.
print('Saving File...')
excel_file = dir_name + "/RankScore_" + file_name + ".xlsx"
sheet_name = 'Sheet1'

writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')
test.to_excel(writer, sheet_name=sheet_name, index=False, header=True)
workbook = writer.book
worksheet = writer.sheets[sheet_name]

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
cols = test.columns.tolist()
cond_cols_search = ['No.','Percentage','SNratio','COUNT'] # columns with these words will be conditionally formatted
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
    worksheet.conditional_format(colrange, {'type': '3_color_scale'}) # conditional_format('A2:A602', {'type': '3_color_scale'})

writer.save()
print('DONE!')