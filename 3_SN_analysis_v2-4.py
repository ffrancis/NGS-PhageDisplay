# v2-4 AC 10/14/2016: take data from 2_xx_v2-2, and remove rename the column name
# v2-3 AC 8/10/2016: instead of listing just the test seq, list all seq for both test and control
# v2-2 Andrew Chang on 6/6/2016: 
# 1. S/N ratio is now the normalized S divided by normalized N. 
# 2. The default percentage for the controls that do not exist in the target sample is set to be the 'minimum percentage value found in the control file'
# v2-1 Andrew Chang on 4/15/2016: remove the percentage calculation, because it is now done in 2_Clean_COUNTcsv_file_v2-1. This script does percentage subtraction only. 
# v2 Andrew Chang on 4/11/2016

# Purpose: calculate the percentage of the test and control(s) sequences, S/N ratio (ratio of absolute count) and S-N (difference of percentage)

# Procedure:
# 1. User set a default number for sequences that do not exist in the control sample
# 2. User select the panned (test) file
# 3. User select the control files 
# 4. Script will create a new folder "ANALYSIS" and color-code all the numbers and save output as .xlss file

import Tkinter
import tkFileDialog
import os
import pandas as pd
import string

#default_count = raw_input("Please set a default number for sequences that do not exist in control sample (1) ") # assign a default count if the seq does not appear in the control sample
root0=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root0.withdraw()
root0.update()
test_file_path = tkFileDialog.askopenfilename(title='Please select the panned file') 
control_file_paths = tkFileDialog.askopenfilenames(title='Please select the control files (multiple files are ok)')
root0.destroy()

test_file_name = os.path.basename(test_file_path).strip('.csv')
test = pd.read_csv(test_file_path, header=0) # test.sort('No. peptide', ascending=False) # put the test protein seq into dataframe

control_file_names = [os.path.basename(j).strip('.csv') for j in control_file_paths] # get filenames without the path (ex of list comprehension) 
controls = []

# import control files
for control_file_path in control_file_paths:
    control = pd.read_csv(control_file_path, header=0) # put all the control .csv files into 'controls'
    controls.append(control)

overlap=test
# compare seq between test and control(s)
for i, control in enumerate(controls):   
    overlap=pd.merge(overlap,control, how='outer', on='Peptide') 

# re-name those column names that contain 'No. peptide' and 'Percentage' to their file name 
# and also fill-in NA rows with: 0 for the test column, and min value for controls
idx = overlap.columns.str.contains('Percentage') # find the column name that contains 'Percentage'
idx_no = overlap.columns.str.contains('No.') # find the column name that contains 'No.'

test_percentage_name = overlap.columns[idx][0]
test_no_name = overlap.columns[idx_no][0]
control_percentage_names = overlap.columns.values[idx][1:]

idx2 = overlap[test_percentage_name].isnull() # find all the NaN cells
overlap.loc[idx2,test_percentage_name] = 0.0 # fill those NaN cells with 0

# calculate the S/N and S-N (normalized, or percentage)      
for control_percentage_name in control_percentage_names:
    idx3 = overlap[control_percentage_name].isnull() # find all the NaN cells in the controls
    overlap.loc[idx3,control_percentage_name] = overlap[control_percentage_name].min() # fill those NaN cells with the min value
    overlap['S/N_'+ control_percentage_name] = overlap[test_percentage_name]/overlap[control_percentage_name] 
    #test['Percentage_S-N_' + control_file_name] = test['Percentage'] - test['Percentage_' + str(i)] # calculate the difference of normalized count

test = overlap.sort_values(test_no_name, ascending=False)

###################################################
# Prepare to color code columns and save the file #
###################################################
# Create a Pandas Excel writer using XlsxWriter as the engine.
if not os.path.exists(os.path.dirname(test_file_path) + "/ANALYSIS"): os.makedirs(os.path.dirname(test_file_path) + "/ANALYSIS")
excel_file = os.path.dirname(test_file_path) + "/ANALYSIS/" + 'Analyzed_' + test_file_name + ".xlsx"
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