# v2-2 Andrew Chang on 10/18/2016: for each column name (percentage and no. peptide..), added the filename at the end of the name. This will help clarify ambiguity for the next step
# v2-1 Andrew Chang on 4/15/2016: added new column that calculates the percentage of each seq. relative to the non-cleaned data. 
# v2 Andrew Chang on 4/11/2016

# Purpose: clean the .csv files in the COUNT folder; remove asterisk, remove frame shift, and add template info

# Procedure:
# 1. User select all .csv files that need to be cleaned
# 2. User input whether to remove sequences with stop codon (but not the amber stop)
# 3. User input whether to remove sequences with frame shift. If yes, then provide the correct last two  amino acids sequences
# 4. User input whether to map each sequence to its template (only Reflexion and CKP are considered here, but the list can be expanded) 

import Tkinter
import tkFileDialog
import os
import pandas as pd

root6=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root6.withdraw()
root6.update()
file_paths = tkFileDialog.askopenfilenames(title='Please select the COUNT files that need to be cleaned (multiple files ok)')
root6.destroy()
rmv_stop = raw_input("Remove sequences contain stop codon (Y or N): ")
rmv_fshift = raw_input("Remove frame-shift (Y or N): ")
if (rmv_fshift == 'Y') or (rmv_fshift == 'y'):
    fshift = raw_input("Last two amino acids when there is no frame-shift (ex: SG for Relfx): ")
frm_search = raw_input("Map the framework (Y or N): ")
if (frm_search == 'Y') or (frm_search =='y'):
    lib = raw_input("Please specify the library: Relexion = 1, CKP = 2 ")
    

for file_path in file_paths:  
    file_name = os.path.basename(file_path).split('.')[0]
    test = pd.read_csv(file_path, header=0)
    test['frmw']='NA' # initialize the new column for storing framework info
    test_totcnt = test['No. peptide'].sum() # total number of sequences in the file before any cleaning
    test['Percentage'] = (test['No. peptide']/test_totcnt)*100 # add a new column of 'normalized counts' or 'percentage' based on the total seq. read in that sample
    
    if (rmv_stop == 'Y') or (rmv_stop =='y'):
        idx=test['Peptide'].str.contains('\*')
        test = test.drop(test[idx].index)
    
    if (rmv_fshift == 'Y') or (rmv_fshift =='y'):
        test = test[test['Peptide'].str[-2:] == fshift]
    
    if (frm_search == 'Y') or (frm_search =='y'):
        Reflexion_dict={}
        CKP_dict = {}
        
        if int(lib) == 1:
            for key, value in Reflexion_dict.iteritems():
                index = test['Peptide'].str.contains(value)
                test.loc[index,'frmw']=key
        
        if int(lib) == 2:
            for key, value in CKP_dict.iteritems():
                index = test['Peptide'].str.contains(value)
                test.loc[index,'frmw']=key
    
    test['L2'] = 'NA' # initialize new column for
    test['L2'] = test['Peptide'].str[-2:] # extract the last two amino acid (check if there is any frame-shift) into a new column
    test = test[['No. peptide','Percentage','Peptide','frmw','L2']]
    test.columns = ['No. peptide_'+file_name,'Percentage_'+file_name, 'Peptide', 'frmw', 'L2']
    test.to_csv(os.path.dirname(file_path) + "/Cleaned_" + file_name + ".csv", index=False) 