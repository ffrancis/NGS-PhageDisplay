# Andrew Chang on 4/11/2016

# Purpose: 
# 1. Extract DNA sequence within the primer recognition sequences and then save data to => DNA folder 
# 2. Translate DNA sequences into protein sequences => PROTEIN folder
# 3. Count the copy number of each protein sequence => COUNT folder
#
# Procedure:
# 1. User select .csv files (multiple files) of interest within the DEMULTIPLEX folder
# 2. User input the forward and reverse primer sequences
# 3. User input whether to use amber stop codon or not
# 4. User set the threshold for the minimum length of the protein sequence, above the threshold the sequence will be saved 

from Bio.Seq import Seq
import os
import sys
import collections
import Tkinter
import tkFileDialog

##### Function Definitions #####
def def_extract(rootdir, file_name_list, fwd_primer, rev_primer, limit):
    """
    Extracts the DNA sequence inserts based on 3' and 5' primer recognition sequences
    """

    print "def_extract"
    
    #Extraction primer dictionary (NEEDS TO BE DELETED BEFORE PUBLICATION)
    stat_extract_dict = {}
    rc_primer  = str(Seq(rev_primer).reverse_complement())
    
    #Opening input files and creating output files
    if not os.path.exists(rootdir + "/DNA"): os.makedirs(rootdir + "/DNA")
    for filenames in file_name_list:
        print filenames
        sequences = open(rootdir + "/DEMULTIPLEXED/" + filenames, "rU")
        output = open(rootdir + "/DNA/" + filenames.strip(".csv") + "-DNA.csv","w")
        stat_extract_dict[filenames] = 0
    
    #Searching primers
        for i,line in enumerate(sequences): #enumerate: adds a counter to each element (i is the line number in this case)
            start = 0
            end = 0
            if line.find(rc_primer) > 0:
                line = Seq(line).reverse_complement()
            start = line.find(fwd_primer)
            end = line.find(rev_primer)
            extract = str(line[start + len(fwd_primer):end])
    
    #Extracting inserts
            if start >= 0 and end >= 0:
            #if len(extract) >= limit_low:
                output.write (extract)
                output.write("\n")
                stat_extract_dict[filenames] += 1
    
            #if i == limit:
                #break
        output.close()
        sequences.close()
    
    #Extraction Statisitcs
    if not os.path.exists(rootdir + "/STATISTICS"): os.makedirs(rootdir + "/STATISTICS")
    stat = open(rootdir + "/STATISTICS/" + "Statistics-extract.txt","w")
    stat.write("\t DNA extracted \n")
    
    for k, v in stat_extract_dict.iteritems():
        stat.write(str(k))
        stat.write("\t")
        stat.write(str(v))
        stat.write("\n")
    stat.close()

def def_translate(rootdir, file_name_list, limit, amber_stop):
    """
    Translates all files in the DNA directory into protein sequences. The output files contain the DNA and the PROTEIN sequence in the order Protein, DNA
    
    rootdir and file_name_list have the same fuction as in def_extraction and are also provided by the main function
    The limit defines the lowest allowed length of the translation product, if the translated sequence is shorter than the limit it will not be considered
    """
    print "def_translate"
    print "limit", limit
    
    base2aa =  {"AAA" : "K", "AAC" : "N", "AAG" : "K", "AAT" : "N", "ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T", "AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S", "ATA" : "I", "ATC" : "I", "ATG" : "M",  "ATT" : "I", "CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H", "CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P", "CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R", "CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L", "GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D", "GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A", "GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G", "GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V", "TAA" : "*", "TAC" : "Y", "TAG" : "*", "TAT" : "Y", "TCA" : "S", "TCC" : "S", "TCG" : "S",  "TCT" : "S", "TGA" : "*", "TGC" : "C", "TGG" : "W", "TGT" : "C", "TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT":"F"}
    
    if (amber_stop == 'T') or (amber_stop == 't'):
        base2aa['TAG']='O' # show 'O' instead of '*' if amber_stop is 'T' true
    
    def seqtranslation(nt):
        nt=nt.upper()
        aa=''
        while len(nt)>=3:
            codon=nt[0:3]
       	    nt=nt[3:]
       	    aa=aa+base2aa[codon]
        return aa
        
    #Opening input files and creating output files
    stat_translate_dict = {}
    if not os.path.exists(rootdir + "/PROTEIN"): os.makedirs(rootdir + "/PROTEIN")
    
    counter=0
    for filenames in file_name_list:
        sequence = open(rootdir + "/DNA/" + filenames.strip(".csv") + "-DNA.csv", "rU")
        output = open(rootdir + "/PROTEIN/" + filenames.strip(".csv") + "-PROTEIN.csv", "w")
        stat_translate_dict[filenames] = 0
        translation = ""
        
        #Translation        
        for line in sequence:
            line = line.replace("\r","").replace("\n","").split(",")
            try:
                translation=seqtranslation(line[0])
                		                
                if len(translation) >= int(limit): # save the protein sequence only when its length is longer than the set limit
                    output.write(str(translation))
                    output.write(",")
                    output.write(str(line[0]))
                    output.write("\n")
                    stat_translate_dict[filenames] += 1
            except:
                continue
        output.close()
        sequence.close()
        counter = counter+1
        length=len(file_name_list)
        sys.stdout.write('\r %d out of %d' %(counter, length))
        sys.stdout.flush()
        
    #Translation Statistics
    if not os.path.exists(rootdir + "/STATISTICS"): os.makedirs(rootdir + "/STATISTICS")
    stat = open(rootdir + "/STATISTICS/" + "Statistics-translation.txt","w")
    stat.write("\t Translated \n")
    for k, v in stat_translate_dict.iteritems():
        stat.write(str(k).strip(".csv"))
        stat.write("\t")
        stat.write(str(v))
        stat.write("\n")
    stat.close()

def def_count(rootdir, file_name_list, limit):
    """
    COUNTS how often a specific amino acid sequence is found in the a given PROTEIN file. First a list with unique sequnce is created form the the protein file.
    Next, it is counted how often a specific sequence occurs. The total number and the protein sequence are then wrtten into the count file.
    
    rootdir and file_name_list have the same fuction as in def_extraction and are also provided by the main function
    The limit defines how often an amino acid sequence has to occur in order to be included in the final count file
    """ 

    print "def_Protein_count"

    #Opening input files and creating output files
    if not os.path.exists(rootdir + "/COUNT"): os.makedirs(rootdir + "/COUNT")
    stat_count_dict = {}
    for filenames in file_name_list:
        sequence = open(rootdir + "/PROTEIN/" + filenames.strip(".csv") + "-PROTEIN.csv", "rU")
        output = open(rootdir + "/COUNT/" + filenames.strip(".csv") + "-COUNT.csv", "w")
        
        #Writing headline in cunt file and creating lists required for counting
        output.write("No. peptide,Peptide\n")
        stat_count_dict[filenames] = 0
        seq_master_file = []
        temp_master_file = []

        #Create temp_masterfile with all amino acid sequences found in the protein file
        for line in sequence: 
            line = line.replace("\r","").replace("\n","").split(",")
            temp_master_file.append(line[0])  # protein sequences in a particular Barcode exp
        sequence.close()

        #Create seq_masterfile that contains only unique sequnces
        for x in set(temp_master_file):
            seq_master_file.append(x) # only the unique seq is collected and saved in seq_master_file
        stat_count_dict[filenames] = len(seq_master_file) # how many unique sequences

        #Counting
        cnt = collections.Counter()  # or can write as: cnt = collections.Counter(temp_master_file). the 'Counter' function already consider the uniqueness 
        for word in temp_master_file:
            cnt[word] += 1

        #Writing count dictionary to file for every sequences that ocurs more often than the limit (limit is set to 1)
        for k,v in cnt.iteritems():
            if v >= int(limit):
                output.write(str(v))
                output.write(",")
                output.write(str(k))
                output.write("\n")
        output.close()

    #Counting Statistics
    if not os.path.exists(rootdir + "/STATISTICS"): os.makedirs(rootdir + "/STATISTICS")
    stat = open(rootdir + "/STATISTICS/" + "Statistics-counting.txt","w")
    stat.write("\t Unique_sequences \n")

    for k, v in stat_count_dict.iteritems():
        stat.write(str(k).strip(".csv"))
        stat.write("\t")
        stat.write(str(v))
        stat.write("\n")
    stat.close()
#########################################
root8=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root8.withdraw()
root8.update()
demutiplex_file_paths = tkFileDialog.askopenfilenames(title='Please select all DEMULTIPLEXED files that need to be analyzed')
root8.destroy()

fwd_primer = raw_input("Please enter sequence of forward extraction primer ")
rev_primer = raw_input("Please enter sequence of reverse extraction primer ")
amber_stop = raw_input("Use amber stop or not? (T/F) ")
# extraction_limit = raw_input("Please enter how many sequences should be analyzed ")
min_length = raw_input("Please enter the length limit for the shortest protein sequence ")
# min_copy_number = raw_input("Please enter the lower copy number limit (5 is recommend) ")

file_name_list = [os.path.basename(j) for j in demutiplex_file_paths]
rootdir = os.path.abspath(os.path.join(os.path.dirname(demutiplex_file_paths[0]),'..'))

def_extract(rootdir, file_name_list, fwd_primer, rev_primer, 0) # extraction limit is set to be ignored
def_translate(rootdir,file_name_list, min_length, amber_stop) 
def_count(rootdir, file_name_list, 1) # minimum copy number is set to be 1. 