# Created by Thomas Hengstl
# Modified by Andrew Chang on 4/11/2016

# Purpose: This script is for demultiplexing barcoded sequences into its own .csv file
#
# Procedure:
# 1. User select the one file for IonTorrent .fastq file and one for the barcode file (excel file)
# 2. This script will look for barcode in both directions within the first X nucleotides and then save all DNA sequences with the same barcode into a .csv file
# 3. A folder called "DEMULTPLEXED" will be created inside the .fastq file directory. All the demultiplexed .csv files will be stored here. 

from pandas import DataFrame as df
from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys
import Tkinter
import tkFileDialog

##### Function Definition #####

def def_demultiplex(rootdir, fastq_file_path, barcode_file_path):
    """
    def_demultiplex produces seperate .fastq DNA files based on their barcodes from multiplexed fastq files and writes the DNA sequences into seperate .csv files
    
    Files need to be in rootdir:
    - Multiplexed.fastq    Fastq sequencing file with sequences containing different barcodes
    - Barcodes.csv         Csv file containing all relevant barcodes and their names in comma seperated format (csv windows style): name, barcode
    
    Output files: 
    - Csv Files with DNA sequences named according to the what's in the barcode file
    - Statistics file called Demultiplexing.txt with the name of the file and the number of demultiplexed sequences
    """
    print "def_demultiplex"
    
    #Create barcode dictionary, statisitcs count dictionary, target fastq files
    barcode_file = open(barcode_file_path,"rU") #format: name , barcode. "r" is for reading file. "U" reads new lines across dif platform from Mindows and Mac
    barcode_dict = {} 
    count_dict = {}
    
    if not os.path.exists(rootdir + "/DEMULTIPLEXED"): os.makedirs(rootdir + "/DEMULTIPLEXED") # create a folder 'DEMULTIPLEXED' if it does not exist in the directory
    
    for barcode_line in barcode_file: #fill the barcode dictionary with filenames (barcodes) and the barcode DNA sequence
        barcode_template = barcode_line.replace("\r","").replace("\n","").split(",") #remove '\r and \n' and seperate the Barcode name from the Barcode sequence
        barcode_dict[barcode_template[0]] = barcode_template[1] #create a dictionary {'barcode_name': 'barcode_sequence',...} 
        count_dict[barcode_template[0]] = 0 #create another dictionary to count the number of reads in that barcode file
        target_file = open(rootdir + "/DEMULTIPLEXED/" + str(barcode_template[0]) + ".csv", "w") # careate a folder 'DEMULTIPLEXED' and generate empty .csv files with barcode names
        target_file.close()
    barcode_file.close()
    
    #Demultiplexing
    print 'It is running in the background, please be patient for couple minutes.\nIf this takes more than 1 hour, then you might consider getting a better computer :)'
    DNA_seq_all = [] # create an empty list that will later store all the DNA sequences in fastq file
    fastq_file = open(fastq_file_path,"rU") 
    for record in SeqIO.parse(fastq_file,"fastq"): # read line-by-line of the.fastq file.
        DNA_seq = str(record.seq) # and then extract only the sequence string
        DNA_seq_all.append(DNA_seq) # and then append each sequence string onto a new list called 'DNA_seq_all' which will then be converted to a dataframe for parsing
        
    fastq_file.close()
    DNA_seq_all = df(DNA_seq_all) # make a DataFrame
    
    counter = 0
    for key, value in barcode_dict.iteritems(): # for each barcode name, get the barcode sequence
        index = DNA_seq_all[0].str.contains(value) | DNA_seq_all[0].str.contains(str(Seq(value).reverse_complement())) # find all the DNA sequences in the fastq file that match this particular barcode sequence or its reverse compliment
        barcode_file = DNA_seq_all[index] # pull out all the DNA sequences that match the barcode seq and rev. compliment
        barcode_file.to_csv(rootdir + "/DEMULTIPLEXED/" + str(key) + ".csv", index=False, header=False) # saving file, overwrite existing file 
        count_dict[key] = barcode_file.size
    
        counter = counter +1 # counter for counting the progress of demuliplexing
        length = len(barcode_dict) # number of barcodes => how many loops have to be run
        # printing the progress of the Demultiplex
        sys.stdout.write('\r %d out of %d' %(counter, length))
        sys.stdout.flush()
    
    #Demultiplexing Statistics
    if not os.path.exists(rootdir + "/STATISTICS"): os.makedirs(rootdir + "/STATISTICS")
    stat = open(rootdir + "/STATISTICS/" + "Statistics-demultiplex.txt","w")
    stat.write("\t")
    stat.write("Demultiplexed")
    stat.write("\n")
    for k, v in count_dict.iteritems():
        stat.write(str(k))
        stat.write("\t")
        stat.write(str(v))
        stat.write("\n")
    stat.close()

#########################################
root=Tkinter.Tk() # if you are using Canopy, please go to 'Preferences > Python > PyLab Backend: Inline(SVG)'
root.withdraw()
root.update()
fastq_file_path = tkFileDialog.askopenfilename(title='Please select the fastq file')
barcode_file_path = tkFileDialog.askopenfilename(title='Please select the Barcode file')
root.destroy()
rootdir = os.path.dirname(fastq_file_path)
def_demultiplex(rootdir, fastq_file_path, barcode_file_path) # run the demuliplex.py