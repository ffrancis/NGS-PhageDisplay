# NGS-PhageDisplay
Two main sections:

1. [Speed Up Antibody Drug Discovery Process via Machine Learning.](https://github.com/ccchang0111/NGS-PhageDisplay/blob/master/ML_PhageDisplay-Demo.ipynb)

2. Also provided set of scripts for analyzing phage display NGS data:
  
  a. Seperate (demultiplex) samples from a single fastq file into multiple .csv files.
  
  b. Translate the DNA sequences into protein sequences and further clean the data. 
  
  c. Calculate the signal to noise (S/N) ratio for each sample(experiment).
  
  d. Compare S/N between multiple samples.
  
  e. Calculate the ranking score for the compared data from Step d.
  
  g. Cluster the sequences via alignment and then visulaize the result.

-[Visualize 'Cluster Scores' across different conditions](https://plot.ly/~ccchang0111/148/?share_key=opqa7axxTHYZYexL0LhbRx)

**_Voila!_** You get your hits! But that is not enough. We can do more by applying machine learning to the massive NGS data and identify hidden patterns in each experiment.

[Speed Up Antibody Drug Discovery Process via Machine Learning.](https://github.com/ccchang0111/NGS-PhageDisplay/blob/master/ML_PhageDisplay-Demo.ipynb
