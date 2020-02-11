# DMS-FastQ-processing
A set of scripts for analyzing deep sequencing data for deep mutational scanning (DMS)

# System requirements
The scripts require a python 3 installation along with the following packages:  
1.pandas  
2.numpy  
3.matplotlib  
4.seaborn  

# Installation
The scripts require no installation. To use the script, set up a working python 3 environment by either  
1. Downloading python 3 as a standalone installation (https://www.python.org/downloads/) and using pip-install command to install the required packages  
  
or  
  
2. Downloading python 3 as a package manager such as Anaconda (https://www.anaconda.com/distribution/), which has all the required packages installed.

# Example on using the scripts
## Setting up the environment
Clone or download this GitHub project to your local machine, which will create a folder with all the scripts inside.

Download the example data sets and extract them into a new folder named 'inputs', in the same directory as all the scripts. The file structure should look something like this:  
DMS-FastQ-processing/  
  
* inputs/  
    * vim2_library_128amp_37c/  
        * AMP_HIGH_37_G1_rep1_S4_R1_001.fastq.gz  
        * ...  
* merge_reads.py  
* ...  

With python installed, each of the scripts can be run directly from the folder where the scripts are located. By opening a command line or terminal in the folder and using the python command to run the script:  
`python name_of_script.py`

## Preparing sample information
Sequencing read files can be organized in the "sample_reference.xlsx" excel workbook with 2 named spreadsheets:

1.samples  
Each row in the spreadsheet corresponds to a single sample with the following information:  
* path - Filepath to the FastQ file. Can be relative or absolute path, with relative paths specified relative to where the scripts are being run.  
* label - Label to use for identifying the general experimental condition of the sample. Usually the name of an antibiotic, nosel (no selection) or wt (no selection, wt DNA).  
* conc - Concentration of antibiotic used in screen.  
* temp - Temperature of screen.  
* group - Gene mutational library group.  
* rep - Replicate number.  
* read - Forward or reverse read. By illumina convention, 1 is the forward read and 2 is the reverse read.  
  
2.trim_info  
Each row in the spreadsheet gives information on the appropriate trimming behavior for each of the mutational library groups.   Each library group indicates a consecutive segment of the target gene that contains mutations. Numbers indicate number of nucleotide positions.  
* gene - Label to indicate the gene. Multiple genes can be placed in the sheet.  
* group - Gene mutational library group.  
* fwd trim - Number of nucleotides to remove from the front of the forward read.  
* rev trim - Number of nucleotides to remove from the front of the reverse read.  
* offset - Keeps track of the starting nucleotide position of the library group. Used for obtaining gene codon and residue positions.  
* length - The number of nucleotides to take from the read after the starting nucleotide. Nucleotides after this length are removed.  
