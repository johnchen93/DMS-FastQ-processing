# DMS-FastQ-processing
A set of scripts for analyzing deep sequencing data for deep mutational scanning (DMS)

# System requirements
The scripts require a python 3 installation along with the following packages:  
1. pandas  
2. numpy  
3. matplotlib  

# Installation
The scripts require no installation. To use the script, set up a working python 3 environment by either  
1. Downloading python 3 as a standalone installation (https://www.python.org/downloads/) and using pip-install command to install the required packages  
  
or  
  
2. Downloading python 3 as a package manager such as Anaconda (https://www.anaconda.com/distribution/), which has all the required packages installed.

# Example on using the scripts
## Setting up the environment
Clone or download this GitHub project to your local machine, which will create a folder with all the scripts inside.

Download the example data sets from the release page (https://github.com/johnchen93/DMS-FastQ-processing/releases/tag/v1.0) and extract them into a new folder named 'inputs', in the same directory as all the scripts. The file structure should look something like this:  
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

1. samples  
Each row in the spreadsheet corresponds to a single sample with the following information:  
* path - Filepath to the FastQ file. Can be relative or absolute path, with relative paths specified relative to where the scripts are being run.  
* label - Label to use for identifying the general experimental condition of the sample. Usually the name of an antibiotic, nosel (no selection) or wt (no selection, wt DNA).  
* conc - Concentration of antibiotic used in screen.  
* temp - Temperature of screen.  
* group - Gene mutational library group.   
* rep - Replicate number.  
* read - Forward or reverse read. By illumina convention, 1 is the forward read and 2 is the reverse read.  
  
2. trim_info  
Each row in the spreadsheet gives information on the appropriate trimming behavior for each of the mutational library groups.   Each library group indicates a consecutive segment of the target gene that contains mutations. Numbers indicate number of nucleotide positions.  
* gene - Label to indicate the gene. Multiple genes can be placed in the sheet.  
* group - Gene mutational library group.  
* fwd trim - Number of nucleotides to remove from the front of the forward read.  
* rev trim - Number of nucleotides to remove from the front of the reverse read.  
* offset - Keeps track of the starting nucleotide position of the library group. Used for obtaining gene codon and residue positions.  
* length - The number of nucleotides to take from the read after the starting nucleotide. Nucleotides after this length are removed.  

## Using the scripts to process sequencing reads into fitness scores
The following is a flowchart of the scripts' working order:  
  
![scripts workflow](https://github.com/johnchen93/DMS-FastQ-processing/blob/master/script_flowchart.png)

Running times for merging paired end reads is the slowest step, and could take a few minutes per million paired end reads. The run time will scale with the number of reads given to the program. Run times for all other scripts are shorter than a few minutes total.

The code is run in the following order:

### 0. FastQ file sorting  
- The "samples_to_csv.py" script is provided as an example on how to organize the FastQ files. This script will need to be modified if the filenames have different naming conventions.  
- Initially, the only folder will be 'input', which contains subfolders that contain individual sequencing files in fastQ format (fwd and rev separate). (see "setting up the environment" above)  
- Use the 'samples_to_csv.py' script to gather all filenames in the directories. The script will also parse the experimental conditions based on the file name. This will provide an output file called "parsed_samples.txt", in tab separated format.
- Copy and paste the output into the excel workbook named "sample_reference.xlsx", and name the sheet "samples". This will be used by later scripts to systematically process all files.  
- The "sample_reference.xlsx" will also have a sheet named "trim_info". As the name implies, it contains information on how reads from each sample should be trimmed, which is specific to our experimental and sequencing set up.  

### 1. Merge fastq fwd and reverse reads
- Run the 'merge_reads.py' script to merge forward and reverse reads. The results will be placed in the newly created "merged_reads" folder.  
- The script uses information in "sample_reference.xlsx" to process all files in the 'input' folder.  
- Results include:  
        a) a file containing the unique DNA sequences observed after quality filtering as well as their count.  
        b) a pdf file with plots showing the distribution of quality scores (including sequences that are filtered out)  
        c) a file containing summary information of the merger process  
        d) a single report file that details the quality parameters of all files processed, useful for finding samples for later scripts  

### 2a. Identify variants

- Run the 'identify_variants.py' script to find both codon and amino acid mutations from the merged fastq reads. Results are placed into the "processed_counts" folder.  
- This script also merges the separate 'groups' and 'replicates' from the experimental design, into a single file for each unique condition (drug, concentration, temperature)  
- Results include:  
        a) a file containing the single variants identified (or wt), the read count, and extra identifying information. Variants with more than one codon mutation are excluded.  
        b) a file that records the distribution of the number of mutated codons. wt has 0 mutated codons, single mutants have 1, etc. This is just for record keeping on how many double and triple mutants there are.  

### 2b. Calculate mutations expected from sequencing errors
    
- Sequencing errors is used to refer to all errors accumulated during the sequencing process, including PCR steps.  
- Run the 'count_wt_errors.py' script, followed by the 'expected_sequencing_error.py' script. Results are placed into the newly created 'wt_errors' folder.  
- The first script will use the fastq reads from sequencing wt dna to count the number of mutations (A>G, A>C, A>T, etc.).   
- The second script will calculate a) the sequencing error rate per position and b) the average chance of specific base changes during mutation (e.g. freq of A>G mutations out of all mutations from A). These numbers are then used to calculate the expected frequency of codon or amino acid mutations (including synonymous) due to single nucleotide substitutions on a wt codon.  

### 3. Filter non-selected library for variants due to noise
    
- Run the 'filter_nosel_library.py' script. Results are saved in the newly created 'scoring' folder.  
- The variants identified in the non-selceted library (2a) are filtered using the expected mutation rates due to sequencing errors (2b). By default, non-selected variants with a count less than 2x the count expected (wt count * expected mutation rate for the wt codon), or a count less or equal to 5 are filtered out.  
- The "filter_nosel_library.py" script has an optional argument 'libs' in the code that can be modified to select which experimental condition should be used for noise filtering.  
- The non-selected library with filtering information is saved in 'scoring', with one file for codon variants and another for amino acid variants.  

### 4. Calculate fitness score

- Run the 'calculate_fitness_score.py' script. Results are saved in the 'scoring' folder.  
- The script calculates fitness scores of a selected condition relative to the non-selected condition.  
- The script produces 3 outputs per input condition for each of codon or amino acid variants:  
        a) a file with the fitness scores of each variants separated by replicate.  
        b) a file with the fitness scores averaged between replicates  
        c) a file containing information on all steps used in scoring. This is just a record.  
