'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import pandas as pd
from libraries.helper import *

def filterNaiveLibrary( libs = ['nosel_37C'] ): # libs is the only input, which is a list of directories in 'processed_counts' from which the variant counts should be used to determine if a variant in a library exceeds the sequencing noise threshhold. Replace or extend the list to add additional conditions

    vim_wt_dna = "".join(['ATGGGATTCAAACTTTTGAGTAAGTTATTGGTCTATTTGACCGCGTCTATCATGGCTATTGCGAGCCCGCTCGCTTTTTCCGTAGATTCTAGCGGAGAATATCCGACAGTCAGCGAA',
                'ATTCCGGTCGGGGAGGTCCGGCTTTACCAGATTGCCGATGGTGTTTGGTCGCATATCGCAACGCAGTCGTTTGATGGCGCAGTCTACCCGTCCAATGGTCTCATTGTCCGTGATGGT',
                'GATGAGTTGCTTTTGATTGATACAGCGTGGGGTGCGAAAAACACAGCGGCACTTCTCGCGGAGATTGAGAAGCAAATTGGACTTCCTGTAACGCGTGCAGTCTCCACGCACTTTCAT',
                'GACGACCGCGTCGGCGGCGTTGATGTCCTTCGGGCGGCTGGGGTGGCAACGTACGCATCACCGTCGACACGCCGGCTAGCCGAGGTAGAGGGGAACGAGATTCCCACGCACTCTCTT',
                'GAAGGACTTTCATCGAGCGGGGACGCAGTGCGCTTCGGTCCAGTAGAACTCTTCTATCCTGGTGCTGCGCATTCGACCGACAACTTAATTGTGTACGTCCCGTCTGCGAGTGTGCTC',
                'TATGGTGGTTGTGCGATTTATGAGTTGTCACGCACGTCTGCGGGGAACGTGGCCGATGCCGATCTGGCTGAATGGCCCACCTCCATTGAGCGGATTCAACAACACTACCCGGAAGCA',
                'CAGTTCGTCATTCCGGGGCACGGCCTGCCGGGCGGTCTTGACTTGCTCAAGCACACAACGAATGTTGTAAAAGCGCACACAAATCGCTCAGTCGTTGAGTAA'])
    wt = {} # make a list of wt codons, used to reference the error rate
    for i in range(0, len(vim_wt_dna), 3):
        pos = i//3 + 1
        wt[pos] = vim_wt_dna[i:i+3]
    wt = pd.DataFrame.from_dict(wt, orient='index', columns=['wt_codon']).reset_index().rename(columns={'index':'pos'})
    
    # load the expected mutation rates due to sequencing
    codon_err = pd.read_csv('wt_errors/codon_noise_rate.txt', sep='\t', index_col=0).reset_index().rename(columns={'index':'wt_codon'}).melt(id_vars='wt_codon',var_name='mut',value_name='rate')
    aa_err = pd.read_csv('wt_errors/aa_noise_rate.txt', sep='\t', index_col=0).reset_index().rename(columns={'index':'wt_codon'}).melt(id_vars='wt_codon',var_name='mut',value_name='rate')
    err = {'codon':codon_err, 'aa':aa_err}
    
    # load the nosel libraries
    path = 'processed_counts'
    noise_threshold = 2 # variants below this multiple of the expected noise will be marked as noise
    min_count = 5 # variants must have greater than the minimum number of counts
    outpath = 'scoring'
    make_sure_path_exists(outpath)
    for l in libs:
        for x in ['codon','aa']:
            lib = pd.read_csv(f"{path}/{l}/{l}_{x}_var.txt", sep='\t', index_col=0).reset_index()
            lib = lib.merge(wt, on='pos', how='left')
            lib = lib.merge(err[x], on=['wt_codon','mut'], how='left')
            
            data_list = [] # keep track of data from separate groups
            for name, data in lib.groupby(['group','rep']):
                wt_count = data[data['mut']=='WT']['count'].values[0]
                data['noise_count'] = data['rate'].multiply(wt_count)
                data['above_noise'] = data.apply(lambda x: x['count'] > x['noise_count']*noise_threshold and x['count'] > min_count, axis=1 )
                data_list.append(data)
            
            out = pd.concat(data_list).set_index('index')
            out.to_csv(f"{outpath}/{l}_{x}_noise_filtered.txt", sep='\t')
            
if __name__=='__main__':
    filterNaiveLibrary()