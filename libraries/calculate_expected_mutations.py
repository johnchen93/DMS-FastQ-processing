'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import pandas as pd
from libraries.helper import *

# codon table
CODON_TABLE = {
        'TTT':'F', 'TCT':'S', 'TAT':'Y', 'TGT':'C',
        'TTC':'F', 'TCC':'S', 'TAC':'Y', 'TGC':'C',
        'TTA':'L', 'TCA':'S', 'TAA':'*', 'TGA':'*',
        'TTG':'L', 'TCG':'S', 'TAG':'*', 'TGG':'W',
        'CTT':'L', 'CCT':'P', 'CAT':'H', 'CGT':'R',
        'CTC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
        'CTA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
        'CTG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
        'ATT':'I', 'ACT':'T', 'AAT':'N', 'AGT':'S',
        'ATC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
        'ATA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
        'ATG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
        'GTT':'V', 'GCT':'A', 'GAT':'D', 'GGT':'G',
        'GTC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
        'GTA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
        'GTG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G'
    }
    
def expected_noise(error_rate, avg_mismatch_rate):
    # calculate the chance of a wt codon changing to the other 9 codons via single amino acid substitutions
    # standard deviation used to be considered, but since it wasn't ever used, it has been omitted
    
    '''
        Inputs:
        error_rate - float. The average positional error rate from sequencing.
        avg_mismatch_rate - dataframe. A dataframe with the mean proportion of specific substitutions (e.g. A>T, A>C, A>G). The index is the starting base, and the column is the resulting base.
        
        Outputs:
        codon_mut_rates_df - data frame with probabilities of transition between codons through single nt sequencing errors
        aa_mut_rates_df - data frames with probabilities of transition between amino acids through single nt sequencing errors
    '''
    sub_dict = {} # stores proprotion of substitutions resulting in specific codon changes
    codon_list = sorted(CODON_TABLE.keys())

    for wt in codon_list:
        sub_dict[wt] = {}
        for mut in codon_list:
            mismatches = differences(wt, mut)
            if sum(mismatches) == 1:
                # error rate only really applies to single mutations
                # higher order mutations are negligible and we would use the higher value
                # to be conservative anyway
                for i in range(len(mismatches)):
                    different = mismatches[i]
                    if different:
                        mismatch_frac = avg_mismatch_rate.loc[wt[i], mut[i]]

                        sub_dict[wt][mut] = mismatch_frac * error_rate
            else:
                sub_dict[wt][mut] = 0

    # turn the collected data to dataframes
    codon_mut_rates = pd.DataFrame.from_dict(sub_dict, orient = "index")
    # codon_mut_rates.to_csv(os.path.join(outpath,f"{name}_expected_codon_error.txt"), sep='\t')

    # collect codon to amino acid error rates, including synonymous mutations
    aa_mut_rates = {}
    for wt in codon_list:
        aa_mut_rates[wt] = {}
        wt_aa = CODON_TABLE[wt]
        for mut in codon_list:
            mut_aa = CODON_TABLE[mut]
            #if mut_aa != wt_aa: # if you wish to exclude synonymous mutations
            if True: # place holder, by pass the logic
                try: # sum the averaged and the variances
                    aa_mut_rates[wt][mut_aa] += codon_mut_rates.loc[wt, mut]
                except:
                    aa_mut_rates[wt][mut_aa] = codon_mut_rates.loc[wt, mut]
            else:
                aa_mut_rates[wt][mut_aa] = 0 #exclude wt to wt changes

    # convert the wt codon to amino acid error rates to dataframes
    aa_mut_rates_df = dict_to_df(aa_mut_rates, orient = "index")
    # aa_mut_rates_df.to_csv(os.path.join(outpath,f"{name}_expected_aa_error.txt"), sep='\t')
    
    return codon_mut_rates, aa_mut_rates_df