'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import pandas as pd
from libraries.count_mismatch import count_mismatch
from libraries.helper import *

def main():
    vim_wt_dna = "".join(['ATGGGATTCAAACTTTTGAGTAAGTTATTGGTCTATTTGACCGCGTCTATCATGGCTATTGCGAGCCCGCTCGCTTTTTCCGTAGATTCTAGCGGAGAATATCCGACAGTCAGCGAA',
                'ATTCCGGTCGGGGAGGTCCGGCTTTACCAGATTGCCGATGGTGTTTGGTCGCATATCGCAACGCAGTCGTTTGATGGCGCAGTCTACCCGTCCAATGGTCTCATTGTCCGTGATGGT',
                'GATGAGTTGCTTTTGATTGATACAGCGTGGGGTGCGAAAAACACAGCGGCACTTCTCGCGGAGATTGAGAAGCAAATTGGACTTCCTGTAACGCGTGCAGTCTCCACGCACTTTCAT',
                'GACGACCGCGTCGGCGGCGTTGATGTCCTTCGGGCGGCTGGGGTGGCAACGTACGCATCACCGTCGACACGCCGGCTAGCCGAGGTAGAGGGGAACGAGATTCCCACGCACTCTCTT',
                'GAAGGACTTTCATCGAGCGGGGACGCAGTGCGCTTCGGTCCAGTAGAACTCTTCTATCCTGGTGCTGCGCATTCGACCGACAACTTAATTGTGTACGTCCCGTCTGCGAGTGTGCTC',
                'TATGGTGGTTGTGCGATTTATGAGTTGTCACGCACGTCTGCGGGGAACGTGGCCGATGCCGATCTGGCTGAATGGCCCACCTCCATTGAGCGGATTCAACAACACTACCCGGAAGCA',
                'CAGTTCGTCATTCCGGGGCACGGCCTGCCGGGCGGTCTTGACTTGCTCAAGCACACAACGAATGTTGTAAAAGCGCACACAAATCGCTCAGTCGTTGAGTAA'])
                
    files = pd.read_csv('merged_reads/fastq_merge_report.txt', sep='\t')
    files = files[files['name'].str.startswith('wt')]['name']
    
    offset_info = pd.read_excel('sample_reference.xlsx', sheet_name='trim_info').set_index('group')
    
    path = 'merged_reads'
    outpath = 'wt_errors'
    make_sure_path_exists(outpath)
    
    bygroup = {}
    mismatch_list = []
    mut_count_list = []
    for f in files:
        sep = f.split('_')
        group = sep[1][-1]
        rep = sep[2][-1]
        dir = '_'.join(sep[:3])
        offset = offset_info.at[int(group),'offset']
        length = offset_info.at[int(group),'length']
        
        df = pd.read_csv(f"{path}/{dir}/{f}.txt", sep='\t', index_col=0)
        wt = vim_wt_dna[offset:offset+length]
        # count the errors
        mismatches, mut_counts, group_info = count_mismatch( df, wt, offset=offset )
        mismatches['total'] = group_info['total_reads']
        mismatches['group'] = group
        mismatches['rep'] = rep
        mut_counts['group'] = group
        mut_counts['rep'] = rep
        
        bygroup[f"G{group}_rep{rep}"] = group_info
        mismatch_list.append(mismatches)
        mut_count_list.append(mut_counts)
    
    # combine mismatches from all groups and replicates into one table
    mdf = pd.concat(mismatch_list).reset_index().rename(columns={'index':'pos'})
    cdf = pd.concat(mut_count_list).reset_index().rename(columns={'index':'mutations','0':'count'})
    
    mdf.to_csv(f"{outpath}/wt_mismatches.txt", sep='\t', index=False)
    cdf.to_csv(f'{outpath}/wt_mut_counts.txt', sep='\t', index=False)
    
if __name__=='__main__':
    main()