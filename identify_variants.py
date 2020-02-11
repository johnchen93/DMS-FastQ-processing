'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from libraries.id_variants import id_single_variants
from libraries.helper import *
import pandas as pd

def id_from_variants( ):
    vim_wt_dna = "".join(['ATGGGATTCAAACTTTTGAGTAAGTTATTGGTCTATTTGACCGCGTCTATCATGGCTATTGCGAGCCCGCTCGCTTTTTCCGTAGATTCTAGCGGAGAATATCCGACAGTCAGCGAA',
                'ATTCCGGTCGGGGAGGTCCGGCTTTACCAGATTGCCGATGGTGTTTGGTCGCATATCGCAACGCAGTCGTTTGATGGCGCAGTCTACCCGTCCAATGGTCTCATTGTCCGTGATGGT',
                'GATGAGTTGCTTTTGATTGATACAGCGTGGGGTGCGAAAAACACAGCGGCACTTCTCGCGGAGATTGAGAAGCAAATTGGACTTCCTGTAACGCGTGCAGTCTCCACGCACTTTCAT',
                'GACGACCGCGTCGGCGGCGTTGATGTCCTTCGGGCGGCTGGGGTGGCAACGTACGCATCACCGTCGACACGCCGGCTAGCCGAGGTAGAGGGGAACGAGATTCCCACGCACTCTCTT',
                'GAAGGACTTTCATCGAGCGGGGACGCAGTGCGCTTCGGTCCAGTAGAACTCTTCTATCCTGGTGCTGCGCATTCGACCGACAACTTAATTGTGTACGTCCCGTCTGCGAGTGTGCTC',
                'TATGGTGGTTGTGCGATTTATGAGTTGTCACGCACGTCTGCGGGGAACGTGGCCGATGCCGATCTGGCTGAATGGCCCACCTCCATTGAGCGGATTCAACAACACTACCCGGAAGCA',
                'CAGTTCGTCATTCCGGGGCACGGCCTGCCGGGCGGTCTTGACTTGCTCAAGCACACAACGAATGTTGTAAAAGCGCACACAAATCGCTCAGTCGTTGAGTAA'])
                
    # paths for input and output 
    path = 'merged_reads'
    outpath = 'processed_counts'
    
    make_sure_path_exists(outpath)
    
    # load list of files
    files = pd.read_csv('merged_reads/fastq_merge_report.txt', sep='\t')
    files = files[~(files['name'].str.startswith('wt'))]['name'] # get all non-wt files
    
    offset_info = pd.read_excel('sample_reference.xlsx', sheet_name='trim_info').set_index('group')
    
    out = {}
    for f in files:
        sep = f.split('_')
        label = '_'.join(sep[:2])
        group = sep[2][-1]
        rep = sep[3][-1]
        offset = offset_info.at[int(group),'offset'] # dna offset
        length = offset_info.at[int(group),'length'] # dna length after offset

        df = pd.read_csv(f"{path}/{f}/{f}.txt", sep='\t', index_col=0)
        wt = vim_wt_dna[offset:offset+length]
        
        results = id_single_variants(df, wt, offset)

        # extract results and add extra info
        single_codon_variants = results['single_codon_variants']
        single_aa_variants = results['single_aa_variants']
        single_codon_variants['group'] = group
        single_aa_variants['group'] = group
        single_codon_variants['rep'] = rep
        single_aa_variants['rep'] = rep

        codon_mut_counts = results['codon_mut_counts']
        codon_mut_counts['group'] = group
        codon_mut_counts['rep'] = rep
        aa_mut_counts = results['aa_mut_counts']
        aa_mut_counts['group'] = group
        aa_mut_counts['rep'] = rep
        
        if out.get(label):
            out[label]['codon_var'].append(single_codon_variants)
            out[label]['aa_var'].append(single_aa_variants)
            out[label]['codon_mut_counts'].append(codon_mut_counts)
            out[label]['aa_mut_counts'].append(aa_mut_counts)
        else:
            out[label] = {'codon_var':[single_codon_variants],'aa_var':[single_aa_variants],'codon_mut_counts':[codon_mut_counts],'aa_mut_counts':[aa_mut_counts]}

    # aggregate and save the data from each label (condition)
    for label, data in out.items():
        prepath = f"{outpath}/{label}"
        make_sure_path_exists(prepath)
        for name, data_group in data.items():
            d = pd.concat(data_group)
            d.to_csv(f'{prepath}/{label}_{name}.txt', sep='\t')
    
if __name__=='__main__':
    id_from_variants()