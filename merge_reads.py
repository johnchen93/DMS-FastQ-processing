'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import os
import pandas as pd
import numpy as np
from libraries.helper import *
import libraries.fq_merge as fqm

def main():
    samples = pd.read_excel('sample_reference.xlsx', sheet_name='samples')
    trim = pd.read_excel('sample_reference.xlsx', sheet_name='trim_info').set_index('group')
    
    # d = samples.merge(trim, on='group',how='left')
    d=samples
    
    # filter settings
    minQ = 10
    max_mismatches = 20
    max_expected_error = 1
    
    outpath='merged_reads'
    out = {} # for saving filtering stats
    for name, data in d.groupby(['label','conc','temp','group','rep']):
        f1 = data[data['read']==1]['path'].values[0] # fwd read
        f2 = data[data['read']==2]['path'].values[0] # rev read
        label=name[0]
        conc=name[1]
        temp=name[2]
        group=name[3]
        rep=name[4]
        
        tr = trim.loc[group,:]
        ftrim=tr['fwd trim']
        rtrim=tr['rev trim']
        length=tr['length']
        
        if label in ['amp','ctx','mem']:
            id = f'{label}{conc}_{temp}C_G{group}_rep{rep}'
        elif label=='nosel':
            id = f'{label}_{temp}C_G{group}_rep{rep}'
        elif label=='wt':
            id = f'{label}_G{group}_rep{rep}'
    
        out[id] = fqm.mergePairedFastQ( f1, f2, f"{outpath}/{id}", f"{id}", f_trim_start=ftrim, r_trim_start=rtrim, minQ=minQ, max_mismatches = max_mismatches, max_expected_error=max_expected_error, trim_length=length )
    
    df = pd.DataFrame.from_dict(out, orient='index')
    df.to_csv(f'{outpath}/fastq_merge_report.txt', sep='\t', index=False)

if __name__ == "__main__":
    main()
