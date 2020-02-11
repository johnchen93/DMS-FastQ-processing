'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import pandas as pd
from libraries.calculate_expected_mutations import expected_noise

def noiseMask():
    '''
        This script calculates the expected mutations for each wt codon based on sequencing errors alone.
    '''
    path = 'wt_errors'
    # load the sequencing errors detected in sequencing wt dna
    m = pd.read_csv(f'{path}/wt_mismatches.txt', sep='\t')
    m = m[['pos','A','T','C','G','wt_base','total']] # only collect the necessary info
    m = m.groupby(['pos','wt_base']).sum()
    
    # calculate the positional error rate, i.e. mismatches per position divided by total reads at the position
    pos_err = m[['A','T','C','G']].sum(axis=1)/m['total']
    mean_pos_err = pos_err.mean()
    
    # calculate the mean proportion of each mismatch outcome, e.g. proportion of A>T transitions, A>C, A>G, etc.
    mut_freq = m[['A','T','C','G']].divide( m[['A','T','C','G']].sum(axis=1), axis=0 )
    mut_freq = mut_freq.groupby('wt_base').mean()[['A','T','C','G']]
    
    codon_noise_rate, aa_noise_rate = expected_noise(mean_pos_err, mut_freq)
    codon_noise_rate.to_csv(f'{path}/codon_noise_rate.txt', sep='\t')
    aa_noise_rate.to_csv(f'{path}/aa_noise_rate.txt', sep='\t')

if __name__=='__main__':
    noiseMask()