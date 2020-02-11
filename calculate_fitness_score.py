'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import pandas as pd
import numpy as np
from libraries.helper import *

def score():
    # center syn
    center_syn = True # option to subtract the mean of the synonymous scores from each group

    # pair each condition with the matching non-selected condition
    pairs = [
                ('amp128_37C','nosel_37C')
            ]
    
    outpath = 'scoring'
    merge_headers = ['index','mut','pos','identity','group','rep']
    
    for pair in pairs:
        for x in ['codon','aa']:
            sel = pd.read_csv(f'processed_counts/{pair[0]}/{pair[0]}_{x}_var.txt', sep='\t', index_col=0).reset_index()
            nosel = pd.read_csv(f'scoring/{pair[1]}_{x}_noise_filtered.txt', sep='\t')
            nosel = nosel[(nosel['above_noise']) | (nosel['mut']=='WT')]
            
            d = nosel.merge(sel, on=merge_headers, how='left', suffixes=['_nosel','_sel'])
            d['count_sel_filled']=d['count_sel'].fillna(1)
            
            col = [] # collection of data frames
            for name, data in d.groupby(['group','rep']):
                data.set_index('index', inplace=True)
                wt_df = data.loc['_wt',:]
                wt_ratio = wt_df['count_sel_filled']/wt_df['count_nosel']
                ''' Scoring formula:
                    fitness score = [ freq variant sel / freq variant nosel] / [freq wt sel / freq wt nosel]
                    freq = count variant / count total
                    
                    note: When the variant and wt freq are measured from the same sequencing sample, they have the same total.
                            In this case the 'count total' cancels out:
                                fitness score = [ count variant sel / count total sel ]  *  [ count wt nosel / count total nosel ]
                                                ----------------------------------------------------------------------------------
                                                [ count variant nosel / count total nosel ]  *  [ count wt sel / count total sel ]
                                                
                                fitness score = [ count variant sel/count variant nosel ] / [ count wt sel/count wt nosel ]
                '''
                data['score'] = data.apply(lambda x: (x['count_sel_filled']/x['count_nosel'])/wt_ratio , axis=1)
                data['log2_score'] = np.log2(data['score']) # transform fitness score to log2
                if center_syn:
                    data['log2_score'] = data['log2_score'] - data[data['identity']=='syn']['log2_score'].mean()
                col.append(data)
            
            out = pd.concat(col).reset_index()
            out.to_csv(f"{outpath}/{pair[0]}_{x}_scoring_record.txt",sep='\t',index=False) # save a version of the file with all calculation steps for convenience
            
            out_headers1 = ['index','mut','pos','identity','group','rep','wt_codon']
            out_headers2 = ['index','mut','pos','identity','group','wt_codon']
            # save a file with just the final log2 scores, but kept in separate replicates
            out[out_headers1+['log2_score']].to_csv(f"{outpath}/{pair[0]}_{x}_replicate_scores.txt",sep='\t',index=False)
            # save a file with the scores averaged
            out[out_headers2+['log2_score']].groupby(out_headers2).mean().reset_index().to_csv(f"{outpath}/{pair[0]}_{x}_averaged_scores.txt",sep='\t',index=False)
            
if __name__=='__main__':
    score()