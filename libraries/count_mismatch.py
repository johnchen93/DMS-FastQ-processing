'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import pandas as pd

def count_mismatch(df, wt_dna, offset=0, mismatch_limit=10, silent=True):
    '''Go through DNA sequences contained in 'df' and compare to the expected WT
        Inputs:
            df - a dataframe with DNA sequences as the index and a column with the count of each sequence. Sequences should have the same length as the wt sequence.
            wt_dna - a string, used as the wildtype dna
            offset - starting position of the wt_dna input, used to calculate the actual position of each nucleotide in the gene
            mismatch_limit - int, sequences with mismatches to wt_dna greater than this limit are excluded from the analysis
            
        Outputs:
            mismatch_df - a dataframe with counts of the number of mismatches to wt_dna at each position separated by the final nucleotide
            mut_count_df - a dataframe with counts of the number of reads with a given number of mismatches to wt_dna
            dict - a dictionary with summary information, includes total number of mismatches, total number of bps considered, total number of reads processed, and the number of reads exceeding the mismatch limit
    '''
    
    tot_mismatch = 0
    tot_bp = 0
    tot_reads = 0
    
    mut_counts = {} # keep track of the number of sequences with a given number of mutations from wt
    mismatched_base = {} # dict for position and type of mismatch
    for j in range(len(wt_dna)):
        pos = j + 1 +offset
        mismatched_base[pos] = {"A":0,"T":0,"C":0,"G":0,"wt_base":wt_dna[j]} # create the dictionary
    
    # extra variable to keep track of frame shifted sequences
    reads_over_mismatch = 0
    
    # loop through every sequence
    for seq, row in df.iterrows():
        length = len(seq)
        count = row["count"]
        diffs = differences(wt_dna,seq)
        mismatches = (sum(diffs))
        # debug output to inform sequences suspected of having frame shifts, then skip to the next sequence
        if mismatches > mismatch_limit or length!=len(wt_dna):
            if not silent:
                print(f"suspected frameshifted read: {mismatches} mismatches, observed {count} times")
                print(wt_dna)
                print( "".join([ "*" if x else " " for x in diffs]) )
                print(seq)
            reads_over_mismatch += count
            continue
            
        tot_mismatch+=mismatches*count
        tot_bp+=length*count
        tot_reads+=count
        try:
            mut_counts[mismatches] += count # add to a list to plot
        except:
            mut_counts[mismatches] = count

        # C: for each mismatch, get the erroneous base and the position where it occurs
        # quantify as percent of all reads. Also cumulative, the same sequence can contribute to multiple mismatch positions and types
        if mismatches > 0:
            for j in range(length):
                pos = j + 1 +offset # off set to get real dna position
                if diffs[j]: # different
                    mismatched_base[pos][seq[j].upper()] += count # index by position, then by base

    mismatch_df = pd.DataFrame.from_dict(mismatched_base, orient = "index")
    mut_count_df = pd.DataFrame.from_dict(mut_counts, orient='index')
    
    return mismatch_df, mut_count_df, {"bp_errors":tot_mismatch,"bp_total":tot_bp, "total_reads":tot_reads, "frame_shifted_sequences":reads_over_mismatch}
    
def differences(s1, s2):
    """Calculate the Hamming distance between two bit strings
    returns the full list of booleans, where true indicates a mismatch"""
    #assert len(s1) == len(s2) # assume, since assert and try except don't work with jit
    return [c1 != c2 for c1, c2 in zip(s1, s2)]