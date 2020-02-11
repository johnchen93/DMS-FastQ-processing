'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from libraries.helper import *
import pandas as pd
import os

def id_single_variants(df, wt_dna, dna_offset = 0):
    '''
    Given a dataframe containing a list of sequences as the index and their counts as values, compare the sequence to a
    wildtype reference and identify variants at the codon and amino acid levels. Has no knowledge of the dataset it was given
    and this should be taken care of outside the function.

    inputs:
        df - Dataframe of DNA sequences and counts, should match the length and intended region of the given wt reference
                assumes that each variant is already aggregated by count and no DNA level duplicates exist
        wt_dna - wildtype reference DNA sequence
        dna_offset - offset in nucleotide positions from the start of the gene, used for numbering purposes

    Output:
        Dictionary with a set of ouputs:
            dataframes with variant identities for both amino acid and codon variants
            dataframes with info on distribution of reads by number of codons mutated, for aa and codons
    '''
    # the output df will have the following columns
    # count - read count of the given variant
    # codon_mut - identity of the mutant
    # codon_mut_pos - codon position of the mutation
    # identity - classify as wildtype(wt), synonymous (syn), or variant (var)
    single_codon_headers = ["count","mut","pos", "identity"]
    single_aa_headers = ["count","mut","pos", "identity"]

    wt_prot = translate(wt_dna)

    single_codon_variants = {} # dict to keep track of the single codon mutants
    single_aa_variants = {} # dict to keep track of the single mutants

    # summary information on counts by number of codons/aa mutated, not strictly required, but interesting to have
    codon_mut_counts = {}
    aa_mut_counts = {}

    for seq, data in df.iterrows():
        var_prot = translate(seq)
        aa_diffs = differences(wt_prot, var_prot)
        aa_mismatches = sum(aa_diffs)
        codon_diffs = codonDifferences(seq, wt_dna)
        codon_mismatches = sum(codon_diffs)

        count = data[0] # the only column is "count"

        # keep track of only the wt count (nucleotide level), and variants mutated at a single codon (synonymous included)
        if codon_mismatches == 0: # has to be wildtype
            single_codon_variants['_wt'] = [count, "WT", None, "WT"]
            single_aa_variants['_wt'] = [count, "WT", None, "WT"]

        elif codon_mismatches == 1 : # single codon, keep
            for i in range(len(codon_diffs)):
                if codon_diffs[i]: # true is mismatched
                    codon_pos = i + dna_offset//3 + 1
                    codon = seq[i*3:i*3+3]
                    c_index = "(c.{pos}.{pre}>{post})".format(
                            pos = codon_pos,
                            pre=wt_dna[i*3:i*3+3],
                            post=codon)

                    aa_index, codon_id = "", ""

                    if aa_mismatches == 0: #codon different but amino acid same, synonymous variant
                        codon_id = "syn"

                        aa_index = "(p.{pos}.=)".format(
                            pos = codon_pos
                            )

                    elif aa_mismatches == 1:
                        codon_id = "var"

                        aa_index = "(p.{pre}.{pos}.{post})".format(
                            pos = codon_pos,
                            pre = AA_CODES[wt_prot[i]],
                            post= AA_CODES[var_prot[i]])

                    # organize and store codon level information
                    codon_data = [count, codon, codon_pos, codon_id]
                    single_codon_variants[c_index] = codon_data

                    # organize the amino acid information
                    aa_data = [count, var_prot[i], codon_pos, codon_id]
                    # need to add to existing mutation if already there, as diff codons can represent the same amino acid
                    try:
                        single_aa_variants[aa_index][0] += count
                    except:
                        single_aa_variants[aa_index] = aa_data

                    break # stop checking if one codon analyzed, only one mismatch expected
        # note distribution of all mutants in the library
        try:
            codon_mut_counts[codon_mismatches] += count
        except:
            codon_mut_counts[codon_mismatches] = count
        try:
            aa_mut_counts[aa_mismatches] += count
        except:
            aa_mut_counts[aa_mismatches] = count

    outdict = {"single_codon_variants":pd.DataFrame.from_dict(single_codon_variants, columns=single_codon_headers, orient = "index"),
                "single_aa_variants":pd.DataFrame.from_dict(single_aa_variants, columns=single_aa_headers, orient = "index"),
                "codon_mut_counts":pd.DataFrame.from_dict(codon_mut_counts, columns=['count'], orient = "index"),
                "aa_mut_counts":pd.DataFrame.from_dict(aa_mut_counts, columns=['count'], orient = "index")  }

    return outdict