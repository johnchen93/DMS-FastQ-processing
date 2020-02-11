'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

# import modules
import gzip
import re
import math
import numpy as np
import pandas as pd
from array import array
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import errno

# translation table for reverse complementing from Enrich2
dna_trans = str.maketrans("actgACTG", "tgacTGAC")

# RE for just the tile info to match two reads
header_tile = re.compile("@(?P<MachineName>.+)"
                            ":(?P<Lane>\d+)"
                            ":(?P<Tile>\d+)"
                            ":(?P<X>\d+)"
                            ":(?P<Y>\d+)")

# fastq object class
# from https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py
class Fastq(object):
    """Fastq object with name and sequence
    """
    
    __slots__ = ('name', 'seq', 'name2', 'qual', 'qbase')

    # qBase of 33 is the default for Sanger Standard fastQ, also illumina 1.8
    # modified __init__ combining original and Enrich2
    def __init__(self, name, seq, name2, qual, convert_quality = True, qbase=33):
        self.name = name
        self.seq = seq
        self.name2 = name2
        # converts the quality to a list of signed chars "b", then iterates through it and gets the integer quality score
        # by subtracting the quality off set from it (33 by default for sanger standard)
        # to reverse the operation, simply use [x + self.qbase for x in self.quality]
        if convert_quality:
            self.qual = [x - qbase for x in array('b', qual.encode()).tolist()]
        else:
            self.qual = qual
        self.qbase = qbase


    def toString(self):
        return (self.name + "\n" + self.seq + "\n" + self.name2 + "\n" +
                     array('b', [int(x) + self.qbase for x in self.qual]).tostring() + "\n")

    def min_quality(self):
        """
        Return the minimum Phred-like quality score.
        """
        return min(self.qual)

    # probably won't use, and stick with the full header if needed
    def getShortname(self, separator):
        if separator:
            self.temp = self.name.split(separator)
            del(self.temp[-1])
            return separator.join(self.temp)
        else:
            return self.name
            
    # probably won't use, plan to process everything internally
    def write_to_file(self, handle):
        handle.write(self.name + "\n" + self.seq + "\n" + self.name2 + "\n" +
                     array('b', [int(x) + self.qbase for x in self.qual]).tostring() + "\n")

    # function to trim, cannot undo
    # the first base is 0
    def trim(self, start, length):
        start = start
        end = start + length
        self.seq = self.seq[start:end]
        self.qual = self.qual[start:end]

    def revcomp(self):
        """
        Reverse-complement the sequence in place. Also reverses the array of
        quality values.
        """
        self.seq = self.seq.translate(dna_trans)[::-1]
        self.qual = self.qual[::-1]

# modified open statement to work with .gz compression files
def myopen(infile, mode="rt"):
    if infile.endswith(".gz"):
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)

# returns a generator, can just hold on to the reference and ask for the next value
# using iter.next()
def fastq_parser(infile):

    buffer_size = 100000
    '''Takes a fastq file infile and returns a fastq object iterator
    '''
    eof = False
    leftover = ""
    with myopen(infile) as f:

        while not eof:
            buf = f.read(buffer_size)

            if len(buf) < buffer_size:
                eof = True

            buf = leftover + buf
            lines = buf.split('\n')
            fastq_count = len(lines) // 4

            if not eof: # handle lines from the trailing partial FASTQ record
                dangling = len(lines) % 4

                if dangling == 0: # quality line (probably) incomplete
                    dangling = 4
                    fastq_count = fastq_count - 1
                # join the leftover lines back into a string
                leftover = '\n'.join(lines[len(lines) - dangling:])

            for i in range(fastq_count):
                yield Fastq(*lines[i * 4:(i + 1) * 4])

def hammingDistance(s1, s2):
    """Calculate the Hamming distance between two bit strings
    returns just the hamming distance"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def differences(s1, s2):
    """Calculate the Hamming distance between two bit strings
    returns the full list of booleans, where true indicates a mismatch"""
    #assert len(s1) == len(s2) # assume, since assert and try except don't work with jit
    return [c1 != c2 for c1, c2 in zip(s1, s2)]

def resolveMismatches(s1, s2, max_mismatches = 20, minQ = 10, expectedError = 1):
    '''
    This function takes in two sequence objects, compares them for mismatches, and tries
    to resolve the mismatch. Where the forward and reverse reads do not match, the one with the higher score
    if selected. If both scores are the same, the base  will be taken from the read where the current position
    is closer to the start of the read.
    Posterior quality scores are calculated based on the Q scores of the forward and reverse reads using a precalculated
    lookup table. The calculation is slightly different depending on whether the cases are matching or not.
    The posterior quality equations are from Edgar and Flycbjerg 2015 "Error Filtering, pair assembly and error
    correction for next-generation sequencing reads"

    Inputs
    *s1 - sequences object with fields name, seq, name2 and qual
    *s2 - sequences object with fields name, seq, name2 and qual
    *max_mismatches - maximum number of mismatches tolerated between s1 and s2,
        the sequence will be ignored if this limit is surpassed
    *minQ - minimum per base Q score, a fastQ object will not be returned if any bases are below this limit
    *expectedError - maximum limit on the sum of all expected errors from Q scores in the sequence,
        a fastQ object will not be returned if the limit is surpassed

    Output
    returns a dictionary with the following key:value pairs
    "status": notes whether the sequence was resolved successfully or whether it failed to pass filters
    "fastq": returns a fastq object if the read was resolved and passes quality filters
    "mismatches": returns the number of mismatches in the read regardless of filtering, for profiling purposes
    "minQ": returns the min quality of this sequence if the mismatches were not too high, regardless of subsequent filtering, for profiling
    "error": returns the expected errors summed from all Q scores in the resolved read, if mismatches are below limit, and regardless of subsequent filtering, for profiling
    '''
    resolved_sequence = list(s1.seq)
    length = len(s1.seq)
    posterior_Q = np.zeros(length,dtype = int)
    diffs = differences(s1.seq, s2.seq)
    mismatches = sum(diffs)

    if mismatches <= max_mismatches:
        q1 = s1.qual
        q2 = s2.qual
        if mismatches == 0:
            for i in range(len(q1)):
                posterior_Q[i] = postQMatch[q1[i]][q2[i]]
            # return Fastq(s1.name, s1.seq, "+", matchingPosteriorQ(s1.qual, s2.qual), False), mismatches
        else:
            for i in range(length):
                isDifferent = diffs[i]

                if isDifferent: # mismatch
                    newBase = "N"
                    if q1[i]>q2[i]:
                        newBase = s1.seq[i]
                    elif q2[i]>q1[i]:
                        newBase = s2.seq[i]
                    else: # both quality scores are equal
                        if i < length/2: # closer to start
                            newBase = s1.seq[i]
                        else: # closer to end
                            newBase = s2.seq[i]
                    resolved_sequence[i] = newBase
                    posterior_Q[i] = postQMismatch[q1[i]][q2[i]]
                else:
                    posterior_Q[i] = postQMatch[q1[i]][q2[i]]
        # calculate regardless for profiling purposes, chances are afterwards each file will be run once anyway, so its not that big of a time loss
        minQ_out = min(posterior_Q)
        err_out = expectedErrorFromQScores(posterior_Q)
        if minQ_out >= minQ:
            if err_out < expectedError:
                return {"status":"passed", "fastq":Fastq(s1.name, "".join(resolved_sequence), "+", posterior_Q, False),
                "mismatches":mismatches, "minQ":minQ_out, "error":err_out}
            else: return  {"status":"high expected error", "fastq":None,
                "mismatches":mismatches, "minQ":minQ_out, "error":err_out}

        else: return {"status":"low base quality", "fastq":None,
                "mismatches":mismatches, "minQ":minQ_out, "error":err_out}

            # return Fastq(s1.name, "".join(resolved_sequence), "+", posterior_Q, False), mismatches
    else: return {"status":"mismatch limit", "fastq":None,
                "mismatches":mismatches, "minQ": -1, "error": -1}

def errorFromQ(q):
    return 10**(-q/10)

def expectedErrorFromQScores(qual):
    err = 0
    for i in qual:
        err += postErrors[i]
    return err

def posteriorQMatching(q1, q2, asQ = False):
    p1 = errorFromQ(q1)
    p2 = errorFromQ(q2)

    postError = (p1*p2/3)/(1-p1-p2 + 4*p1*p2/3)

    if asQ:
        return -10*math.log10(postError)
    else:
        return postError

def posteriorQMismatch(q1, q2, asQ = False):
    p1 = errorFromQ(q1)
    p2 = errorFromQ(q2)

    if p2<p1:
        p1, p2 = p2, p1

    postError = p1*(1-p2/3)/(p1+p2-4*p1*p2/3)

    if asQ:
        return -10*math.log10(postError)
    else:
        return postError

# make posterior Q lookup array to save some recalculating, probably not that big of a deal
# order of indexing doesn't matter, was made to be symmetric
# use numpy clip function to clamp between 0 and 40, other quality scores are not compatible
postQMatch = [[  int(round(np.clip((posteriorQMatching(x,y,True)),0,40) )) for y in range(0,41)  ] for x in range(0,41)]
postQMismatch = [[  int(round(np.clip((posteriorQMismatch(x,y,True)),0,40) )) for y in range(0,41)  ] for x in range(0,41)]
postErrors = [ 10**(-x/10)  for x in range (0,41) ]

# actual matching process
# assumes both fastq files have the same number of reads and that the reads are in the same order

# vim2 group2 expected WT just for benchmarking vim wt libraries
_WT = 'ATTCCGGTCGGGGAGGTCCGGCTTTACCAGATTGCCGATGGTGTTTGGTCGCATATCGCAACGCAGTCGTTTGATGGCGCAGTCTACCCGTCCAATGGTCTCATTGTCCGTGATGGT'

# need to set the proper trim start for each set
# the trim length is the same for NDM and VIM CCMs

## main function for calling to merge reads
def mergePairedFastQ(f_in, r_in, f_out_path, f_out_name ,f_trim_start=0, r_trim_start = 0, minQ = 10,
                     max_mismatches = 20 , max_expected_error = 1,trim_length=117):
    '''
    Takes a forward and reverse paired end sequencing fastq files and tries to merge them and get sequence counts.
    The workflow is as follows. The forward and reverse fastq files are read in parallel, and the headers are checked
    for a match. If the headers match, the two sequences are trimmed to remove flanking sequences and the reverse read
    is reverse complemented. The trimmed sequences are merged and quality filtered. If the mismatches exceeds the
    set limit, the sequence is currently discarded ( plan to align ). The merged fastq sequence is added to a dict of
    counts, if it passes quality filtering. The counts and the result of the merger process are outputted.

    Input:
    *f_in* the full path and filename of the forward fastq file including the .fastq and .gz extensions.

    *r_in* the full path and filename of the reverse fastq file including the .fastq and .gz extensions

    *f_out_name* the full path and filename of the output files excluding any extensions. Variantions will
    be added to this name for the different file outputs.

    *f_trim_start* all bases before this position in the forward fastq will be discarded. The first base is counted as 0.

    *r_trim_start* all bases before this position in the reverse fastqwill be discarded. The first base is counted as 0.

    *minQ* the minimum allowed posterior Q score at each basecall in the resolved sequence

    *max_mismatches* the max number of mismatches allowed between forward and reverse reads before trying to merge.
    Empirically determined to contain the majority of properly aligned sequences, with decent quality.

    *max_expected_error* the expected error of a resolved sequence is the sum of the error probabilities across
    all positions in the sequence, based on posterior Q scores of the merged sequence. If the expected error
    for a full sequence is greater than this limit it is discarded.

    *trim_length* starting from the f_trim_start and r_trim_start, the length of sequence to keep after wards.
    Bases after the end of the length is discarded.

    Output:
    A HDF5 store and matching tsv file for holding the identities of the resolved sequences and the number of
    each sequence observed.
    A pdf with graphs noting the distribution of the number of mismatches observed before tring to merge sequences,
    and the distribution of minimum posterior Q scores and expected errors.
    A text file with a report on the number of sequences encountered, and how many of those were merged successfully
    or failed to pass quality filters.

    '''
    make_sure_path_exists(f_out_path)
    f_out_filepath = os.path.join(f_out_path, f_out_name)

    header_mismatch = 0
    mismatch_over_limit = 0
    low_quality = 0
    high_expected_error = 0
    final_sequences = 0
    total_sequences = 0

    mismatch_list = []
    min_qual_list = []
    expected_err_list = []

    # creates a generator for each file
    f_seq = fastq_parser(f_in)
    r_seq = fastq_parser(r_in)

    counts = {}
    # loops through all fastq lines
    while(True):
        # get the fastQ objects
        try:
            curr_f = next(f_seq)
            curr_r = next(r_seq)
        # no sequences left
        except:
            break

        total_sequences += 1
        # check if they match in header abd sequence, should be a string to string comparison
        same_name = header_tile.search(curr_f.name).group() == header_tile.search(curr_r.name).group()

        resolved = None # reset the variable
        # the two sequences match exactly
        if same_name:
            # trim forward
            curr_f.trim(f_trim_start, trim_length)
            # trim (has to be done first to make sure trim is correct) and reverse complement
            curr_r.trim(r_trim_start, trim_length)
            curr_r.revcomp()

            # resolve the sequences, gets a dict of the analysis info and fastq if it passed cut off
            resolved = resolveMismatches(curr_f, curr_r, max_mismatches, minQ, max_expected_error)
            # reads are processed, quality filter here
            status = resolved["status"]
            if status == "passed":
                final_sequences += 1
                try: counts[resolved["fastq"].seq] += 1
                except: counts[resolved["fastq"].seq] = 1

            elif status == "low base quality": low_quality += 1
            elif status == "high expected error": high_expected_error += 1
            elif status == "mismatch limit": mismatch_over_limit += 1

            # note the number of mismatches, minimum quality, and expected error
            mismatch_list.append(resolved["mismatches"])
            min_qual_list.append(resolved["minQ"])
            expected_err_list.append(resolved["error"])
        else:
            header_mismatch += 1

    # report the results
    print(f"{f_out_name}: {final_sequences}/{total_sequences} passed")

    df = pd.DataFrame.from_dict(counts, columns=['count'], orient="index")
    df.to_csv(f_out_filepath + ".txt", sep='\t',na_rep="NA")

    # create a report file with the outcome of the merge
    with open(f_out_filepath + "_fqMerge_report.txt", 'w') as file:
        file.write("Quality Filters:\n\tMax Mismatches before merge <= "+ str(max_mismatches)
                    + "\n\tMin Posterior Q score of any position >= "+ str(minQ)
                    + "\n\tExpected Error per sequence < "+ str(max_expected_error) + "\n")
        file.write("header mismatch: "+ str(header_mismatch) + "\n")
        file.write("mismatches over limit: "+ str(mismatch_over_limit) + "\n")
        file.write("low quality: "+ str(low_quality) + "\n")
        file.write("high expected error: "+ str(high_expected_error) + "\n")
        file.write("sequences passed: "+ str(final_sequences) + "\n")
        file.write("sequences total: "+ str(total_sequences) + "\n")

    pp = PdfPages(f_out_filepath + "_fqmerge_graphs.pdf")
    make_hist(pp, mismatch_list, "Distribution of Mismatches between Fwd and Rev Reads",
    "mismatches", "counts", 100)
    make_hist(pp, min_qual_list, "Minimum posterior Q score (-1 means mismatches > " + str(max_mismatches) + ")",
    "posterior Q score", "counts", 40)
    make_hist(pp, expected_err_list, "Expected error per read (-1 means mismatches > " + str(max_mismatches) + ")",
    "Expected error per read", "counts", 40)
    pp.close()
    #return counts, mismatch_list, min_qual_list, expected_err_list

    # return the merging stats as a dictionary
    return {"name":f_out_name, "max mismatches allowed":max_mismatches, "min quality allowed":minQ, "max expected error allowed":max_expected_error,
            "header mismatch":header_mismatch, "mismatches over limit":mismatch_over_limit, "low quality":low_quality,"high expected error":high_expected_error,
            "sequences passed":final_sequences, "sequences total":total_sequences}

def make_hist( pdfPages , values, title = "histogram", x_label = "x axis", y_label = "y axis", bins = 40, weights = None):
    f, ax = plt.subplots()
    ax.hist(values, bins=bins, density=False, facecolor='green', alpha=0.75, weights = weights)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)

    f.savefig(pdfPages, format='pdf')
    plt.close(f)

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

