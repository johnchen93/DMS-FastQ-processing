'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import errno
import os
import string
import pandas as pd
import gzip

# single letter list of amino acids, sorted by type
aa_list = ['H', 'K', 'R',                # (+)
           'D', 'E',                     # (-)
           'C', 'M', 'N', 'Q', 'S', 'T', # Polar-neutral
           'A', 'G', 'I', 'L', 'P', 'V', # Non-polar
           'F', 'W', 'Y',                # Aromatic
           '*']

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

#: Conversions between single- and three-letter amino acid codes
AA_CODES = {
        'Ala' : 'A', 'A' : 'Ala',
        'Arg' : 'R', 'R' : 'Arg',
        'Asn' : 'N', 'N' : 'Asn',
        'Asp' : 'D', 'D' : 'Asp',
        'Cys' : 'C', 'C' : 'Cys',
        'Glu' : 'E', 'E' : 'Glu',
        'Gln' : 'Q', 'Q' : 'Gln',
        'Gly' : 'G', 'G' : 'Gly',
        'His' : 'H', 'H' : 'His',
        'Ile' : 'I', 'I' : 'Ile',
        'Leu' : 'L', 'L' : 'Leu',
        'Lys' : 'K', 'K' : 'Lys',
        'Met' : 'M', 'M' : 'Met',
        'Phe' : 'F', 'F' : 'Phe',
        'Pro' : 'P', 'P' : 'Pro',
        'Ser' : 'S', 'S' : 'Ser',
        'Thr' : 'T', 'T' : 'Thr',
        'Trp' : 'W', 'W' : 'Trp',
        'Tyr' : 'Y', 'Y' : 'Tyr',
        'Val' : 'V', 'V' : 'Val',
        'Ter' : '*', '*' : 'Ter',
        '???' : '?', '?' : '???'
}

# a table for pdb secondary structure definitions
SS_TABLE = {
        'H' : 'alpha-helix',
        'B' : 'isolated beta-strand',
        'E' : 'bonded beta-strand',
        'G' : '3/10-helix',
        'I' : 'pi-helix',
        'T' : 'H-bond turn',
        'S' : 'bend',
        ' ' : 'not classified'
}

# fastq object class
# taken from https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py

# __init__ is kind of slow, could just use a dict and make the methods general
class Fastq(object):
    """Fastq object with name and sequence
    """

    # idea to use slots from Enrich2
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
            self.qual = [x - qbase for x in array('b', qual).tolist()]
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

    # original function to trim, cannot undo, inspired by Enrich2
    # the first base is 0
    def trim(self, start, length):
        start = start
        end = start + length
        self.seq = self.seq[start:end]
        self.qual = self.qual[start:end]

    # reverse complement function taken from Enrich2
    def revcomp(self):
        """
        Reverse-complement the sequence in place. Also reverses the array of
        quality values.
        """
        self.seq = self.seq.translate(dna_trans)[::-1]
        self.qual = self.qual[::-1]

    # minimum quality function taken from Enrich2

def myopen(infile, mode="r"):
    if infile.endswith(".gz"):
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)

# returns a generator, can just hold on to the reference and ask for the next value
# using iter.next()
def fasta_parser(infile, buffer_size = 100000, lines_per_record = 2, fastq = False):

    # assumes that the fasta files holds the sequence on a single line or fixed number of lines like alignments
    #buffer_size = 100000
    if fastq:
        lines_per_record = 4
    """Takes a fastq file infile and returns a fastq object iterator
    """
    eof = False
    leftover = ""
    with myopen(infile) as f:

        while not eof:

            buf = f.read(buffer_size)

            if len(buf) < buffer_size:
                eof = True

            buf = leftover + buf
            lines = buf.split('\n')
            record_count = len(lines) // lines_per_record

            if not eof: # handle lines from the trailing partial FASTQ record
                dangling = len(lines) % lines_per_record

                if dangling == 0: # quality line (probably) incomplete
                    dangling = lines_per_record
                    record_count = record_count - 1
                # join the leftover lines back into a string
                leftover = '\n'.join(lines[len(lines) - dangling:])

            for i in xrange(record_count):
                if fastq:
                    yield Fastq(*lines[i * lines_per_record:(i + 1) * lines_per_record])
                else:
                    yield lines[i * lines_per_record:(i + 1) * lines_per_record]

# general fasta parser that takes multiline fastas and returns the sequence as a value in a dictionary
# where the key is the header
def fasta_parser_simple(infile):
    sequences = {}
    with open(infile,'r') as f_in:
        currheader = ""
        for line in f_in:
            if line.startswith(">"):
                currheader = line[1:].strip()
                sequences[currheader] = []
            else:
                sequences[currheader].append(line.rstrip())

    for k,v in sequences.items():
        # make multiline sequences a single entry
        sequences[k] = "".join(v)

    return sequences

# translation table for reverse complementing from Enrich2
dna_trans = str.maketrans("actgACTG", "tgacTGAC")

def revcomp(seq):
        """
        Reverse-complement the sequence
        """
        return seq.translate(dna_trans)[::-1]


def rStr(num,x):
    return str(round(num,x))

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

def codonDifferences(s1, s2):
    """Calculate codon differences

    Input:
    s1 - a string of DNA, divisible by 3, in frame with intended reading frame
    s2 - similar to s1, must be same length
    """
    return [s1[i:i+3] != s2[i:i+3] for i in range(0, len(s1),3)]


def translate(dna):
    return "".join(CODON_TABLE[dna[i:i+3]] for i in range(0, len(dna),3))

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def get_files_in_dir(mypath):
    f = []
    for (dirpath, dirnames, filenames) in os.walk(mypath):
        f.extend(filenames)
        break
    return filenames

def make_sure_path_exists(path):
    # recursively creates a series of directories if it does not already exist
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def html_report(template_dir,template_name,temp_var,output):
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template( template_name )
    html_out = template.render(temp_var)

    with open(output, 'w') as f_out:
        f_out.write(html_out)

def dict_to_df( dict, columns = None, orient = 'index', sort_index = True):
    df = pd.DataFrame.from_dict(dict, orient=orient)
    if sort_index:
        df = df.sort_index()
    if columns:
        df.columns = columns
    return df


def propagateSD(avg1, sd1, avg2, sd2):
    prop_var = (sd1*sd2)**2 + (sd1*avg2)**2 + (sd2*avg1)**2
    prop_sd = (prop_var)**(0.5)
    return prop_sd

def expected_aa_variants(wt_dna):
    prot = translate(wt_dna)
    total = 0
    for aa in prot:
        if aa in ['M','W']:
            total += 20
        else:
            total += 21

    return total
