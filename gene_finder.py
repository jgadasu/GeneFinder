
# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: John-Edwin Gadasu
This carefully created Python Code has many functions that can be used for the manipulation of strings of nucleotides'
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
import doctest

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide=='':
        return ''
    val=''
    for i in range(len(nucleotide)):
        b=nucleotide[i]
        if b=='A':
            val+='T'
        elif b=='G':
            val+='C'
        elif b=='T':
            val+='A'
        elif b=='C':
            val+='G'
    return val
#doctest.run_docstring_examples(get_complement,globals(), verbose=True)




def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement('')
    ''
    """
    retval=get_complement(dna)
    return retval[::-1]
#doctest.run_docstring_examples(get_reverse_complement,globals(), verbose=True)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop1='TAG'
    stop2='TAA'
    stop3='TGA'             #The stop codons
    import re
    dna_array=re.findall('...',dna)         #Splits string into a list of triple-stringed codons

    if stop1 in dna_array:
        return dna[0:(dna_array.index(stop1))*3]
    elif stop2 in dna_array:
        return dna[0:(dna_array.index(stop2))*3]
    elif stop3 in dna_array:
        return dna[0:(dna_array.index(stop3))*3]
    else:
        return dna
#doctest.run_docstring_examples(rest_of_ORF,globals(), verbose=True)


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    >>> find_all_ORFs_oneframe('ATGGTATATATACACACCTAGATAGTGTCG')
    ['ATGGTATATATACACACC']


    """
    orfList = []
    i = 0
    while i < len(dna):
        #frame = ""
        tempCodon = dna[i:(i+3)]
        if tempCodon == 'ATG':
            frame = rest_of_ORF(dna[i:])
            orfList.append(frame)
            skip = len(frame)
            i += len(frame)
        else:
            i+=3
    return orfList
#doctest.run_docstring_examples(find_all_ORFs_oneframe,globals(), verbose=True)



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    dna=dna[dna.find('ATG'):]
    lists=[]
    c=(get_reverse_complement(dna))
    c=c[c.find('ATG'):]
    a=find_all_ORFs_oneframe(dna)
    b=(find_all_ORFs_oneframe(c))
    lists=a+b
    return lists
#doctest.run_docstring_examples(find_all_ORFs_both_strands,globals(), verbose=True)


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    list=[]
    var=find_all_ORFs_both_strands(dna)
    for i in var:
        list.append(len(i))
    return var[list.index(max(list))
]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    lists=[]
    length_list=[]
    i=0
    while i <=num_trials:
        a=shuffle_string(dna)
        lists.append(a)
        i+=1
    for i in lists:
        a=longest_ORF(i)
        length_list.append(longest_ORF(a))
    newstr=''
    for i in length_list:
        if len(i)>len(newstr):
            newstr=i
    return newstr



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    lists=[]
    strand=[]
    from amino_acids import aa_table
    import re
    lists = re.findall('...',dna)  #turns a string into an array of codons
    for elem in lists:
        strand.append(aa_table[elem])
    return ''.join(strand)
#doctest.run_docstring_examples(coding_strand_to_AA,globals(), verbose=True)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold=len(longest_ORF(dna))
    a=find_all_ORFs_both_strands(dna)
    lists=[]
    for i in a:
        if len(i)>=threshold:
            lists.append(coding_strand_to_AA(i))
    return lists


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    from load import load_seq
    dna=load_seq("./data/X73525.fa")
    print (gene_finder(dna))
