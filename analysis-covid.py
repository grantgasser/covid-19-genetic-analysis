import os
# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC

# the replication machinery in bacteria 
# can only copy a DNA sequence at a maximum 
# speed of about 1000 nucleotides per second
NUCS_COPIED_PER_SEC = 1000
PROMOTER_SEQ = 'tataaaa'

# file paths
DATA_PATH = 'data/'

"""
Use biopython moving forward
"""

def read_genome(path):
    """Reads genome from text file into a string."""
    with open(path, 'r') as f:
        genome = f.read().replace('\n', '').lower()

    return genome

def count_nucleotide(seq):
    """Count the nucleotides of each base (A, C, G, U or T?) of the RNA"""
    nucleotide_count = {}

    for nucleotide in seq:
        nucleotide_count[nucleotide] = nucleotide_count.get(nucleotide, 0) + 1
    
    return nucleotide_count


def complement(seq):
    """Returns the complement of the sequence"""
    complements = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_seq = ''

    for nucleotide in seq:
        complement_seq += complements[nucleotide]

    return complement_seq

def reverse_seq(seq):
    return seq[::-1]

def reverse_complement(seq):
    seq = complement(seq)
    return reverse_seq(seq)

def main():
    covid_genome = read_genome('covid-19-genome')

    # 30318, 1 byte per character, so 30KB
    print('\nLength of gene sequence:', len(covid_genome))

    # count nucleotides (isn't covid-19 a single-strand RNA virus?)
    # so RNA should be A, C, G, and U instead of T?
    print('\nCounts:', count_nucleotide(covid_genome))
    print('Note that A ~= T and C ~= G which is a bit unexpected because its single stranded right?')
    print("'AT' rich... (A + T) > (G + C)")

    # complement
    comp = complement(covid_genome)
    print('\nComplement:', comp[:10], '...')

    # reverse complement (switch to 5'--------3')
    reverse_comp = reverse_complement(covid_genome)
    print('\nReverse Complement:', reverse_comp[-11:-1], '...')

    # test 1
    assert(reverse_seq(reverse_comp) == comp)

    # test 2
    test = 'atgtagttagctca'
    answer = 'tgagctaactacat'
    assert(reverse_complement(test) == answer)

    # QUESTION: which is the coding strand and which is the template strand?
    if PROMOTER_SEQ in reverse_comp:
        print('\nReverse complement is the coding strand')
        print('The genome is the template strand.')
    else:
        print('\nNot sure.')

    

    # read homo sapiens fragment
    sapiens_genome = read_genome('covid-19-genome')
    print('\nHomo Sapies count:', count_nucleotide(sapiens_genome))

main()