import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

"""
the replication machinery in bacteria can only copy a DNA sequence at a maximum 
speed of about 1000 nucleotides per second
"""
NUCS_COPIED_PER_SEC = 1000
PROMOTER_SEQ = 'TATAAAA'
CODON_LEN = 3
"""
codons => amino acids, basically tRNAs
NOTE: AUG: M is the start, UAG, UAA, and UGA stop
    - since translation is energy expensive, having 3 stops means more likely to stop 
    if RNA is scrambled by an insertion
""" 
CODON_TO_AA = {'AUG': 'M', 'GCG': 'A', 'UCA': 'S', 'GAA': 'E', 'GGG': 'G', 'GGU': 'G', 'AAA': 'K', 
    'GAG': 'E', 'AAU': 'N', 'CUA': 'L', 'CAU': 'H', 'UCG': 'S', 'UAG': 'STOP', 'GUG': 'V', 'UAU': 'Y', 
    'CCU': 'P', 'ACU': 'T', 'UCC': 's', 'CAG': 'Q', 'CCA': 'P', 'UAA': 'STOP', 'AGA': 'R', 'ACG': 'T', 
    'CAA': 'Q', 'UGU': 'C', 'GCU': 'A', 'UUC': 'F', 'AGU': 'S', 'AUA': 'I', 'UUA': 'L', 'CCG': 'P', 
    'AUC': 'I', 'UUU': 'F', 'CGU': 'R', 'UGA': 'STOP', 'GUA': 'V', 'UCU': 'S', 'CAC': 'H', 'GUU': 'V', 
    'GAU': 'D', 'CGA': 'R', 'GGA': 'G', 'GUC': 'V', 'GGC': 'G', 'UGC': 'C', 'CUG': 'L', 'CUC': 'L', 
    'CGC': 'R', 'CGG': 'R', 'AAC': 'N', 'GCC': 'A', 'AUU': 'I', 'AGG': 'R', 'GAC': 'D', 'ACC': 'T', 
    'AGC': 'S', 'UAC': 'Y', 'ACA': 'T', 'AAG': 'K', 'GCA': 'A', 'UUG': 'L', 'CCC': 'P', 'CUU': 'L', 'UGG': 'W'
}

# file paths
DATA_PATH = 'data/'

def read_genome(path):
    """Reads genome from text file into a string."""
    with open(path, 'r') as f:
        genome = f.read().replace('\n', '').upper()
    return genome

def count_nucleotide(seq):
    """Count the nucleotides of each base (A, C, G, U or T?) of the RNA"""
    nucleotide_count = {}
    for nucleotide in seq:
        nucleotide_count[nucleotide] = nucleotide_count.get(nucleotide, 0) + 1
    return nucleotide_count


def complement(seq):
    """Returns the complement of the sequence"""
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = ''
    for nucleotide in seq:
        complement_seq += complements[nucleotide]
    return complement_seq

def reverse_seq(seq):
    return seq[::-1]

def reverse_complement(seq):
    seq = complement(seq)
    return reverse_seq(seq)

def transcribe(seq):
    """Transcribes original coding to RNA (T => U) """
    return seq.replace('T', 'U')

def translate(rna_seq):
    """Translates RNA sequence into amino acids (protein sequence)"""
    protein = []
    i = 0
    while i < len(rna_seq):
        codon = rna_seq[i:i+CODON_LEN]
        aa = CODON_TO_AA[codon]
        if aa != 'STOP':
            rotein.append(CODON_TO_AA[codon])
            i += CODON_LEN   
        else:
            break 
    return protein

def main():
    covid_genome = read_genome('covid-19-genome')

    # 30318, 1 byte per character, so 30KB
    print('\nLength of gene sequence:', len(covid_genome))
    print('NOTE: not a multiple of 3, might be a problem.')

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
    test = 'ATGTAGTTAGCTCA'
    answer = 'TGAGCTAACTACAT'
    assert(reverse_complement(test) == answer)

    # QUESTION: which is the coding strand and which is the template strand?
    # if PROMOTER_SEQ in covid_genome:
    #     print('\ngenome is the coding strand')
    #     print('The reverse complement is the template strand.')
    # else:
    #     print('\nNot sure.')

    # read homo sapiens fragment
    sapiens_genome = read_genome('covid-19-genome')
    print('\nHomo Sapies count:', count_nucleotide(sapiens_genome))

    # BioPython (test if hand-written functions produce same result)
    seq = Seq(covid_genome, IUPAC.unambiguous_dna)
    assert(seq.reverse_complement() == reverse_comp)

    # transcription
    covid_transcribed = seq.transcribe()
    assert(covid_transcribed == transcribe(covid_genome))

    # translate to proteins (20 amino acids)
    # 4 letter nucleotide alphabet => 20 amino acid alphabet
    # 1:1 is impossible since 4 < 16
    # 2 at a time is not quite enough 4^2 = 16 < 20
    # need to do 3 at a time since 4^3 = 64 > 20, a codon

    # TODO: check if genome file is correct/complete (not a multiple of 3)
    # seq.translate() 


main()