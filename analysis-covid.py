import os

"""
Use biopython moving forward
"""

def count_base(seq):
    """Count the bases (A, C, G, U or T?) of the RNA"""
    base_count = {}

    for base in seq:
        if base is not '\n':
            base_count[base] = base_count.get(base, 0) + 1
    
    return base_count


def complement(seq):
    """Returns the complement of the sequence"""
    complements = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_seq = ''

    for base in seq:
        complement_seq += complements[base]

    return complement_seq

def reverse_seq(seq):
    return seq[::-1]

def reverse_complement(seq):
    seq = complement(seq)
    return reverse_seq(seq)

def main():
    f = open('covid-19-2-genome', 'r')
    covid_genome = f.read()
    f.close()

    # 30318, 1 byte per character, so 30KB
    print('\nLength of gene sequence:', len(covid_genome))

    # count bases (isn't covid-19 a single-strand RNA virus?)
    # so RNA should be A, C, G, and U instead of T?
    print('\nCounts:', count_base(covid_genome))

    # complement
    comp = complement(covid_genome)
    print('\nComplement:', comp[:10], '...')

    # reverse complement (switch to 5'--------3')
    reverse_comp = reverse_complement(covid_genome)
    print('\nReverse Complement:', reverse_comp[-11:-1], '...')

    assert(reverse_seq(reverse_comp) == comp)

    verify_str = 'tataaaa'
    if verify_str in reverse_comp:
        print('\nReverse complement computed correctly.')
    else:
        print('\nError. Reverse complement computed incorrectly. No existance of' + verify_str)


main()
    
