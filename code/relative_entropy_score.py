'''This function computes conservation score
as a difference in relative entropy scores.

created: 20th Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import util
import json
import numpy as np


def calculate_entropy(patterns):
    freq = util.calc_frequency(patterns)
    if '-' in freq:
        del freq['-']


def do(patterns, pattern, changed_pattern):
    '''Finds the difference in the entropy in patterns
    when a pattern is changed to changed pattern.

    Args:
        patterns: list of patterns
        pattern & changed_pattern: whose entropy diff
                    we want to calculate.

    Return Value:
        difference in entropy score

    '''

    # case 1: for original amino acid
    patterns.append(pattern)
    entropy1 = calculate_entropy(patterns)

    # case 2: for mutant amino acid
    patterns.pop() # change the original amino acid to mutant type
    patterns.append(changed_pattern)
    entropy2 = calculate_entropy(patterns)

    return (entropy2 - entropy1)


def main():
    filename = 'alignments/Q92481'
    position = 72
    changed_pattern = 'R'

    import load_alignments
    import uniquify_alignments
    import check_mutation_position
    import correct_alignment_position

    proteins = load_alignments.do(filename)
    sorted_list = sorted(proteins, key=itemgetter('match_percentage'))
    proteins = sorted_list
    proteins.reverse()

    proteins = uniquify_alignments.do(proteins)
    query,prots = util.fetch_query_protein_in_alignments(proteins)

    mod_position = correct_alignment_position.do(query['alignment'], position)
    pattern = query['alignment'][mod_position]

    patterns = check_mutation_position.do([prot['alignment'] for prot in proteins], mod_position)
    print(do(patterns, pattern, changed_pattern))


if __name__ == '__main__':
    main()
