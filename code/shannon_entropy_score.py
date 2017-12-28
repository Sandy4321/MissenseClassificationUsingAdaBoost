'''Computes the difference in Shannon Entropy.

created on: 26th June, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import math
from operator import itemgetter
from itertools import chain
import util
import numpy as np
import check_mutation_position
import landgraf_sequence_weights


def calculate_entropy_with_gaps(aa, sw = []):
    '''
    '''

    N = len(aa)
    W = 0 # effective weight of all sequences

    if len(sw) == 0:
        sw = list(np.ones(N))

    prob = {'A': 0,
    'R': 0,
    'N': 0,
    'D': 0,
    'C': 0,
    'Q': 0,
    'E': 0,
    'G': 0,
    'H': 0,
    'I': 0,
    'L': 0,
    'K': 0,
    'M': 0,
    'F': 0,
    'P': 0,
    'S': 0,
    'T': 0,
    'W': 0,
    'Y': 0,
    'V': 0,
    'U': 0,
    '-': 0
    }

    for i in range(N):
        if aa[i] == '-' or aa[i] == 'X':
            continue

        elif aa[i] == 'B':
	    aa[i] = 'N'
	elif aa[i] == 'Z':
            aa[i] = 'Q'

        prob[aa[i]] += sw[i]
        W += sw[i]

    # sharing gap frequency
    gap_freq = float(prob['-']) / 20
    del prob['-']

    entropy = 0
    for k,v in prob.iteritems():
        prob[k] = float(v + gap_freq) / W
        if prob[k] == 0:
            continue
        entropy -= prob[k] * math.log(prob[k], 2)

    return entropy


def calculate_entropy(aa, sw = []):

    N = len(aa)
    W = 0 # effective weight of all sequences
    
    if len(sw) == 0:
	sw = list(np.ones(N))

    prob = {'A': 0,
    'R': 0,
    'N': 0,
    'D': 0,
    'C': 0,
    'Q': 0,
    'E': 0,
    'G': 0,
    'H': 0,
    'I': 0,
    'L': 0,
    'K': 0,
    'M': 0,
    'F': 0,
    'P': 0,
    'S': 0,
    'T': 0,
    'W': 0,
    'Y': 0,
    'V': 0,
    'U': 0
    }

    for i in range(N):
	if aa[i] == '-' or aa[i] == 'X':
	    continue

	elif aa[i] == 'B':
	    aa[i] = 'N'
	elif aa[i] == 'Z':
	    aa[i] = 'Q'

	prob[aa[i]] += sw[i]
	W += sw[i]

    entropy = 0
    for k,v in prob.iteritems():
	if prob[k] == 0:
	    continue
	prob[k] = float(v) / W
	entropy -= prob[k] * math.log(prob[k], 2)

    return entropy


def do(sequences, query_seq, position, pattern, changed_pattern, SW = 0):
    '''Finds the difference in the entropy in patterns
    when a pattern is changed to changed pattern.
    
    '''

    # case 1: for original amino acid
    sequences.append(query_seq)
    patterns = check_mutation_position.do(sequences, position)
    if SW:
	sw = landgraf_sequence_weights.do(sequences)
    else:
	sw = list(np.ones(len(sequences)))
    entropy1 = calculate_entropy(patterns, sw)

    # case 2: for mutant amino acid
    sequences.pop() 
    new_query_seq = util.strsub(query_seq, position, changed_pattern)
    sequences.append(new_query_seq)
    patterns = check_mutation_position.do(sequences, position)
    if SW:
        sw = landgraf_sequence_weights.do(sequences)
    else:
	sw = list(np.ones(len(sequences)))
    entropy2 = calculate_entropy(patterns, sw)

    return (entropy2 - entropy1)


def main():
    filename = '../alignments/Q92481'
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

    sequences = [prot['alignment'] for prot in prots]

    patterns = check_mutation_position.do(sequences, mod_position)
    print(do(sequences, query['alignment'], mod_position, pattern, changed_pattern, 1))


if __name__ == '__main__':
    main()

