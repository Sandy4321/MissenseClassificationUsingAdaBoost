'''Computes Sum of Pairs scores,
as defined in capra2007predicting, 
and valdar2002scoring.

created: 21st Sept, 2015
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


def calculate_sop(aa, sw):
    blosum62_json = util.fetch_blosum62_json()

    N = len(aa)
    W = 0
    for i in range(N):
	if aa[i] == 'B':
	    aa[i] = 'N'
	elif aa[i] == 'Z':
	    aa[i] = 'Q'

    sp_score = 0
    for i in range(N):
	for j in range(i+1,N):
	    if aa[i]=='-' or aa[j]=='-' or aa[i]=='X' or aa[j]=='X':
		continue
	    W += sw[i]*sw[j]
	    sp_score += sw[i]*sw[j]*blosum62_json[aa[i]][aa[j]] 

    return float(sp_score)/W


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
    entropy1 = calculate_sop(patterns, sw)

    # case 2: for mutant amino acid
    sequences.pop()
    new_query_seq = util.strsub(query_seq, position, changed_pattern)
    sequences.append(new_query_seq)
    patterns = check_mutation_position.do(sequences, position)
    if SW:
        sw = landgraf_sequence_weights.do(sequences)
    else:
        sw = list(np.ones(len(sequences)))
    entropy2 = calculate_sop(patterns, sw)

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
    print(do(list(sequences), query['alignment'], mod_position, pattern, changed_pattern, 0))
    print(do(list(sequences), query['alignment'], mod_position, pattern, changed_pattern, 1))


if __name__ == '__main__':
    main()
