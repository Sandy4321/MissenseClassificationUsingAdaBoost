'''Distributes frequency of gaps among
the 20 amino acids.

created: 21th Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import util
from operator import itemgetter
from itertools import chain
import numpy as np


def do(aa, orig_aa, mut_aa, sw=[]):
    '''
    '''

    # optional parameter sw: if not set - set it to ones
    if len(sw) == 0:
        sw = list(np.ones(len(aa)))

    N = len(aa) # number sequences
    W = 0 # sum of the weights of sequences w/o 'X'

    # to account for query amino acid
    o_gf_score = 0
    m_gf_score = 0

    for i in range(N):
	if aa[i] == 'X':
	    continue

	elif aa[i] == 'B':
	    aa[i] = 'N'
	elif aa[i] == 'Z':
	    aa[i] = 'Q'

	elif aa[i] == '-':
	    o_gf_score += float(sw[i])/20
	    m_gf_score += float(sw[i])/20

	if aa[i] == orig_aa:
	    o_gf_score += sw[i]

	elif aa[i] == mut_aa:
	    m_gf_score += sw[i]

	W += sw[i]

    o_gf_score = float(o_gf_score)/W
    m_gf_score = float(m_gf_score)/W

    return o_gf_score, m_gf_score


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

    patterns = check_mutation_position.do([prot['alignment'] for prot in proteins], mod_position)

    print(do(patterns, pattern, changed_pattern))


if __name__ == '__main__':
    main()
