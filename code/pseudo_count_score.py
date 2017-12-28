'''Adds pseudo-count of 1 to account for each type of
missing amino acid. Takes log to base 20 of the final
value computed.

created: 21st Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import util
from operator import itemgetter
from itertools import chain

def do(aa, orig_aa, mut_aa, PC={}):
    '''
    '''

    if len(PC) == 0:
	PC = {
		'A': 1,
		'R': 1,
		'N': 1,
	        'D': 1,
		'C': 1,
		'Q': 1,
		'E': 1,
		'G': 1,
		'H': 1,
		'I': 1,
		'L': 1,
		'K': 1,
		'M': 1,
		'F': 1,
		'P': 1,
		'S': 1,
		'T': 1,
		'W': 1,
		'Y': 1,
		'V': 1,
		'U': 1
	    }
    PC_total = sum([v for k,v in PC.iteritems()])

    for i in range(len(aa)):
	if aa[i] == 'X' or aa[i] == '-':
	    continue
	if aa[i] == 'B':
	    aa[i] = 'N'
	elif aa[i] == 'Z':
	    aa[i] = 'Q'

    freqs = util.calc_frequency(aa)

    N = sum([val for k,val in freqs.iteritems()]) # number of valid sequences

    probs = {}
    for k,v in freqs.iteritems():
	probs[k] = float(v)/N

    o_ps_score = (float(freqs[orig_aa]) + PC[orig_aa]) / (N + PC_total)
    m_ps_score = (float(freqs[mut_aa]) + PC[mut_aa]) / (N + PC_total)

    return o_ps_score, m_ps_score


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
