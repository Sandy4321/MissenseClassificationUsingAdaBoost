'''Computes sequence weights according to
S.J.Valder; 2002; Scoring Residue Conservation.

created: 21st Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import util
import check_mutation_position
import numpy as np


def do(sequences):

    L = len(sequences[0]) # length of aligned sequences
    N = len(sequences) # number of sequences

    weights = []
    for i in range(L):
	aa = check_mutation_position.do(sequences, i)
        freq = util.calc_frequency(aa)

        uniq_aa = util.unique(aa)

	# compute sequence weights
	w = []
	for i in range(N):
	    if aa[i] == '-' or aa[i] == 'X':
		w.append(0)
	    else:
		if aa[i] == 'B':
		    aa[i] = 'N'
		elif aa[i] == 'Z':
		    aa[i] = 'Q'
		w.append(1.0/(len(uniq_aa)*freq[aa[i]]))

	weights.append(w) # N x L matrix, for each position find the weights

    # compute average of w over all positions
    avg_weight = np.zeros(N)
    for i in range(L):
	avg_weight += np.array(weights[i])
    avg_weight *= 1.0/L

    return avg_weight
