'''Computes window based score
as the weighted sum of individual
scores of neighbouring positions.

created: 22nd Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import load_alignments
import util
import uniquify_alignments
import correct_alignment_position
import check_mutation_position
from operator import itemgetter
from itertools import chain
import simple_frequency_score
import pseudo_count_score
import landgraf_sequence_weights
import gapped_frequency_score

def do(file1, file2, winsize=1):
    # fetch all mutations
    mutations = []
    with open(file1) as fp:
        for line in fp:
            parts = line.strip().split(',')
            mutations.append([parts[0],parts[1],parts[2],parts[3],parts[4]])

#    ofp = open(file2,'w')

    for mut in mutations:
        # fetch reqd info
        protein = mut[0]
        position = int(mut[1]) - 1
        orig_aa = mut[2]
        mut_aa = mut[3]

        # if protein not aligned - pass this
        alignment_file = '../alignments/%s'%protein
        if not os.path.isfile(alignment_file):
            continue

        print(mut)

        # load the aligned sequences
        alignments = load_alignments.do(alignment_file)

        # remove duplicate sequences in aligned sequences
        alignments = uniquify_alignments.do(alignments)

        # fetch query sequence
        query, alignments = util.fetch_query_protein_in_alignments(alignments)
        if len(alignments) == 0:
            continue

        # fetch actual position of mutation in aligned query sequence
        try:
            actual_pos = correct_alignment_position.do(query['alignment'], position)
        except Exception as e:
            print(query['alignment'])
            continue

	o_win_score = 0
	for w in range(winsize+1):
	    
        
	
	scores = [o_score, m_score, o_ps_score, m_ps_score, o_sw_score, m_sw_score, o_gf_score, m_gf_score]
        print(scores)

	#ofp.write(','.join([ str(item) for item in list(chain.from_iterable([mut, scores])) ])+'\n')

        break

    #ofp.close()

def main():
    args = sys.argv

    if len(args) > 1:
        f1 = args[1] # where mutation info is present
        f2 = args[2] # where we store the scores
	winsize = args[3]

    else:
        f1 = 'datasets/test/unseen_proteins.csv'
        f2 = 'datasets/test/unseen_proteins_freq_scores.csv'
	winsize = 3

    do(f1, f2, winsize)


if __name__ == '__main__':
    main()
