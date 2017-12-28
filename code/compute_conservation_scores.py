'''Combines computation of all different 
conservation scores into one program.

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
import shannon_entropy_score
import von_neumann_entropy_score
import sum_of_pairs_score


def do(file1, file2):
    # fetch all mutations
    mutations = []
    with open(file1) as fp:
        for line in fp:
            parts = line.strip().split(',')
            mutations.append([parts[0],parts[1],parts[2],parts[3],parts[4]])

    ofp = open(file2,'w')

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

        # prune the proteins 
        alignments = util.prune_proteins_list(alignments)

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
        
        # fetch the corresponding sequences
        sequences = [a['alignment'] for a in alignments]

        # compute frequency of amino acid at desired position
        aa = check_mutation_position.do([k['alignment'] for k in alignments], actual_pos)

        # compute simple frequency of original & mutant amino acid
        o_score, m_score = simple_frequency_score.do(list(aa), orig_aa, mut_aa)

        # compute score using pseudo-counts in order to account for missing aa
        o_ps_score, m_ps_score = pseudo_count_score.do(list(aa), orig_aa, mut_aa)

        # compute simple sequence-weighted frequency score
        sequence_weights = landgraf_sequence_weights.do([a['alignment'] for a in alignments])
        o_sw_score, m_sw_score = simple_frequency_score.do(list(aa), orig_aa, mut_aa, sequence_weights)
        # using gap frequencies
        o_gf_score, m_gf_score = gapped_frequency_score.do(list(aa), orig_aa, mut_aa, sequence_weights)

		# calculate shannon entropy score w/o sequence weights
        shannon = shannon_entropy_score.do(list(sequences), query['alignment'], actual_pos,  orig_aa, mut_aa)
        # shannon entropy score with sequence weights
        shannon_weighted = shannon_entropy_score.do(list(sequences), query['alignment'], actual_pos, orig_aa, mut_aa, 1)

        # calculate von-neumann entropy score
        von_neumann_score = von_neumann_entropy_score.do(list(aa), orig_aa, mut_aa)

	# calculate sum-of-pairs scores
        sop = sum_of_pairs_score.do(list(sequences), query['alignment'], actual_pos, orig_aa, mut_aa, 0) # wo seq wg
        sop_wg = sum_of_pairs_score.do(list(sequences), query['alignment'], actual_pos, orig_aa, mut_aa, 1) # w seq wg

        # append all scores together and write to file
        scores = [o_score, m_score, o_ps_score, m_ps_score, o_sw_score, m_sw_score, o_gf_score, m_gf_score, shannon, shannon_weighted, von_neumann_score, sop, sop_wg]
        #print(scores)

        ofp.write(','.join([ str(item) for item in list(chain.from_iterable([mut, scores])) ])+'\n')

        #break

    ofp.close()

def main():
    args = sys.argv

    if len(args) > 1:
        f1 = args[1] # where mutation info is present
        f2 = args[2] # where we store the scores

    else:
        f1 = 'datasets/training/valid_humvar.csv'
        f2 = 'datasets/training/humvar_conserv_scores.csv'

    do(f1,f2)


if __name__ == '__main__':
    main()
