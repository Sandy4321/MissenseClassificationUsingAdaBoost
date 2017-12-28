'''Computes conservation scores of the neighbours.

created: 29th Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import util
import check_mutation_position
import correct_alignment_position
import load_alignments
from operator import itemgetter
from itertools import chain
import landgraf_sequence_weights
import gapped_frequency_score


def do(file1, file2, winsize):
    mutations = []
    with open(file1) as fp:
	for line in fp:
	    mutations.append(line.strip().split(','))

    with open(file2, 'w') as ofp:
        for mut in mutations:

	    print mut
	    flag = 0

	    protein = mut[0]
	    position = int(mut[1]) - 1
	    orig_aa = mut[2]
	    mut_aa = mut[3]

	    fasta_seq = util.read_sequence('../fasta/%s.fasta'%protein)

	    alfile = '../alignments/%s'%protein
	    if not os.path.isfile(alfile):
		continue

	    proteins = load_alignments.do(alfile)
	    proteins = util.prune_proteins_list(proteins)
	    query,proteins = util.fetch_query_protein_in_alignments(proteins)
	    query_seq = query['alignment']

	    alignments = [a['alignment'] for a in proteins]
	    try:
		sequence_weights = landgraf_sequence_weights.do(alignments)
	    except Exception as e:
		print str(e)
		flag = 1

	    if flag:
		continue

	    scores = []
    
	    for w in range(winsize+1): # 0,1,2,3 for winsize=3, hence the +1
		try:
		    # what score to use?
		
		    if w==0:
			mod_pos = correct_alignment_position.do(query_seq, position)
			aa = check_mutation_position.do(alignments, mod_pos)
			o, m =  gapped_frequency_score.do(list(aa), orig_aa, mut_aa, sequence_weights) # mod_pos
			scores.append(m)

		    else:
			mod_pos = correct_alignment_position.do(query_seq, position-w)
			aa = check_mutation_position.do(alignments, mod_pos)
			aa_in_fasta = fasta_seq[position-w]
			o, m =  gapped_frequency_score.do(list(aa), aa_in_fasta, mut_aa, sequence_weights) # mod_pos
			scores.append(o) # left neighbour at position w from mutation position

			mod_pos = correct_alignment_position.do(query_seq, position+w)
			aa = check_mutation_position.do(alignments, mod_pos)
			aa_in_fasta = fasta_seq[position+w]
			o, m =  gapped_frequency_score.do(list(aa), aa_in_fasta, mut_aa, sequence_weights) # mod_pos
			scores.append(o) # right neighbour at position w from mutation position

		except Exception:
		    flag = 1
		    break

	    if flag:
		continue

	    #print(scores)

	    ofp.write(','.join([ str(item) for item in list(chain.from_iterable([mut, scores])) ])+'\n')	
	    #break

	
def main():
    args = sys.argv

    if len(args) > 1:
	file1 = args[1]
	file2 = args[2]
	winsize = args[3]

    else:
	file1 = 'datasets/test/unseen_mutations.csv'
	file2 = 'datasets/test/temp.csv'
	winsize = 3

    do(file1,file2,winsize)

if __name__ == '__main__':
    main()
