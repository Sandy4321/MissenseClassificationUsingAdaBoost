'''Computes cumulative frequency of sliding neighbour.

created: 22nd Sept, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import util
import os
import sys
import load_alignments
import correct_alignment_position
import uniquify_alignments
from operator import itemgetter
from itertools import chain


# directory paths
aligndir = '../alignments'
fastadir = '../fasta'


def calc_hd_score(sps, qp, relpos, hd):
    '''Computes sliding window score.

    Args:
	sqp: sequences parts
	qp: query sequence parts
	hd: hamming distance (number of mis-hits)

    Return Values:
	conservation_score

    '''

    score = 0

    for seq in sps:
	shd = hd
	flag = 0
        for l in range(len(seq)):
	    if seq[l] != qp[l]:	# if there is a mismatch in character
		if l==relpos:	# if position is mutation position
		    flag = 1
		    break
		elif l!= relpos: # if some other position in the window
		    if shd == 0: # if number of mis-hits are exhausted
			flag = 1    # break
			break
		    else:
			shd -= 1 # else minus 1 from number of mis-hits	

	if flag:
	    continue
	else:
	    score += 1

    if len(sps) == 0:
	return 0

    return float(score)/len(sps)


def do(file1, file2, nbsize):
    mutations = []
    with open(file1) as fp:
	for line in fp:
	    mutations.append(line.strip().split(','))

    mutations.reverse()

    with open(file2,'w') as ofp:

	# ---------------------------------------------------------
	# write the heading to the output file

	heading = ['Protein','Position','Orig_aa','Mut_aa','Label']
	# neighbour part
	for n in range(1, nbsize+1):
	    for i in range(n+1):
		if n > 3:
		    maxhd = 2
		elif n > 2:
		    maxhd = 1
		else:
		    maxhd = 0

		for k in range(maxhd+1):
		    # hamming distance k
		    heading.append('nbwidth:%s srtpos:%s hd:%s orig' %(str(n), str(-i), str(k)))
		    heading.append('nbwidth:%s srtpos:%s hd:%s mut' %(str(n), str(-i), str(k)))

	#print heading
	ofp.write(','.join(heading)+'\n')

	# ---------------------------------------------------------
	# compute scores for mutations

	for mut in mutations:
	    protein = mut[0]
	    position = int(mut[1]) - 1
	    orig_aa = mut[2]
	    mut_aa = mut[3]

	    alignment_path = '%s/%s' %(aligndir, protein)
	    if not os.path.exists(alignment_path):
		continue

	    #print(mut)
	    query_flag = 0 # to capture IndexError in query protein sequence

	    # fetch alignments, and query protein
	    alignments = load_alignments.do(alignment_path)
	    query, alignments = util.fetch_query_protein_in_alignments(alignments)
	    alignments = util.prune_proteins_list(alignments)

	    # fetch all sequences
	    sequences = [a['alignment'] for a in alignments]
	    query_seq = query['alignment']
 
	    # all neighbour scores
	    scores = []

	    # go into neighbour mode now
	    for n in range(1, nbsize+1):

		if n > 3:
		    maxhd = 2
		elif n > 2:
		    maxhd = 1
		else:
		    maxhd = 0

		for i in range(n+1):

		    # i is the relative position of mutation position
		    # in nb_pos list

		    start_pos = position - i

		    # fetch the indices of neighbours for this round
		    nb_pos = []
		    try:
		    	for j in range(n+1):
			    # correct the position in the alignment
			    cpos = correct_alignment_position.do(query_seq, start_pos+j)
			    nb_pos.append(cpos)
		    except Exception:
			query_flag = 1
			break

		    # fetch corresponding positions of aligned sequences
		    seq_parts = []
		    for seq in sequences:
			s = list(seq)
			row = []
			flag = 0
			for pos in nb_pos:
			    try:
				# do correction for 'B' and 'Z' amino acids
				if s[pos] == 'B':
				    row.append('N')
				elif s[pos] == 'Z':
				    row.append('Q')
				else:
				    row.append(s[pos])
			    except IndexError:
				flag = 1
				break
			if flag:
			    continue
			seq_parts.append(row)

	   		if len(seq_parts) == 0:
			    query_flag = 1
			    break

		    if query_flag:
			break			
 
		    # fetch corresponding positions of query sequence
		    query_seq_part = []
		    s = list(query_seq)
		    for pos in nb_pos:
			try:
			    # do correction for 'B' and 'Z' amino acids
			    if s[pos] == 'B':
				query_seq_part.append('N')
			    elif s[pos] == 'Z':
				query_seq_part.append('Q')
			    else:
				query_seq_part.append(s[pos])
			except IndexError:
			    query_flag = 1
			    break

		    if query_flag:
			break
			    
		    # construct the mutated query sequence's parts
		    mut_query_seq_part = list(query_seq_part)
		    mut_query_seq_part[i] = mut_aa

		    for k in range(maxhd+1):
			#print('%s %s %s - original' %(str(n),str(i),str(k)))
			# compute hd-k orig score
			scores.append(calc_hd_score(seq_parts, query_seq_part, i, k))
			#print('%s %s %s - mutant' %(str(n),str(i),str(k)))
			# compute hd-k mut score
			scores.append(calc_hd_score(seq_parts, mut_query_seq_part, i, k))

		if query_flag:
		    break

	    if query_flag:
		continue

	    #print scores
	    ofp.write(','.join([ str(item) for item in list(chain.from_iterable([mut, scores])) ])+'\n')

	    #break

	    
		
def main():
    args = sys.argv

    if len(args) > 1:
        f1 = args[1] # where mutation info is present
        f2 = args[2] # where we store the scores
		nbsize = args[3]

    else:
        f1 = 'datasets/Training\ dataset\unique_total_training_valid.csv'
        f2 = 'features/training_slng_scores.csv'
		nbsize = 15

    do(f1, f2, nbsize)


if __name__ == '__main__':
    main()
