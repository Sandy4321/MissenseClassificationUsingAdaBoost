''' This module computes conservation score as the difference
in the Shannon Entropy values after the mutation.

created on: 26th July, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import multiprocessing as mp
# user-defined modules
import check_mutation_position
import load_alignments
import shannon_entropy_score
import von_neumann_entropy_score
import sum_of_pairs_score
import correct_alignment_position
from operator import itemgetter
from itertools import chain
import uniquify_alignments
import util

CHUNK_SIZE = 100

align_dirs= [os.path.join(os.path.pardir,'alignments/')]
#align_dirs = ['alignments/', 'taxonomy/humans/alignments', 'taxonomy/rodents/alignments', 'taxonomy/mammals/alignments', 'taxonomy/vertebrates/alignments']


def prune_proteins_list(proteins, percent=None):
    '''Let us remove alignments that have match percentage BELOW the
    desired percent.

    Args:
        percent: cut-off percentage

    Return Value:
    proteins: that have match percentage above the desired cut-off

    '''

    proteins = uniquify_alignments.do(proteins)

    # apply some filtering on proteins here
    #    sorted_result = sorted(proteins, key=operator.itemgetter(1))
    sorted_list = sorted(proteins, key=itemgetter('match_percentage'))
    proteins = sorted_list
    proteins.reverse()

    # TBD: filter them based on a cut-off
    if percent:
        print('List will be filtered on percent')

    return proteins


def work(line):
    '''This method computes all kinds of scores,
    for a given mutation, and returns the computed scores.

    Args:
    line: which contains mutation info.
        <PROTEIN_ID, POSITION, ORIGINAL_AMINO_ACID, MUTANT_AMINO_ACID>

    Return Values:
    [parts, scores] flattend list

    '''

    line = line.strip()
    parts = line.split(',')

    protein = parts[0]
    position = int(parts[1]) - 1 # correction for indexing python lists
    mut_aa = parts[3]

    result = []

    print(protein)

    for alignment_dir in align_dirs:

        alignment_file = os.path.join(alignment_dir, protein)
        if not os.path.isfile(alignment_file):
            return []

        fasta = util.read_sequence('../fasta/%s.fasta' %protein)

        proteins = load_alignments.do(alignment_file) # types of proteins required

        # sort & prune the list of proteins
        proteins = prune_proteins_list(proteins)

        # fetch the record pertaining to the query protein
        query, p = util.fetch_query_protein_in_alignments(proteins) # p does not contain query sequence

        prots = [ prot["id"] for prot in proteins ]
        types = [ prot["type"] for prot in proteins ]
        match_percents = [ prot["match_percentage"] for prot in proteins ]
        alignments = [ prot["alignment"] for prot in proteins ]
        proteins = prots   #replacing proteins dicitonary with only ids

        result.append(len(prots)+1)

        orig_aa = fasta[position]
        # correct the position of query protein wrt to alignment
        mod_pos = correct_alignment_position.do(query['alignment'], position)
        aa = check_mutation_position.do(alignments, mod_pos)
        
	# calculate shannon entropy score w/o sequence weights
        result.append(shannon_entropy_score.do(list(alignments), query['alignment'], mod_pos,  orig_aa, mut_aa))
	# shannon entropy score with sequence weights
        result.append(shannon_entropy_score.do(list(alignments), query['alignment'], mod_pos, orig_aa, mut_aa, 1))

	# calculate von-neumann entropy score
	result.append(von_neumann_entropy_score.do(list(aa), orig_aa, mut_aa))

	# calculate relative entropy score
	#result.append(relative_entropy_score.do(list(aa), orig_aa, mut_aa))

	# calculate jensen-shannon divergence score
	#result.append(jensen_shannon_divergence_score(list(aa), orig_aa, mut_aa))

	# calculate sum-of-pairs scores
	result.append(sum_of_pairs_score.do(list(alignments), query['alignment'], mod_pos, orig_aa, mut_aa, 0)) # wo seq wg
	result.append(sum_of_pairs_score.do(list(alignments), query['alignment'], mod_pos, orig_aa, mut_aa, 1)) # w seq wg

    # return mutation information and scores for recording in a file
    return [ str(item) for item in list(chain.from_iterable([parts, result])) ]


def do(file1, file2):
    '''This method does single/multi processing.

    Args:
    file1: file which contains the mutations
    file2: file where we write the result

    Return Values:
    returns nothing.

    '''

    iter = 0
    res = []

    with open(file1) as ifp: #, open(file2, 'w') as ofp:
        # create the header to the result file
        fields = ["Protein ID", "position", "original AA", "mutated AA", "label"]

        for i in range(len(align_dirs)):
            fields.append("num_seqs")
            fields.append("shannon_entropy_diff")
            fields.append("vonneumann_entropy_diff")

        # write the header to the result file
        #ofp.write(','.join(fields)+'\n')

#       # using multiprocessing
#       pool = mp.Pool()
#       for res in pool.imap_unordered(work, ifp, CHUNK_SIZE):
#           if result:
#               ofp.write(','.join(res)+'\n')

        # using single processing
        for line in ifp:
            res = []
            try:
                #iter = iter + 1
                res = work(line)
		print res
		break
                #if res:
                #    ofp.write(','.join(res)+'\n')
                #if iter == 2:
		#    break
            except ZeroDivisionError as e:
                print(line)
                pass
            except IndexError as ie:
                print(line)
                print("index error")
                pass
            except UnboundLocalError as ue:
                print(line)
                print('Unbounded Local Error')
                pass
            #except TypeError as te:
            #    print(line)
            #    print('Human, etc not present for this protein')
            #    pass
            except Exception as ee:
                print(res)
                print(line)
                print(str(ee))
                raise ee

if __name__ == "__main__":
    args = sys.argv

    if len(args) > 1:
        file1 = args[1]
        file2 = args[2]

    else:
        file1 = os.path.join(os.getcwd(), 'datasets/test/unseen_proteins.csv')
        file2 = os.path.join(os.getcwd(), 'datasets/test/unseen_proteins_entropy_scores.csv')

    do(file1, file2)
