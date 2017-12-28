'''This script generates difference in
physico-chemical properties.

values.csv has values for each amino acid
and the mean of the values in that
alignment column.

17th Jan 2017.
'''

mutfile = "../datasets/training/unique_total_training_valid.csv"
destfile = "../scores/traiining/total_training_physco2.csv"
# file containing physico-chemical property values for each amino acid
f1 = "../scores/physco_values.csv"

# load the physico-chemical properties
import json
values = json.load(open(f1))
properties = ["Bulkiness","Compressibility","Flexibility0","Flexibility1","Flexibility2","Helix_probability","Hydrophobicity","Ionisation","Isoelectric_point","Polarisability","Sheet_probability","Steric","Volume"]

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

with open(mutfile) as fp, open(destfile,"w") as ofp:
    for line in fp:
        print(line)

        parts = line.strip().split(',')
        # fetch reqd info
        protein = mut[0]
        position = int(mut[1]) - 1
        orig_aa = mut[2]
        mut_aa = mut[3]

        # if protein not aligned - pass this
        alignment_file = '../alignments/%s'%protein
        if not os.path.isfile(alignment_file):
            continue

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

        # compute mean value in the column
        mean = {}
        cnt = 0
        for key in values.keys(): # for each physco property
            first = 1
            for z in aa:
                if first:
                    if z != '-':
                        mean[key] = values[key][z]
                        cnt += 1
                        first = 0
                    else:
                        continue
                else:
                    if z != '-':
                        mean[key] += values[key][z]
                        cnt += 1
                    else:
                        continue
                mean[key] = mean[key] / cnt # compute the mean value

        # compute values of mutant
        mutvalue = {}
        for key in values.keys():
            mutvalue[key] = values[key][mut_aa]

        # compute the difference
        scores = []
        for key in properties:
            scores.append((mean[key] - mutvalue[key])*1.0/mean[key])

        row = ",".join(parts + [str(s) for s in scores]) + "\n"
        ofp.write(row)
