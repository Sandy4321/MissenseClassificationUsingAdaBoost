'''Contains commonly used methods.

created on: 18th Sept, 2015
author: Rashmi Balasubramyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import json
import uniquify_alignments
from operator import itemgetter
from itertools import chain


def unique(arr):
    '''returns unique elements in array.

    Args:
	arr: array containing duplicates

    Return Values:
	unique array elements

    '''

    return list(set(arr))


def fetch_query_protein_in_alignments(proteins):
    '''Fetches the query protein, deletes it from the
    list of protein alignments passed as input.

    Args:
        proteins: Protein alignments in following format
                <id, type, match_percentage, alignment>

    Return Values:
        query: query protein
        proteins: list of aligned proteins - query protein

    '''

    index = 0    #will refer to the record of query protein
    for elem in proteins:
        if elem['type'] == 'query':
            break
        else:
            index = index + 1
    query = proteins[index]
    del proteins[index]

    return query, proteins


def read_sequence(filename):
    '''Reads fasta file and stores the data in an array.

    Args:
    filename: where the fasta sequence is stored

    Return Value:
    seq: sequence as array of lines in the fasta format
        if the fasta sequence file was found to be empty,
        it returns nothing.

    '''

    seq = []
    try:
        with open(filename) as fp:
            for line in fp:
                seq.append(line.strip())
        del seq[0]
        return ''.join(seq)
    except Exception as e:
        print(str(e))
        pass


def calc_frequency(aa):
    '''This file creates a hash containing the count of the
    amino acids at the mutation position of alignment of protein
    sequences.

    Args:
    aa: list containing all amino acids at the mutation position
        of the aligned protein sequences

    Return Value:
    freq: hash map of 'amino_acid': count

    '''

    freq = {'A': 0,
    'R': 0,
    'N': 0,
    'D': 0,
    'C': 0,
    'Q': 0,
    'E': 0,
    'G': 0,
    'H': 0,
    'I': 0,
    'L': 0,
    'K': 0,
    'M': 0,
    'F': 0,
    'P': 0,
    'S': 0,
    'T': 0,
    'W': 0,
    'Y': 0,
    'V': 0,
    'U': 0,
    '-': 0
    }

    try:
        for a in aa:
            # ignore if 'X' character
            if a == 'X':
                continue

            elif a == 'B':
                freq['N'] = freq['N'] + 1
	    elif a == 'Z':
                freq['Q'] = freq['Q'] + 1

            else:
                freq[a] = freq[a]+1

        return freq

    except Exception as e:
        print str(e)
        pass


def calc_frequency_2(aa):
    '''This file creates a hash containing the count of the
    amino acids at the mutation position of alignment of protein
    sequences. This method computes frequency of "-" as that of 
    "-", "X", "B" and "Z". [Zhang et. al. 2008, Estimating residue evolutionary conservation by introducing von Neumann entropy and a novel gap-treating approach]

    Args:
    aa: list containing all amino acids at the mutation position
        of the aligned protein sequences

    Return Value:
    freq: hash map of 'amino_acid': count

    '''

    freq = {'A': 0,
    'R': 0,
    'N': 0,
    'D': 0,
    'C': 0,
    'Q': 0,
    'E': 0,
    'G': 0,
    'H': 0,
    'I': 0,
    'L': 0,
    'K': 0,
    'M': 0,
    'F': 0,
    'P': 0,
    'S': 0,
    'T': 0,
    'W': 0,
    'Y': 0,
    'V': 0,
    'U': 0,
    '-': 0
    }

    try:
        for a in aa:
            # ignore if 'X' character
            if a == 'X' or a == '-' or a == 'B' or a == 'Z':
		freq['-'] += 1

            else:
                freq[a] += 1

        return freq

    except Exception as e:
        print str(e)
        pass


def fetch_blosum62_json():
    '''
    '''

    with open('blosum62.json') as fp:
	obj = json.load(fp)

    blosum62_json = {}
    for k1,v1 in obj.iteritems():
	row = {}
	for k2,v2 in v1.iteritems():
	    row[k2] = float(v2)
	blosum62_json[k1] = row

    return blosum62_json


def fetch_blosum62_matrix():
    '''
    '''

    blosum62_json = fetch_blosum62_json()

    blosum62_matrix = []
    order = []

    for k,v in blosum62_json.iteritems():
	if k=='-' or k=='X' or k=='B' or k=='Z':
	    continue
	row = []
	for k1,v1 in v.iteritems():
	    if k1=='-' or k1=='X' or k1=='B' or k1=='Z':
		continue
	    row.append(float(v1))
	blosum62_matrix.append(row)
	order.append(k)

    return blosum62_matrix, order


def strsub(word, position, newchar):
    word = list(word)
    word[position] = newchar
    new_word = ''.join(word)
    return new_word


def prune_proteins_list(proteins, percent=None):
    '''Let us remove alignments that have match percentage BELOW the
    desired percent.

    Args:
        percent: cut-off percentage

    Return Value:
    proteins: that have match percentage above the desired cut-off

    '''

    # remove duplicate sequences in aligned sequences
    proteins = uniquify_alignments.do(proteins)

    # apply some filtering on proteins here
    # sorted_result = sorted(proteins, key=operator.itemgetter(1))
    sorted_list = sorted(proteins, key=itemgetter('match_percentage'))
    proteins = sorted_list
    proteins.reverse()

    # TBD: filter them based on a cut-off
    if not percent:
	percent = 30.0
    
    reqd_proteins = [prot for prot in proteins if float(prot['match_percentage'])>percent]

    return reqd_proteins
