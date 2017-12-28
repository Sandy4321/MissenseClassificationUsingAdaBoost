'''This method takes in the alignment and proteins list,
and counts the frequency of each of the 20 amino acids at 
the mutation position.

created on: 7th May, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys

def aa_in_pos(alignment, position):
    '''This method returns the list of all amino acids
    (including gap) at the mutation position of the 
    aligned sequences.

    Args:
	alignment: of matched protein sequences
	position: of interest (mutation position)

    Return Values:
	aa: from all sequences, present in mutation position

    '''
    
    aa = []
    for line in alignment:
	line = line.strip()
	aa.append(line[position])

    return aa

    
def do(alignment, position):
    '''Performs the counting operation.

    Args:
	alignment: of matching protein sequences
	position: of mutation in query protein

    Return Value:
	freq: returned by "calc_frequency" method

    '''

    aa = aa_in_pos(alignment, position)

    return aa

if __name__ == '__main__':
    args = sys.argv

    if len(args) > 1:
	alignment = args[1]
	position = int(args[2])
    
    else:
	import load_alignments
	alignment = load_alignments.do(os.path.join(os.getcwd(), 'hh'))[0]
	position = 20
   
    print do(alignment, position)
