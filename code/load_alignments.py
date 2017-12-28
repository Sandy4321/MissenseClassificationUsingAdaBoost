'''This program loads the alignments into list of proteins.

created on: 7th May, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys

def do(filename):
    '''This method does all the hard work! This method loads the protein
    alignment from the filename passed as argument.

    Args:
	filename: of the file containing the alignments

    Return Values:
	alignment: of proteins
	proteins: list of proteins corresponding the sequences in alignment

    '''
    proteins = []
    seq = []
    row = {}
    iter = 0
	
    with open(filename) as fp:	
        for line in fp:
            if line[0] == '>':
                line = line.strip()
                parts = line.split(' ')
                row = {"type":parts[0][1:], "id":parts[1], "match_percentage":parts[2]}
                proteins.append(row)
                if seq:
                    proteins[iter]["alignment"] = ''.join(seq)
                    iter = iter + 1
                    seq = []
            else:
                seq.append(line.strip())
                # for the last sequence
                proteins[iter]["alignment"] = ''.join(seq)

    return proteins

if __name__ == '__main__':
    args = sys.argv

    if len(args) > 1:
        filename = args[1]
    else:
        filename = os.path.join(os.getcwd(), 'test/testalign')
        
    print(do(filename))
