'''This program computes the correct position number of
query sequence in the alignments.

created on: 17th June, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys


def do(sequence, position):
    '''Corrects the value of desired position
    of the query sequence in the sequences passed
    as parameter.

    Args:
        sequence: aligned sequence with gaps in which 
	    we want to find letter of desired position
        position: desired position of query protein

    Return Value:
        mod_position: modified position

    '''

    position = position + 1

    # compute the real position in the aligned sequence 
    for i in range(len(sequence)):
        if not (sequence[i] == "-"):
            position = position - 1
            if position == 0:
                mod_position = i
                break
    
    return mod_position

def main():
    args = sys.argv

    if len(args) > 1:
        # alignment filename
        filename = args[1]
        position = int(args[2])

    else:
        filename = "P26439_short.align"
        position = 60

    import load_alignments
    sequences = load_alignments.do(filename)

    print(do(sequences, position))        

if __name__ == "__main__":
    main()
