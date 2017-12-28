'''Computes optimum threshold from the result of
thresholding.py method, and evaluates this against
the given list of datasets.

We use G-means for computing optimum threshold.

created: 2nd Oct, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import util
import math


def do(file1, file2, file3, index, direction=1):
    data = []
    with open(file1) as fp:
	for l in fp:
	    p = l.strip().split(',')
	    k = [float(a) for a in p]
	    data.append(k)

    # data[0] = threshold value
    # data[1] = tp
    # data[2] = tn
    # data[3] = fp
    # data[4] = fn

    max_gmean = 0
    opt_idx = 0

    for i in range(len(data)):
	d = data[i]

	gmean = math.sqrt((d[1]*d[2])/((d[1]+d[4])*(d[2]+d[3])))
	if gmean > max_gmean:
	    max_gmean = gmean
	    opt_idx = i

    optimum = data[opt_idx]
    opt_threshold = optimum[0]

    print "Training data"
    print "Optimum Threshold: %s" %opt_threshold
    print "TP: %s" %optimum[1]
    print "TN: %s" %optimum[2]
    print "FP: %s" %optimum[3]
    print "FN: %s" %optimum[4]

    # test set 1 
    vX = []
    vY = []
    with open(file2) as fp:
	for l in fp:
	    p = l.strip().split(',')
	    vX.append(float(p[index]))
	    vY.append(int(p[-1]))
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for i in range(len(vX)):
	if direction == 1:
	    if vX[i] < opt_threshold:
		if vY[i] == 1:
		    tp += 1
		else:
		    fp += 1
	    else:
		if vY[i] == 1:
		    fn += 1
		else:
		    tn += 1
	else:
	    if vX[i] >= opt_threshold:
		if vY[i] == 1:
		    tp += 1
		else:
		    fp += 1
	    else:
		if vY[i] == 1:
		    fn += 1
		else:
		    tn += 1
    print "Test Set: 1"
    print "TP: %s" %tp
    print "TN: %s" %tn
    print "FP: %s" %fp
    print "FN: %s" %fn

    # test set 2
    tX = []
    tY = []
    with open(file3) as fp:
        for l in fp:
            p = l.strip().split(',')
            tX.append(float(p[index]))
            tY.append(int(p[-1]))
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for i in range(len(tX)):
        if direction == 1:
            if tX[i] < opt_threshold:
                if tY[i] == 1:
                    tp += 1
                else:
                    fp += 1
            else:
                if tY[i] == 1:
                    fn += 1
                else:
                    tn += 1
        else:
            if tX[i] >= opt_threshold:
                if tY[i] == 1:
                    tp += 1
                else:
                    fp += 1
            else:
                if tY[i] == 1:
                    fn += 1
                else:
                    tn += 1
    print "Test Set: 2"
    print "TP: %s" %tp
    print "TN: %s" %tn
    print "FP: %s" %fp
    print "FN: %s" %fn

 

def main():
    args = sys.argv
  
    if len(args) > 1:
	file1 = args[1]	# data from which we find optimum threshold
	file2 = args[2]	# test set 1
	file3 = args[3]	# test set 2
	index = int(args[4])
	if len(args) == 6:
	   direction = int(args[5])

	do(file1, file2, file3, index, direction)

    else:
	print "Please enter threshold's filename, training data filename, and 2 test data filenames"
	exit


if __name__ == '__main__':
    main()
