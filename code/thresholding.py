'''This method finds the best threshold in a binary classifier,
such that the system has optimum Specificity and Sensitivity.

created on: 10th June, 2015
author: Rashmi Balasubramanyam <rashmi.bmanyam@gmail.com>

'''

# required packages and modules
import os
import sys
import numpy as np
import util

def calc_metrics(data, thresh):
    '''
    '''

    tp = 0
    tn = 0
    fp = 0
    fn = 0
    
    for item in data:
	if item[0] >= thresh:
	    ypred = 1
	else:
	    ypred = 0

	if ypred == 1:
	    if item[1] == 1:
		tp = tp + 1
	    else:
		fp = fp + 1
	else:
	    if item[1] == 1:
		fn = fn + 1
	    else:
		tn = tn + 1

    return tp, tn, fp, fn


def find_optimum_threshold(X, y):
    '''This method finds optimum threshold as required.

    Args:
	X: score
	y: label

    Return Values:
	threshold: the optimum value as desired

    '''

    # save the data into 2-d array so that we can sort them
    # on the X values
    data = np.array([X, y])
    data = np.transpose(data)

    # metrics
    tps = []
    tns = []
    fps = []
    fns = []
    thresholds = []
    allthresh = util.unique(list(data[:,0]))
    allthresh.sort()

    thresh = allthresh[0] - 0.000001
    tp, tn, fp, fn = calc_metrics(data, thresh)
    tps.append(tp)
    tns.append(tn)
    fps.append(fp)
    fns.append(fn)
    thresholds.append(thresh)
   
 
    for i in range(len(allthresh)-1):
	thresh = (allthresh[i] + allthresh[i+1]) / 2
	tp, tn, fp, fn = calc_metrics(data, thresh)
	tps.append(tp)
	tns.append(tn)
	fps.append(fp)
	fns.append(fn)
	thresholds.append(thresh)

    thresh = allthresh[-1] + 0.000001
    tp, tn, fp, fn = calc_metrics(data, thresh)
    tps.append(tp)
    tns.append(tn)
    fps.append(fp)
    fns.append(fn)
    thresholds.append(thresh)


    return thresholds, tps, tns, fps, fns


def do(file1, file2, index):
    '''This method is a wrapper method, that computes the threshold.

    Args:
	file1: containing the data and classification labels
	file2: we write metrics

    Return Value:
	threshold: value which produces optimum Specificity and Sensitivity
    
    '''

    X = []
    y = []

    # load the data
    with open(file1) as fp:
	lines = fp.readlines()
	
	for line in lines:
	    # data is in this format
	    # we are interested in label, and score in desired index

	    line = line.strip()
	    parts = line.split(",")

	    X.append(float(parts[index]))
	    y.append(int(parts[-1]))

	thresholds, tps, tns, fps, fns = find_optimum_threshold(X, y)

    with open(file2, 'w') as ofp:
	for i in range(len(thresholds)):
	    ofp.write("%s,%s,%s,%s,%s\n" %(thresholds[i],tps[i],tns[i],fps[i],fns[i]))


def main():
    args = sys.argv

    if len(args) > 1:
	file1 = args[1]
	file2 = args[2]
	index = int(args[3])

    else:
	file1 = "datasets/training/humvar_conserv_scores_matlab.csv"
	file2 = "datasets/training/consscores/humvar_sf_orig.csv"
	index = 0 

    print do(file1, file2, index)

if __name__ == "__main__":
    main()
