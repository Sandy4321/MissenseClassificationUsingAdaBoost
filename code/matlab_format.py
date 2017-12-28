'''This script converts the generated
scores files into the format that can
be read by Matlab.

<scores, label>

author  : Rashmi Balasubramanyam
          <rashmi.bmanyam@gmail.com>
created : 7th jan, 2017
'''

#f1 = "../scores/training/total_training_conserv_scores_with_window_numseqs.csv"
#f2 = "../scores/training/total_training_conserv_scores_with_window_numseqs_Matlab.csv"
f1 = "../scores/test/TM_conserv_scores_with_window_numseqs.csv"
f2 = "../scores/test/TM_conserv_scores_with_window_numseqs_Matlab.csv"

#for i in [1,2,3,4,5]:
#    file1 = f1%i
#    file2 = f2%i
with open(f1) as fp, open(f2,"w") as ofp:
    for line in fp:
        parts = line.strip().split(",")
        reqd = parts[5:]
        label = parts[4]
        key = ",".join(reqd + [label])
        ofp.write(key+"\n")
