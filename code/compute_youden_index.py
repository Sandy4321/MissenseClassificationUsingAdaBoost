'''This script computes
the best Youden index from
the scores.

18th Jan, 2017
'''

# -----------------------------------------------------------------------------
# inputs

# conserv scores as features

'''
infiles = ["../scores/training/cs/o_gf.csv","../scores/training/cs/m_sf.csv","../scores/training/cs/o_ps.csv","../scores/training/cs/m_ps.csv","../scores/training/cs/o_sw.csv","../scores/training/cs/m_sw.csv","../scores/training/cs/o_gf.csv","../scores/training/cs/m_gf.csv","../scores/training/cs/shannon.csv","../scores/training/cs/shannon_wg.csv","../scores/training/cs/vonNeumann.csv","../scores/training/cs/sop.csv","../scores/training/cs/sop_wg.csv","../scores/training/cs/windows1.csv","../scores/training/cs/windows2.csv"]
outfiles = ["../scores/training/cs/o_gf_ROC.csv","../scores/training/cs/m_sf_ROC.csv","../scores/training/cs/o_ps_ROC.csv","../scores/training/cs/m_ps_ROC.csv","../scores/training/cs/o_sw_ROC.csv","../scores/training/cs/m_sw_ROC.csv","../scores/training/cs/o_gf_ROC.csv","../scores/training/cs/m_gf_ROC.csv","../scores/training/cs/shannon_ROC.csv","../scores/training/cs/shannon_wg_ROC.csv","../scores/training/cs/vonNeumann_ROC.csv","../scores/training/cs/sop_ROC.csv","../scores/training/cs/sop_wg_ROC.csv","../scores/training/cs/windows1_ROC.csv","../scores/training/cs/windows2_ROC.csv"]
testfile = "../scores/test/TM_conserv_scores_with_window.csv"
indices = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
directions = [1,0,1,0,1,0,1,0,1,1,1,0,0,1,1]

# physco-features
infiles = ["../scores/training/physco/1.csv","../scores/training/physco/2.csv","../scores/training/physco/3.csv","../scores/training/physco/4.csv","../scores/training/physco/5.csv","../scores/training/physco/6.csv","../scores/training/physco/7.csv","../scores/training/physco/8.csv","../scores/training/physco/9.csv","../scores/training/physco/10.csv","../scores/training/physco/11.csv","../scores/training/physco/12.csv","../scores/training/physco/13.csv"]
outfiles = ["../scores/training/physco/1_ROC.csv","../scores/training/physco/2_ROC.csv","../scores/training/physco/3_ROC.csv","../scores/training/physco/4_ROC.csv","../scores/training/physco/5_ROC.csv","../scores/training/physco/6_ROC.csv","../scores/training/physco/7_ROC.csv","../scores/training/physco/8_ROC.csv","../scores/training/physco/9_ROC.csv","../scores/training/physco/10_ROC.csv","../scores/training/physco/11_ROC.csv","../scores/training/physco/12_ROC.csv","../scores/training/physco/13_ROC.csv"]
testfile = "../scores/test/TM_physco.csv"
indices = [5,6,7,8,9,10,11,12,13,14,15,16,17]
directions = [1,1,1,1,1,1,1,1,1,1,1,1,1]
'''
# physco-features - 2 (from the web site exPASY)
infiles = ["../PhyscoScores/thresholding/5.csv","../PhyscoScores/thresholding/6.csv","../PhyscoScores/thresholding/7.csv","../PhyscoScores/thresholding/8.csv","../PhyscoScores/thresholding/9.csv","../PhyscoScores/thresholding/10.csv","../PhyscoScores/thresholding/11.csv","../PhyscoScores/thresholding/12.csv","../PhyscoScores/thresholding/13.csv","../PhyscoScores/thresholding/14.csv","../PhyscoScores/thresholding/15.csv","../PhyscoScores/thresholding/16.csv","../PhyscoScores/thresholding/17.csv","../PhyscoScores/thresholding/18.csv","../PhyscoScores/thresholding/19.csv","../PhyscoScores/thresholding/20.csv","../PhyscoScores/thresholding/21.csv","../PhyscoScores/thresholding/22.csv","../PhyscoScores/thresholding/23.csv","../PhyscoScores/thresholding/24.csv","../PhyscoScores/thresholding/25.csv"]
#outfiles = ["../scores/training/physco2/1_ROC.csv","../scores/training/physco2/2_ROC.csv","../scores/training/physco2/3_ROC.csv","../scores/training/physco2/4_ROC.csv","../scores/training/physco2/5_ROC.csv","../scores/training/physco2/6_ROC.csv","../scores/training/physco2/7_ROC.csv","../scores/training/physco2/8_ROC.csv","../scores/training/physco2/9_ROC.csv","../scores/training/physco2/10_ROC.csv","../scores/training/physco2/11_ROC.csv","../scores/training/physco2/12_ROC.csv","../scores/training/physco2/13_ROC.csv","../scores/training/physco2/14_ROC.csv","../scores/training/physco2/15_ROC.csv","../scores/training/physco2/16_ROC.csv","../scores/training/physco2/17_ROC.csv"]
testfile = "../PhyscoScores/TM_physco.csv"
indices = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
directions = [1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1]
# -----------------------------------------------------------------------------
import math

for i in range(len(infiles)):
    f1 = infiles[i]
    #f2 = outfiles[i]
    f3 = testfile
    d = indices[i]
    direction = directions[i]

    #print("\n*************************************")
    #print(f1.strip().split("/")[-1].split(".")[0])

    # best finds
    bJ = 0
    bS = 0
    btp = 0
    btn = 0
    bfp = 0
    bfn = 0

    # does all the computation magic
    with open(f1) as fps:#, open(f2,"w") as ofp:
        for line in fps:

            parts = line.strip().split(",")
            if parts[0] == 'NaNNaN':
                continue
            score = float(parts[0])
            tp = int(parts[1])
            tn = int(parts[2])
            fp = int(parts[3])
            fn = int(parts[4])

            if direction == 0:
                t = fp
                fp = tn
                tn = t
                t = fn
                fn = tp
                tp = t

            sensitivity = tp*1.0/(tp+fn)
            specificity = tn*1.0/(tn+fp)
            J = sensitivity + specificity - 1
            #J = (tp+tn)/(tp+tn+fp+fn); 

            if J > bJ:
                bJ = J
                bS = score
                btp = tp
                btn = tn
                bfp = fp
                bfn = fn

        #ofp.write(",".join([str(s) for s in [score, sensitivity, 1.0-specificity]])+"\n")

    # print the best youden results and the metrics
    try:
        accuracy = (btp+btn)*1.0/(btp+btn+bfp+bfn)
        sensitivity = btp*1.0/(btp+bfn)
        specificity = btn*1.0/(btn+bfp)
        precision = btp*1.0/(btp+bfp)
        npv = btn*1.0/(btn+bfn)
        mcc = (btp*btn - bfp*bfn)*1.0/math.sqrt((btp+bfp)*(btp+bfn)*(btn+bfp)*(btn+bfn))
    except ZeroDivisionError as e:
        print([btp,btn,bfp,bfn])
        print(bS)
        print(bJ)
        print(line)
        #raise(e)

    #print("%f" %bS)
    #print("Best threshold: %f" %bS)
    #print("Confusion matrix values: %d, %d, %d, %d" %(btp,btn,bfp,bfn))
    #print("Best Youden Index: %f" %bJ)
    #print("Values of 6 metric values on training dataset are as follows:")
    print(i+1)
    print("%f,%d,%f,%f,%f,%f,%f,%f" %(bS,direction,accuracy,sensitivity,specificity,precision,npv,mcc))

    # -----------------------------------------------------------------------------
    # test out the Benchmark dataset
    tp = 0
    tn = 0
    fp = 0
    fn = 0

    with open(f3) as fps: 
        for line in fps:
            parts = line.strip().split(",")

            y = int(parts[4])
            if 'NaN' in line:
                continue
            score = float(parts[int(d)])

            if direction == 1:
                if score > bS:
                    if y == 1:
                        tp += 1
                    else:
                        fp += 1
                else:
                    if y == 1:
                        fn += 1
                    else:
                        tn += 1
            else:
                if score <= bS: 
                    if y == 1:
                        tp += 1
                    else:
                        fp += 1
                else:
                    if y == 1:
                        fn += 1
                    else:
                        tn += 1

    accuracy = (tp+tn)*1.0/(tp+tn+fp+fn)
    sensitivity = tp*1.0/(tp+fn)
    specificity = tn*1.0/(tn+fp)
    precision = tp*1.0/(tp+fp)
    npv = tn*1.0/(tn+fn)
    mcc = (tp*tn - fp*fn)*1.0/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

    #print("Values of 6 evaluation metrics on Benchmark datasets are as follows:")
    print("%f,%f,%f,%f,%f,%f" %(accuracy,sensitivity,specificity,precision,npv,mcc))
