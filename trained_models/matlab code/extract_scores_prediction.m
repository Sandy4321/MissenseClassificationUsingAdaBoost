% -----------------------------------------------------------------------------
% load the learnt model
load('adamodel.mat')
finalmodel

% -----------------------------------------------------------------------------
% write alphas
fp = fopen('modelparams/alphas.csv','w')
fprintf(fp,'%12.8f\n',finalmodel.alpha)
fclose(fp)

% -----------------------------------------------------------------------------
% save which classifier was chosen in the rounds
fp = fopen('modelparams/classifiers.csv','w')
fprintf(fp,'%d\n',finalmodel.dimension)
fclose(fp)

% -----------------------------------------------------------------------------
% get the predictions and adaboost score for training data
data = csvread('../../../datasets/slng/total_training_slng_scores_Matlab.csv');
d = size(data,2) - 1 % ignore label, it is not a feature
dataX = data(:, 1:d);
dataY = data(:, d+1);
dataY = 2*dataY - 1;

% write predictions
[ypred, confidence] = adaboost('apply', dataX, finalmodel);
ypred = (ypred+1)/2;

% write this to file
fp = fopen('modelparams/training_labels.csv','w')
fprintf(fp, '%d\n', ypred)
fclose(fp)

% write confidence to file
fp = fopen('modelparams/training_confidence.csv','w')
fprintf(fp,'%12.8f\n',confidence)
fclose(fp)

% -----------------------------------------------------------------------------
% get the predictions and adaboost score for benchmark data
trdata = csvread('../../../datasets/slng/TM2_slng_scores_Matlab.csv');
td = size(trdata,2) - 1 % ignore label, it is not a feature
trX = trdata(:, 1:td);
trY = trdata(:, td+1);
trY = 2*trY - 1;

% write predictions
[typred, tconfidence] = adaboost('apply', trX, finalmodel);
typred = (typred+1)/2;

% write this to file
fp = fopen('modelparams/TM2_labels.csv','w')
fprintf(fp, '%d\n', typred)
fclose(fp)

% write confidence to file
fp = fopen('modelparams/TM2_confidence.csv','w')
fprintf(fp,'%12.8f\n',tconfidence)
fclose(fp)

% -----------------------------------------------------------------------------