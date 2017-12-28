% This script does the learning
% with k-fold cross validation
% 
% Method employed does 5-fold 
% cross-validation on the dataset.
%
% Parameter being tuned is the
% number of iterations 
%
% author  : Rashmi Balasubramanyam
%			<rashmi.bmanyam@gmail.com>
% created : 7th Jan, 2017

part1 = csvread('datasets/cs_slng/total_training_cs_slng_1_Matlab.csv');
part2 = csvread('datasets/cs_slng/total_training_cs_slng_2_Matlab.csv');
part3 = csvread('datasets/cs_slng/total_training_cs_slng_3_Matlab.csv');
part4 = csvread('datasets/cs_slng/total_training_cs_slng_4_Matlab.csv');
part5 = csvread('datasets/cs_slng/total_training_cs_slng_5_Matlab.csv');

% number of features
d = size(part1,2)-1;
n = size(part1,1) + size(part2,1) + size(part3,1) + size(part4,1) + size(part5,1);

% choosing from following values of hyper-parameters
T = [10,20,50,70,80,100,200,300,400,500,600,700,800,900,1000,1200,1500,1600,1800,2000]; %[2100:500:5000];%

fileid = fopen('model/cs_slng/errors.txt','w');

trainerror = zeros(6,1);
valerror = zeros(6,1);
for t = T
	t % display the iteration number
	trainkerror = zeros(5,1);
	valkerror = zeros(5,1);
	for k = 1:5 % each fold
		% make training and validation set
		if k == 1
			vy = part1(:,d+1);
			vx = part1(:,1:d);
			ty = [part2(:,d+1); part3(:,d+1); part4(:,d+1); part5(:,d+1)];
			tx = [part2(:,1:d); part3(:,1:d); part4(:,1:d); part5(:,1:d)];
		elseif k == 2
			vy = part2(:,d+1);
			vx = part2(:,1:d);
			ty = [part1(:,d+1); part3(:,d+1); part4(:,d+1); part5(:,d+1)];
			tx = [part1(:,1:d); part3(:,1:d); part4(:,1:d); part5(:,1:d)];
		elseif k == 3
			vy = part3(:,d+1);
			vx = part3(:,1:d);
			ty = [part2(:,d+1); part1(:,d+1); part4(:,d+1); part5(:,d+1)];
			tx = [part2(:,1:d); part1(:,1:d); part4(:,1:d); part5(:,1:d)];
		elseif k == 4
			vy = part4(:,d+1);
			vx = part4(:,1:d);
			ty = [part2(:,d+1); part3(:,d+1); part1(:,d+1); part5(:,d+1)];
			tx = [part2(:,1:d); part3(:,1:d); part1(:,1:d); part5(:,1:d)];
		else
			vy = part5(:,d+1);
			vx = part5(:,1:d);
			ty = [part2(:,d+1); part3(:,d+1); part4(:,d+1); part1(:,d+1)];
			tx = [part2(:,1:d); part3(:,1:d); part4(:,1:d); part1(:,1:d)];
		end
		ty = 2*ty - 1;
		vy = 2*vy - 1;
		% display size of training and validation data
		size(tx)
		size(ty)
		size(vx)
		size(vy)
		% learn the model
		[y, model] = adaboost('train',tx, ty, t);
		clear y
		% compute training error, also display it
		clear predY
		predY = adaboost('apply', tx, model);
		trainkerror(k) = calculate_error(ty, predY)
		% compute validation error also display it
		clear predY
		predY = adaboost('apply', vx, model);
		valkerror(k) = calculate_error(vy, predY)
	end
	trainerror(t) = mean(trainkerror);
	valerror(t) = mean(valkerror);

	% save the errors, can be used for plotting
	fprintf(fileid, 'iteration t = %d\n', t);
	fprintf(fileid, 'mean training error: %12.8f\n',trainerror(t));
	fprintf(fileid, 'mean validation error: %12.8f\n',valerror(t));
end
fprintf(fileid, '---------------------\n\n');

bestT = find(valerror == max(valerror),1);
%bestT = T(bestTindex(1));
fprintf(fileid, 'optimum number of iterations T = %d\n\n', bestT);

% train the final model with best value of T
data = csvread('datasets/cs_slng/total_training_cs_slng_Matlab.csv');

dataX = data(:,1:d);
dataY = data(:,d+1);
dataY = 2*dataY - 1;

[y, finalmodel] = adaboost('train', dataX, dataY, bestT);
clear y
save 'model/cs_slng/adamodel.mat' finalmodel;

clear predY
predY = adaboost('apply', dataX, finalmodel);
terror = calculate_error(dataY, predY);
fprintf(fileid, 'final training error: %12.8f\n', terror);
matrix = confusion_matrix(dataY, predY);
fprintf(fileid, 'tp = %d, tn = %d, fp = %d, fn = %d \n\n', matrix(1), matrix(2), matrix(3), matrix(4));
clearvars predY;

% calculate error on benchmark dataset
benchmark = csvread('datasets/cs_slng/TM2_cs_slng_Matlab.csv');
bX = benchmark(:,1:d);
bY = benchmark(:,d+1);
bY = 2*bY - 1;

clear predY
predY = adaboost('apply', bX, finalmodel);
berror = calculate_error(bY, predY);
fprintf(fileid, 'error on benchmark dataset: %12.8f\n', berror);
matrix = confusion_matrix(bY, predY);
fprintf(fileid, 'tp = %d, tn = %d, fp = %d, fn = %d \n\n', matrix(1), matrix(2), matrix(3), matrix(4));
clearvars predY;
fclose(fileid);
