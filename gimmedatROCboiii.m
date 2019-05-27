%%%%%%%%%%%%%%%%%%%
%  Plot ROC curve and compute AUC for 4 different classification algorithms using fNIRS data to predict stimulation signal:
%		GLM with logistic regression
%		SVM
%		Classification trees
%%%%%%%%%%%%%%%%%%%

%% Load data and prepare arrays
load('/Users/e/Desktop/RDM/09_Data_after_cleaning/S01_pilot_alternateSpeaker/S01_pilot_stimuliNIRSproc.mat');

responses = squeeze(nirs{1,1}(:,:,1));
predictors = squeeze(stim{1,1}(:,:,1));

resp         = rangenorm(responses(1,:));
pred         = rangenorm(predictors(1,:));

respSVM   = rangenorm(responses(:,:));

%% Fit a logistic regression model to estimate the posterior probabilities for FNIRS signal to match STIMULUS signal 
mdl         = fitglm(pred,resp,'Distribution','binomial','Link','logit');
score_log = mdl.Fitted.Probability; % Probability estimates

% Compute the standard ROC curve using the probabilities for scores.
[Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,1);

%% Train an SVM classifier on the same sample data. Standardize the data.
mdlSVM = fitcecoc(pred',resp);

% Compute the posterior probabilities (scores).
[~,score_svm] = resubPredict(mdlSVM);
% The second column of score_svm contains the posterior probabilities

% Compute the standard ROC curve using the scores from the SVM model.
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_log,1);

% Plot the ROC curves on the same graph.
figure
plot(Xlog,Ylog, 'r')
hold on
plot(Xsvm,Ysvm, 'b')
legend('Logistic Regression','Support Vector Machines','Location','Best')
xlabel('False positive rate'); ylabel('True positive rate');
title('ROC Curves for Logistic Regression and SVM')
hold on
% Although SVM produces better ROC values for higher thresholds, logistic
% regression is usually better at distinguishing the FNIRS from STIMULUS

% Compare the area under the curve for both classifiers.
AUClog
AUCsvm

%% Compute classification trees
load fisheriris
%The column vector, species, consists of iris flowers of three different species: setosa, versicolor, virginica. 
%The double matrix meas consists of four types of % measurements on the flowers: sepal length, sepal width, petal length, and petal width. All measures are in centimeters.
% Train a classification tree using the sepal length and width as the predictor variables. It is a good practice to specify the class names.

Model = fitctree(meas(:,1:2),species, ...
    'ClassNames',{'setosa','versicolor','virginica'});

% Predict the class labels and scores for the species based on the tree Model.
[~,score] = resubPredict(Model);
%The scores are the posterior probabilities that an observation (a row in the data matrix) belongs to a class. The columns of score correspond to the classes specified % by 'ClassNames'. 
%So, the first column corresponds to setosa, the second corresponds to versicolor, and the third column corresponds to virginica.

% Compute the ROC curve for the predictions that an observation belongs to versicolor, given the true class labels species. 
%Also compute the optimal operating point and % y values for negative subclasses. Return the names of the negative classes.
diffscore = score(:,2) - max(score(:,1),score(:,3));
[X,Y,T,~,OPTROCPT,suby,subnames] = perfcurve(species,diffscore,'versicolor');
%X, by default, is the false positive rate (fallout or 1-specificity) and Y, by default, is the true positive rate (recall or sensitivity). 
%The positive class label is %versicolor. Because a negative class is not defined, perfcurve assumes that the observations that do not belong to the positive class are in one class. The function % %accepts it as the negative class.
OPTROCPT
suby
subnames
%OPTROCPT =
%
%    0.1000    0.8000
%subnames =
%
%  1x2 cell array
%
%    {'setosa'}    {'virginica'}

% Plot the ROC curve and the optimal operating point on the ROC curve.
figure
plot(X,Y)
hold on
plot(OPTROCPT(1),OPTROCPT(2),'ro')
xlabel('False positive rate')
ylabel('True positive rate')
title('ROC Curve for Classification by Classification Trees')
hold on

% Find the threshold that corresponds to the optimal operating point.
T((X==OPTROCPT(1))&(Y==OPTROCPT(2)))

%Specify virginica as the negative class and compute and plot the ROC curve for versicolor.
%Again, you must supply perfcurve with a function that factors in the scores of the negative class. An example of a function to use is .
diffscore = score(:,2) - score(:,3);
[X,Y,~,~,OPTROCPT] = perfcurve(species,diffscore,'versicolor', ...
    'negClass','virginica');
OPTROCPT
figure, plot(X,Y)
hold on
plot(OPTROCPT(1),OPTROCPT(2),'ro')
xlabel('False positive rate')
ylabel('True positive rate')
title('ROC Curve for Classification by Classification Trees')
hold on

%% Rangenorm
function out = rangenorm(data)

out = (data-min(data))/(max(data)-min(data));
end
