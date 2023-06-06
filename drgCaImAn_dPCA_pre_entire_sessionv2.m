function handles_out=drgCaImAn_dPCA_pre_entire_sessionv2(handles_choices)
%This program trains the SVZ with the post odorant and then determines what happens throughout the entire timecouse 

if exist('handles_choices')==0
    clear all
    close all

    %Load file
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');

    processing_algorithm=2;

    %For some reason for post_time=5 we get artificailly high prediction with
    %k_fold=10 even when we post_shift -5 or 20

    %         post_time=5; %The SVZ will be trained with all points 5 sec following odor on
    %         k_fold=10; %Note: make sure that this does not result in points masked=0 --->overtrained
    %         post_shift=0;

    %         post_time=5; %The SVZ will be trained with all points 20-255 sec following odor on
    %         k_fold=10; %Note: make sure that this does not result in points masked=0 --->overtrained
    %         post_shift=20;
    %

    post_time=5; %The SVZ will be trained with all points in 5 sec intervals 0-25 sec following odor on
    k_fold=5; %Note: make sure that this does not result in points masked=0 --->overtrained
    post_shift=20;
    
    trial_time_from=-10;
    trial_time_to=30;

    %     post_time=20; %The SVZ will be trained with all points 20 sec following odor on
    %     k_fold=30; %Ran with 10 in the past, 63 is loo
    %     post_shift=0;


    %     post_time=20; %Because of post_shift=10 the SVZ will be trained with all points 10-30 sec following odor on
    %     k_fold=30; %Ran with 10 in the past, 63 is loo
    %     post_shift=10; %This shifts the post training by this time
    %
    pre_time=5;
    MLalgo_to_use=6;

else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
    processing_algorithm=handles_choices.processing_algorithm;
    post_time=handles_choices.post_time;
    k_fold=handles_choices.k_fold;
    post_shift=handles_choices.post_shift;
    MLalgo_to_use=handles_choices.MLalgo_to_use;
    pre_time=handles_choices.pre_time;
end

 
moving_mean_n=20;
show_figures=1;

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
fprintf(1, 'post_time = %d, k_fold= %d, post_shift= %d\n',post_time,k_fold,post_shift);


classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';
classifier_names{6}='Binomial glm';

handles_out.pre_perFileName=pre_perFileName;
handles_out.pre_perPathName=pre_perPathName;
handles_out.post_time=post_time;
handles_out.k_fold=k_fold;
handles_out.post_shift=post_shift;
handles_out.pre_time=pre_time;
handles_out.MLalgo_to_use=MLalgo_to_use;

figNo=0;

%time has the time for the dF/F traces(ROI,time)
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
    hold on

    % Determine the y spacing of the traces
    y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));

    %Plot the traces and do z normalization
    %For S+ and S- plot odor on and reinforcement
    for epoch=1:handles.dropcData.epochIndex
        %Epoch 2 is odor on, 3 is odor off
        plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
        if plot_epoch
            if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                    '-r','LineWidth',1)
            else
                plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                    '-b','LineWidth',1)
            end
        end
    end
 
    for trNo=1:no_traces
        % for trNo=1:20
        plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        traces(trNo,:)=traces(trNo,:)/std(traces(trNo,:));
    end

    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    xlabel('time(sec)')
end
%epochs is a vector of the length of time that gives information on
%behavior
% 1=Final Valve
% 6=Hit (on for the duration of odor on)
% 7=Miss
% 8=FA
% 9=CR

%For example Hit||Miss shows S+ odor application times (red)
%and FA||CR gives S- (blue)

%Post points
Nall=size(traces,1);
dt=time(2)-time(1);

no_points_before=floor(-trial_time_from/dt);
no_points_after=floor(trial_time_to/dt);
measurements_per_trial=[];
training_decisions=[];
training_decisions_pre=[];

%Do S+
at_end=0;
this_ii=0;
ii_post=0;
ii_pre=0;
ii=0;
trNo=0;


%training_decisions is 1 if S+ and 2 if S-
while (at_end==0)
    next_ii_sp=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    delta_ii_sp=find(epochs(this_ii+next_ii_sp:end)~=epochs(this_ii+next_ii_sp),1,'first');

    next_ii_sm=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    delta_ii_sm=find(epochs(this_ii+next_ii_sm:end)~=epochs(this_ii+next_ii_sm),1,'first');
    
    if (isempty(next_ii_sp))&(isempty(next_ii_sm))
        at_end=1;
    else
        
        if isempty(next_ii_sm)
            next_ii=next_ii_sp;
            if (no_points_after+this_ii+next_ii<length(epochs))&(this_ii+next_ii-no_points_before>0)
                trNo=trNo+1;
                measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_before+no_points_after-1)';
                this_ii=this_ii+next_ii+delta_ii_sp;
                training_decisions(trNo)=1;
                if trNo==1
                    training_decisions_pre(trNo)=1; %This is arbitrary
                else
                    training_decisions_pre(trNo)=training_decisions(trNo-1);
                end
            else
                at_end=1;
            end
        end
        
        if isempty(next_ii_sp)
            next_ii=next_ii_sm;
            if (no_points_after+this_ii+next_ii<length(epochs))&(this_ii+next_ii-no_points_before>0)
                trNo=trNo+1;
                measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_before+no_points_after-1)';
                this_ii=this_ii+next_ii+delta_ii_sm;
                training_decisions(trNo)=2;
                if trNo==1
                    training_decisions_pre(trNo)=1; %This is arbitrary
                else
                    training_decisions_pre(trNo)=training_decisions(trNo-1);
                end
            else
                at_end=1;
            end
        end
        
        if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
            if next_ii_sm<next_ii_sp
                next_ii=next_ii_sm;
                if (no_points_after+this_ii+next_ii<length(epochs))&(this_ii+next_ii-no_points_before>0)
                    trNo=trNo+1;
                    measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_before+no_points_after-1)';
                    this_ii=this_ii+next_ii+delta_ii_sm;
                    training_decisions(trNo)=2;
                    if trNo==1
                        training_decisions_pre(trNo)=1; %This is arbitrary
                    else
                        training_decisions_pre(trNo)=training_decisions(trNo-1);
                    end
                else
                    at_end=1;
                end
            else
                next_ii=next_ii_sp;
                if (no_points_after+this_ii+next_ii<length(epochs))&(this_ii+next_ii-no_points_before>0)
                    trNo=trNo+1;
                    measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_before+no_points_after-1)';
                    this_ii=this_ii+next_ii+delta_ii_sp;
                    training_decisions(trNo)=1;
                    if trNo==1
                        training_decisions_pre(trNo)=1; %This is arbitrary
                    else
                        training_decisions_pre(trNo)=training_decisions(trNo-1);
                    end
                else
                    at_end=1;
                end
            end
        end
        
        
    end
    
end

%Here we shift the training decisions to the decision of the prior trial
training_decisions=training_decisions_pre;

handles_out.Nall=Nall;
handles_out.dt=dt;
handles_out.no_points_before=no_points_before;
handles_out.no_points_after=no_points_after;
handles_out.measurements_per_trial=measurements_per_trial;
handles_out.training_decisions=training_decisions;


%This section transfers the data to the dPCA input format
%
% The data should be
% joined in three arrays of the following sizes (for the Romo-like task):
%
% trialNum: N x S
% firingRates: N x S x T x maxTrialNum
% firingRatesAverage: N x S x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard
% drive as a sparse matrix), then
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

N = size(traces,1);%100;   % number of neurons
T = size(measurements_per_trial,1);%20;     % number of time points
S = 2;%7;       % number of stimuli


time=(1:T)*dt+trial_time_from;

numComp=10; %Components to be extracted from data by dpca
lambda=0.00001; %Note: if I do not make this larger than thisLambda = 1e-10; in line 131 of dpca
%I get Error using  /

% setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % false for Romo's data
% change this to simulate simultaneous
% recordings (they imply the same number
% of trials for each neuron)

E=max([sum(training_decisions==1) sum(training_decisions==2)]);

firingRates=NaN(N,S,T,E);

%Tranfer the data
% firingRates: N x S x T x maxTrialNum
trial_number=zeros(1,2);
for trNum=1:length(training_decisions)
    for s=1:S
        if s==1
            if training_decisions(trNum)==1
                trial_number(1)=trial_number(1)+1;
            end
        else
            if training_decisions(trNum)==2
                trial_number(2)=trial_number(2)+1;
            end
        end
        for n=1:N
            if s==1
                if training_decisions(trNum)==1
                    firingRates(n,s,:,trial_number(1))=measurements_per_trial(:,n,trNum);
                end
            else
                if training_decisions(trNum)==2
                    firingRates(n,s,:,trial_number(2))=measurements_per_trial(:,n,trNum);
                end
            end
            
        end
    end
end


% trialNum: N x S
for n=1:N
    for s=1:S
        if s==0
            trialNum(n,s)=sum(training_decisions==1);
        else
            trialNum(n,s)=sum(training_decisions==2);
        end
    end
end


% computing PSTHs
firingRatesAverage = nanmean(firingRates, 4);



%% Define parameter grouping

% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; here we ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};

combinedParams = {{1, [1 2]}, {2}};
margNames = {'Stimulus', 'Stimulus-independent'};
margColours = [23 100 171; 187 20 25]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = 0;

% check consistency between trialNum and firingRates
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
            assert(isempty(find(isnan(firingRates(n,s,:,1:trialNum(n,s))), 1)), 'Something is wrong!')
    end
end
 
%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);
sgtitle('Regular PCA, minimal plot')

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);
 
% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);

sgtitle('Regular PCA, informative plot')

%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
   'combinedParams', combinedParams);

sgtitle('PCA in each marginalization')
%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ...
    'combinedParams', combinedParams,'lambda',lambda);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);

sgtitle('Demixed PCA')
%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, numComp, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

sgtitle('dPCA with regularization')
%% Optional: estimating "signal variance"

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

sgtitle('dPCA with regularization with signal dFF variance')
% %% Optional: decoding
% 
% decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
% 
% accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'simultaneous', ifSimultaneousRecording, ...
%     'numRep', 5, ...        % increase to 100
%     'filename', 'tmp_classification_accuracy.mat');
% 
% dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
% 
% accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'simultaneous', ifSimultaneousRecording, ...
%     'numRep', 5, ...        % increase to 100
%     'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
%     'filename', 'tmp_classification_accuracy.mat');
% 
% dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
% 
% componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16,                ...
%     'componentsSignif', componentsSignif);

pffft=1;