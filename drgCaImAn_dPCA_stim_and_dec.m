function handles_out=drgCaImAn_dPCA_stim_and_dec(handles_choices)
%This program performs a demixed principal component analysis (dPCA) 
%according to Kobak et al eLife 2016 DOI: 10.7554/eLife.10989
%
% the input is either no input in which case you will be asked to enter a pre_per file 
% or a handles_choices variable with the pre_per file location and all
% other choices
%
% To find where you setup the dPCA search for "dPCA_input=="

close all

if exist('handles_choices')==0
    clear all

    %Load file
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');
    do_dPCA=1; %If do_dPCA the program performs dPCA
    show_figures=1; %Show the figures

    dPCA_input=1; %If dPCA_input==1 the program uses the data in the input file
    %if dPCA input==0 the program uses simulated Romo data

else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
    show_figures=handles_choices.show_figures;
    do_dPCA=handles_choices.do_dPCA;
    dPCA_input=handles_choices.dPCA_input;

end


post_shift=0;

trial_time_from=-10; %-10
trial_time_to=20; %30




%     MLalgo_to_use=6;
p_threshold=1.1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test

convert_z=2; %0, do not convert, 1= z=x/std, 2= per trial z=(x-x_pre)/std
pre_time=[-5 -1]; %Used to calcuate x_pre



ifSimultaneousRecording=false;

moving_mean_n=20;
% show_figures=1;
 
load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
fprintf(1, 'p threshold = %d\n',p_threshold);

% 
% classifier_names{1}='Linear Discriminant';
% classifier_names{2}='Support Vector Machine';
% classifier_names{3}='Naive Bayes Classifier';
% classifier_names{4}='Neural Network';
% classifier_names{5}='Decision tree';
% classifier_names{6}='Binomial glm';

handles_out.trial_time_from=trial_time_from;
handles_out.trial_time_to=trial_time_to;
handles_out.pre_perFileName=pre_perFileName;
handles_out.pre_perPathName=pre_perPathName;
handles_out.post_shift=post_shift;
handles_out.pre_time=pre_time;
% handles_out.MLalgo_to_use=MLalgo_to_use;

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

hit_per_trial=[];
miss_per_trial=[];
cr_per_trial=[];
fa_per_trial=[];

if convert_z==2
    std_traces=zeros(1,size(traces,1));
    for trace_no=1:size(traces,1)
        std_traces(trace_no)=std(traces(trace_no,:));
    end
    std_traces=std_traces';
end

%Post points
Nall=size(traces,1);
dt=time(2)-time(1);
ii_pre=ceil(pre_time/dt);
ii_p_threshold=ceil(dt_p_threshold/dt);
no_points_before=floor(-trial_time_from/dt);
no_points_after=floor(trial_time_to/dt);
measurements_per_trial=[];
sp_or_sm=[];
mouse_licked=[];

%Do S+
at_end=0;
this_ii=0;
ii_post=0;

ii=0;
trNo=0;
ii_sp_post=0;
ii_sm_post=0;

dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI=[];
dFFs_sm_per_trial_per_ROI=[];

while (at_end==0)
    next_ii=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    if ~isempty(next_ii)
        %This is S+
        if epochs(this_ii+1+next_ii)==6
            %Hit
            hit_per_trial=[hit_per_trial 1];
            miss_per_trial=[miss_per_trial 0];
            mouse_licked=[mouse_licked 1];
        else
            %Miss
            hit_per_trial=[hit_per_trial 0];
            miss_per_trial=[miss_per_trial 1];
            mouse_licked=[mouse_licked 0];
        end
        cr_per_trial=[cr_per_trial 0];
        fa_per_trial=[fa_per_trial 0];

        if (no_points_after+this_ii+next_ii<length(epochs))&(this_ii+next_ii-no_points_before>0)
            trNo=trNo+1;
             if convert_z==2
                mean_trace_before=mean(traces(:,this_ii+next_ii+ii_pre(1):this_ii+next_ii+ii_pre(2)),2);
                mean_trace_before=repmat(mean_trace_before,1,size(traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1),2));
                these_std_traces=repmat(std_traces,1,size(traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1),2));
                measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=((traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1)...
                    -mean_trace_before)./these_std_traces)';
            else
                measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1)';
            end
%             measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_before+no_points_after-1)';
            ii_sp_post=ii_sp_post+1;
            dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,this_ii+next_ii:this_ii+next_ii+ii_p_threshold);        
            this_ii=this_ii+next_ii+no_points_after;
            sp_or_sm(trNo)=1;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

if convert_z==1
    for trace_no=1:size(traces,1)
        traces(trace_no,:)=traces(trace_no,:)/std(traces(trace_no,:));
    end
end

 

%Do S-
at_end=0;
this_ii=0;
ii=0;


while (at_end==0)
    next_ii=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    if ~isempty(next_ii)
        %This is S-
        if epochs(this_ii+1+next_ii)==9
            %CR
            cr_per_trial=[cr_per_trial 1];
            fa_per_trial=[fa_per_trial 0];
            mouse_licked=[mouse_licked 0];
        else
            %FA
            cr_per_trial=[cr_per_trial 0];
            fa_per_trial=[fa_per_trial 1];
            mouse_licked=[mouse_licked 1];
        end
        hit_per_trial=[hit_per_trial 0];
        miss_per_trial=[miss_per_trial 0];
        if (no_points_after+this_ii+next_ii<length(epochs))&(this_ii+next_ii-no_points_before>0)
            trNo=trNo+1;
            if convert_z==2
                mean_trace_before=mean(traces(:,this_ii+next_ii+ii_pre(1):this_ii+next_ii+ii_pre(2)),2);
                mean_trace_before=repmat(mean_trace_before,1,size(traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1),2));
                these_std_traces=repmat(std_traces,1,size(traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1),2));
                measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=((traces(:,this_ii+next_ii-no_points_before:this_ii+next_ii+no_points_after-1)...
                    -mean_trace_before)./these_std_traces)';
            else
                measurements_per_trial(1:no_points_before+no_points_after,:,trNo)=traces(:,this_ii+next_ii-no_points_before:this_ii+next_iie+no_points_after-1)';
            end
            ii_sm_post=ii_sm_post+1;
            dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,this_ii+next_ii:this_ii+next_ii+ii_p_threshold);
            this_ii=this_ii+next_ii+no_points_after;
            sp_or_sm(trNo)=0;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

%Calculate percent correct
handles_out.percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);
fprintf(1, ['\nPercent correct ' num2str(handles_out.percent_correct) '\n\n']);

p_values=ones(1,size(dFFs_sm_per_trial_per_ROI,2));
for iiROI=1:size(dFFs_sm_per_trial_per_ROI,2)
    dFF_sm=zeros(size(dFFs_sm_per_trial_per_ROI,1),size(dFFs_sm_per_trial_per_ROI,3));
    dFF_sm(:,:)=dFFs_sm_per_trial_per_ROI(:,iiROI,:);
    dFF_sp=zeros(size(dFFs_sp_per_trial_per_ROI,1),size(dFFs_sp_per_trial_per_ROI,3));
    dFF_sp(:,:)=dFFs_sp_per_trial_per_ROI(:,iiROI,:);
    
    [h,p_values(iiROI)]=ttest2(mean(dFF_sp,2),mean(dFF_sm,2));
end

p_value_mask=logical(p_values<=p_threshold);

%Trim the number of ROIs in all matrices
noROIs_before_trimming=size(measurements_per_trial,2);
% dFF_per_trial_sp=dFF_per_trial_sp(:,p_value_mask,:);
% dFF_per_trial_sm=dFF_per_trial_sm(:,p_value_mask,:);
measurements_per_trial=measurements_per_trial(:,p_value_mask,:);
traces=traces(p_value_mask,:);
no_traces=size(traces,1);


%Plot the trimmed traces
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
    end
    
    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    xlabel('time(sec)')
    title(['dFF timecourses after p value trimming ' num2str(size(measurements_per_trial,2)) ' ROIs'])
end

fprintf(1, ['Demixed PCA run with %d ROIs (original no ROIs %d)...\n'],size(measurements_per_trial,2),noROIs_before_trimming);


handles_out.Nall=Nall;
handles_out.dt=dt;
handles_out.no_points_before=no_points_before;
handles_out.no_points_after=no_points_after;
handles_out.measurements_per_trial=measurements_per_trial;
handles_out.sp_or_sm=sp_or_sm;


%This section transfers the data to the dPCA input format
%
% The data should be
% joined in three arrays of the following sizes (for the Romo-like task):
%
% trialNum: N x S
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
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


numComp=10; %Components to be extracted from data by dpca
lambda=0.00001; %Note: if I do not make this larger than thisLambda = 1e-10; in line 131 of dpca
%I get Error using  /

switch dPCA_input
    case 1
        %Use the data from the file
        N = size(traces,1);%100;   % number of neurons
        T = size(measurements_per_trial,1);% number of time points
        S = 2; % number of stimuli
        D = 2; %Number of decisions


        time=(1:T)*dt+trial_time_from;



        % setting random number of repetitions for each neuron and condition
        % ifSimultaneousRecording = true;  % false for Romo's data
        % change this to simulate simultaneous
        % recordings (they imply the same number
        % of trials for each neuron)

        % trialNum: N x S x D
        trial_number=zeros(S,D);
        for trNum=1:length(sp_or_sm)
            s=sp_or_sm(trNum)+1;
            d=mouse_licked(trNum)+1;
            trial_number(s,d)=trial_number(s,d)+1;
        end


        E=max(trial_number(:)); %Maximal number of trial repetitions

        firingRates=NaN(N,S,D,T,E);

        %Tranfer the data
        % firingRates: N x S x D x T x maxTrialNum
        % trialNum: N x S x D
        trialNum=zeros(N,S,D);
        for trNum=1:length(sp_or_sm)
            s=sp_or_sm(trNum)+1;
            d=mouse_licked(trNum)+1;
            for n=1:N
                trialNum(n,s,d)=trialNum(n,s,d)+1;
                firingRates(n,s,d,:,trialNum(n,s,d))=measurements_per_trial(:,n,trNum);
            end
        end

    case 2
        %Simulated Romo data

        N = 100;%100;   % number of neurons
        T = 20;%20;     % number of time points
        S = 7;%7;       % number of stimuli
        D = 2;          % number of decisions
        E = 20;%20;     % maximal number of trial repetitions

        % generating three components
        time = (1:T) / 10;
        component{1} = bsxfun(@times, ones(1,S,D,T,E), shiftdim(sin(time), -2)) * 20;
        component{2} = bsxfun(@times, ones(1,S,D,T,E), (1:S)-ceil(S/2)) .* bsxfun(@times, ones(1,S,D,T,E), shiftdim(time, -2));
        component{3} = bsxfun(@times, ones(1,S,D,T,E), shiftdim([-1 1],-1)) .* ...
            bsxfun(@times, ones(1,S,D,T,E), shiftdim([time(1:end/2)*0 time(1:end/2)], -2)) * 6;

        % mixing them randomly along non-orthogonal axes
        V = randn(N, length(component));
        V(:,2) = V(:,2) + V(:,1)/3;
        V(:,3) = V(:,3) + V(:,2)/4;
        V = bsxfun(@times, V, 1./sqrt(sum(V.^2)));
        for c = 1:length(component)
            components(:,c) = component{c}(:,:);
        end
        firingRates = V*components';
        clear component components

        % adding noise and normalizing
        firingRates = firingRates + randn(size(firingRates));
        firingRates = firingRates - min(firingRates(:)) + 10;
        firingRates = reshape(firingRates, [N S D T E]);

        % setting random number of repetitions for each neuron and condition
        ifSimultaneousRecording = true;  % false for Romo's data
        % change this to simulate simultaneous
        % recordings (they imply the same number
        % of trials for each neuron)
        if ~ifSimultaneousRecording
            trialNum = randi([5 E], [N S D]);
        else
            trialNum = repmat(randi([5 E], [1 S D]), [N 1 1]);
        end

        for n = 1:N
            for s = 1:S
                for d = 1:D
                    firingRates(n,s,d,:,trialNum(n,s,d)+1:E) = nan;
                end
            end
        end

    case 3
        %Simulated go-no go data
        %The number of ROIs and trials is the same as in the actual data

        %Use the data from the pre_per file
        N = size(traces,1);%100;   % number of neurons
        T = size(measurements_per_trial,1);% number of time points
        S = 2; % number of stimuli
        D = 2; %Number of decisions

        perCorr=handles_out.percent_correct;

        time=(1:T)*dt+trial_time_from;
        odor_response=1;
        decision_response=0;
        SD_odor=0.5;


        % setting random number of repetitions for each neuron and condition
        % ifSimultaneousRecording = true;  % false for Romo's data
        % change this to simulate simultaneous
        % recordings (they imply the same number
        % of trials for each neuron)

        % trialNum: N x S x D
        trial_number=zeros(S,D);
        for trNum=1:length(sp_or_sm)
            s=sp_or_sm(trNum)+1;
            d=mouse_licked(trNum)+1;
            trial_number(s,d)=trial_number(s,d)+1;
        end


        E=max(trial_number(:)); %Maximal number of trial repetitions

        firingRates=NaN(N,S,D,T,E);

        %Different numbers of neurons are engaged as stimulus neurons or
        %decision neurons

        %Learning
%         fraction_of_neurons=[0.01 0.05 0.2];

        %No learning
        fraction_of_neurons=[0.2 0.2 0.2];

        NI=ceil(0.2*N);

        %Note: make sure NS+ND+NI<1

        if perCorr<65%
            NS=ceil(fraction_of_neurons(1)*N);
            ND=ceil(fraction_of_neurons(1)*N);
        else
            if perCorr>=80
                NS=ceil(fraction_of_neurons(3)*N);
                ND=ceil(fraction_of_neurons(3)*N);
            else
                NS=ceil(fraction_of_neurons(2)*N);
                ND=ceil(fraction_of_neurons(2)*N);
            end
        end


        % Simulate the data
        % firingRates: N x S x D x T x maxTrialNum
        % trialNum: N x S x D
        trialNum=zeros(N,S,D);
        for trNum=1:length(sp_or_sm)
            s=sp_or_sm(trNum)+1;
            d=mouse_licked(trNum)+1;
            for n=1:N
                trialNum(n,s,d)=trialNum(n,s,d)+1;
                %All cells have some noise
                firingRates(n,s,d,:,trialNum(n,s,d))=SD_odor*randn(1,T);

                %Stimulus sensitive neurons
                if  n<=NS
                    %These are stimulus sensitive neurons
                    if s==1
                        %S+
                        firingRates(n,s,d,(time>0)&(time<=5),trialNum(n,s,d))=firingRates(n,s,d,(time>0)&(time<=5),trialNum(n,s,d))+odor_response;
                    end
                end

                %Decision neurons
                if (n>NS)&(n<NS+ND)
                    if d==1
                        %Mouse licked
                        firingRates(n,s,d,(time>1)&(time<=6),trialNum(n,s,d))=firingRates(n,s,d,(time>1)&(time<=6),trialNum(n,s,d))+decision_response;
                    end
                end


                %Odor and decision independent neurons
                if (n>NS+ND)&(n<NS+ND+NI)
                    firingRates(n,s,d,(time>-2)&(time<=20),trialNum(n,s,d))=firingRates(n,s,d,(time>-2)&(time<=20),trialNum(n,s,d))+odor_response;
                end

            end
        end
end

%Save the dPCA inputs in the output variable
handles_out.N = N;
handles_out.T = T;
handles_out.S = S;
handles_out.D = D;
handles_out.time = time;
handles_out.E=E;
handles_out.firingRates=firingRates;
handles_out.trialNum=trialNum;

% computing PSTHs
% firingRatesAverage: N x S x D x T
firingRatesAverage = nanmean(firingRates, 5);
firingRatesAverage(isnan(firingRatesAverage))=0; %This is key when there are zero trials for some s,d combinations.
handles_out.firingRatesAverage=firingRatesAverage;

if do_dPCA==1
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
 
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

    % Time events of interest (e.g. stimulus onset/offset, cues etc.)
    % They are marked on the plots with vertical lines
    timeEvents = 0;

    % check consistency between trialNum and firingRates
    for n = 1:size(firingRates,1)
        for s = 1:size(firingRates,2)
            for d = 1:size(firingRates,3)
                assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
            end
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

    % This is the core function.
    % W is the decoder, V is the encoder (ordered by explained variance),
    % whichMarg is an array that tells you which component comes from which
    % marginalization
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
        'legendSubplot', 16,          ...
        'displaySumStimComps', 0);

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

    %% Plot the different components and 3D for each trial

    % drg_dpca_plot_per_trial(firingRatesAverage, firingRates, W, V, @dpca_plot_default, ...
    %     'explainedVar', explVar, ...
    %     'marginalizationNames', margNames, ...
    %     'marginalizationColours', margColours, ...
    %     'whichMarg', whichMarg,                 ...
    %     'time', time,                        ...
    %     'timeEvents', timeEvents,               ...
    %     'timeMarginalization', 3,           ...
    %     'legendSubplot', 16,          ...
    %     'displaySumStimComps', 0);
end

pffft=1;