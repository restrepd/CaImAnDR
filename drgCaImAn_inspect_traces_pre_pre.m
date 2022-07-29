function handles_out=drgCaImAn_inspect_traces_pre_pre(handles_choices)
%This program trains several decoding algorithms with the post odorant and then determines what happens throughout the entire timecouse
%The user enters the choices entered under exist('handles_choices')==0
%
% processing_algorithm= 1 and 2 were used for troublehsooting and do not
% produce reliable results because of overtraining, use
% processing_algoritm=3, that was vetted for our manuscript
%
%
% the input is a pre_per file version 2
if exist('handles_choices')==0
    clear all
    close all

    %Load file
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');

    processing_algorithm=3; %Use 3
    k_fold=5; %Only used for processing_algorithm=2,
    post_time=5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on
    post_shift=0; %Set to 0 if you want to train with odor on points
    pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
    MLalgo_to_use=[1]; %Vector with the decoding algorithms you want to use
    ii_cost=3;
    p_threshold=1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
    dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures

    handles_choices.pre_perFileName=pre_perFileName;
    handles_choices.pre_perPathName=pre_perPathName;
    handles_choices.processing_algorithm=processing_algorithm;
    handles_choices.post_time=post_time;
    handles_choices.k_fold=k_fold;
    handles_choices.post_shift=post_shift;
    handles_choices.MLalgo_to_use=MLalgo_to_use;
    handles_choices.pre_time=pre_time;
    handles_choices.p_threshold=p_threshold;
    handles_choices.dt_p_threshold=dt_p_threshold;
    handles_choices.show_figures=show_figures;
    handles_choices.ii_cost=ii_cost;

else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
    processing_algorithm=handles_choices.processing_algorithm;
    post_time=handles_choices.post_time;
    k_fold=handles_choices.k_fold;
    post_shift=handles_choices.post_shift;
    MLalgo_to_use=handles_choices.MLalgo_to_use;
    pre_time=handles_choices.pre_time;
    p_threshold=handles_choices.p_threshold;
    dt_p_threshold=handles_choices.dt_p_threshold;
    show_figures=handles_choices.show_figures;
    ii_cost=handles_choices.ii_cost;
end

warning('off')

%Restart random seeds
rng('shuffle');

convert_z=1; %Convert dFF traces to z
dt_span=40; %Seconds shown before and after odor on in the figures
moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
no_shuffles=10; %Number of shuffles for per trial shuffling

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
fprintf(1, 'post_time = %d, p_threshold= %d, post_shift= %d, cost %d\n',post_time,p_threshold,post_shift, ii_cost);

if convert_z==1
    for trace_no=1:size(traces,1)
        traces(trace_no,:)=traces(trace_no,:)/std(traces(trace_no,:));
    end
end

classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';
classifier_names{6}='Binomial glm';

handles_out.handles_choices=handles_choices;
handles_out.pre_perFileName=pre_perFileName;
handles_out.pre_perPathName=pre_perPathName;
handles_out.post_time=post_time;
handles_out.k_fold=k_fold;
handles_out.post_shift=post_shift;
handles_out.pre_time=pre_time;
handles_not_out.MLalgo_to_use=MLalgo_to_use;

figNo=0;

%time has the time for the dF/F traces(ROI,time)

%If for_grant is on this generates the figure for the grant
for_grant=0;
if for_grant==1
    to_traces=17;
else
    to_traces=no_traces;
end




figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on

% Determine the y spacing of the traces
y_shift=6*(prctile(traces(:),95)-prctile(traces(:),5));

%Plot the traces and do z normalization
%For S+ and S- plot odor on and reinforcement

for epoch=1:handles.dropcData.epochIndex
    %Epoch 2 is odor on, 3 is odor off
    plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
    if plot_epoch
        if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                '-r','LineWidth',1.5)
        else
            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                '-b','LineWidth',1.5)
        end
    end
end

for trNo=1:to_traces
    % for trNo=1:20
    plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1.5)
end

if for_grant==1
    ylim([5 145])
    xlim([0 800])
else
    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
end

xlabel('time(sec)')
title(['All dFF timecourses ' num2str(size(traces,1)) ' ROIs'])


%Calculate the crosscorrelations
croscorr_traces=abs(corrcoef(traces')); %please note that I am using the absolute value

%Set autocorrelations to zero
for ii=1:size(croscorr_traces,1)
    croscorr_traces(ii,ii)=0;
end
Z = linkage(croscorr_traces,'complete','correlation');
figNo=figNo+1;
[H,T,outperm]=dendrogram(Z,0,'Orientation','left');
set(H,'LineWidth',2)
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.05 .1 .1 .8])

%re-sort the matrix
for ii=1:size(croscorr_traces,1)
    for jj=1:size(croscorr_traces,1)
        perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
    end
end

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .8 .8])
hold on
pcolor(perm_croscorr_traces)
colormap fire

caxis([0    0.6])
title(['Cross correlations for ' pre_perFileName])

%Plot rainbow
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


prain=[0:0.6/99:0.6];
pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%             colormap jet
colormap fire
shading interp
ax=gca;
set(ax,'XTickLabel','')

%Plot sorted traces
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
hold on



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
    plot(time,traces(outperm(trNo),:)+y_shift*trNo,'-k','LineWidth',1)
end

ylim([-y_shift*0.2 (no_traces+2)*y_shift])
xlabel('time(sec)')
title(['All dFF timecourses hierarchical order ' num2str(size(traces,1)) ' ROIs'])


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
ii_p_threshold=ceil(dt_p_threshold/dt);
no_points_post=floor(post_time/dt);
no_points_post_shift=floor(post_shift/dt);
no_points_pre=floor(pre_time/dt);
measurements_post=[];
measurements_pre=[];
epochs_sp_post=zeros(1,length(time));
epochs_sp_pre=zeros(1,length(time));
epochs_sm_post=zeros(1,length(time));
epochs_sm_pre=zeros(1,length(time));
which_model_for_traces_loo=no_odor_trials*ones(1,size(traces,2));

%Do both S+ and S-
at_end=0;
this_ii=0;
ii_post=0;
ii_pre=0;
ii=0;
trial_no=0;
ii_sp_post=0;
ii_sm_post=0;
ii_which_model=0;
dt_post_which_model=floor(20/dt); %Points that model will be used beyond the training period
ii_span=ceil(dt_span/dt);
dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI=[];
dFFs_sm_per_trial_per_ROI=[];
hit_per_trial=[];
miss_per_trial=[];
cr_per_trial=[];
fa_per_trial=[];

%training_decisions is 1 if S+ and 2 if S-
while (at_end==0)
    %6 is hit, 7 is miss, 8 FA and 9 CR
    next_ii_sp=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    next_ii_sm=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');

    if (isempty(next_ii_sp))&(isempty(next_ii_sm))
        at_end=1;
    else

        if isempty(next_ii_sm)
            %This is S+
            if epochs(this_ii+1+next_ii_sp)==6
                hit_per_trial=[hit_per_trial 1];
                miss_per_trial=[miss_per_trial 0];
            else
                hit_per_trial=[hit_per_trial 0];
                miss_per_trial=[miss_per_trial 1];
            end
            cr_per_trial=[cr_per_trial 0];
            fa_per_trial=[fa_per_trial 0];

            next_ii=next_ii_sp;
            if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)&(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))
                measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                training_decisions_post(ii_post+1:ii_post+no_points_post,:)=1;
                ii_sp_post=ii_sp_post+1;
                dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;
                trial_no=trial_no+1;
                decisions_per_trial(trial_no)=1;
                which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                this_ii=this_ii+next_ii+no_points_post;
                ii_post=ii_post+no_points_post;
                ii_pre=ii_pre+no_points_pre;
            else
                if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||(no_points_post_shift+this_ii+next_ii-ii_span>0)&(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
                    at_end=1;
                else
                    this_ii=this_ii+next_ii+no_points_post;
                end
            end
        end

        if isempty(next_ii_sp)
            %This is S-
            if epochs(this_ii+1+next_ii_sm)==9
                cr_per_trial=[cr_per_trial 1];
                fa_per_trial=[fa_per_trial 0];
            else
                cr_per_trial=[cr_per_trial 0];
                fa_per_trial=[fa_per_trial 1];
            end
            hit_per_trial=[hit_per_trial 0];
            miss_per_trial=[miss_per_trial 0];

            next_ii=next_ii_sm;
            if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)...
                    &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                training_decisions_post(ii_post+1:ii_post+no_points_post,:)=0;
                ii_sm_post=ii_sm_post+1;
                dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;
                trial_no=trial_no+1;
                decisions_per_trial(trial_no)=0;
                which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                this_ii=this_ii+next_ii+no_points_post;
                ii_post=ii_post+no_points_post;
                ii_pre=ii_pre+no_points_pre;
            else
                if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||(no_points_post_shift+this_ii+next_ii-ii_span>0)&(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
                    at_end=1;
                else
                    this_ii=this_ii+next_ii+no_points_post;
                end
            end
        end

        if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
            if next_ii_sm<next_ii_sp
                %This is S-
                next_ii=next_ii_sm;

                if epochs(this_ii+1+next_ii_sm)==9
                    cr_per_trial=[cr_per_trial 1];
                    fa_per_trial=[fa_per_trial 0];
                else
                    cr_per_trial=[cr_per_trial 0];
                    fa_per_trial=[fa_per_trial 1];
                end
                hit_per_trial=[hit_per_trial 0];
                miss_per_trial=[miss_per_trial 0];

                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post,:)=0;
                    ii_sm_post=ii_sm_post+1;
                    dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                    dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;
                    trial_no=trial_no+1;
                    decisions_per_trial(trial_no)=0;
                    which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                    ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                    ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                    epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                    measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                    epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                    this_ii=this_ii+next_ii+no_points_post;
                    ii_post=ii_post+no_points_post;
                    ii_pre=ii_pre+no_points_pre;
                else
                    if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||(no_points_post_shift+this_ii+next_ii-ii_span>0)&(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
                        at_end=1;
                    else
                        this_ii=this_ii+next_ii+no_points_post;
                    end
                end
            else
                %This is S+
                next_ii=next_ii_sp;

                if epochs(this_ii+1+next_ii_sp)==6
                    hit_per_trial=[hit_per_trial 1];
                    miss_per_trial=[miss_per_trial 0];
                else
                    hit_per_trial=[hit_per_trial 0];
                    miss_per_trial=[miss_per_trial 1];
                end
                cr_per_trial=[cr_per_trial 0];
                fa_per_trial=[fa_per_trial 0];

                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post,:)=1;
                    ii_sp_post=ii_sp_post+1;
                    dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                    dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;
                    trial_no=trial_no+1;
                    decisions_per_trial(trial_no)=1;
                    which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                    ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                    ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                    epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                    measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                    epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                    this_ii=this_ii+next_ii+no_points_post;
                    ii_post=ii_post+no_points_post;
                    ii_pre=ii_pre+no_points_pre;
                else
                    if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||(no_points_post_shift+this_ii+next_ii-ii_span>0)&(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
                        at_end=1;
                    else
                        this_ii=this_ii+next_ii+no_points_post;
                    end
                end
            end
        end


    end

end

%Calculate percent correct
handles_out.percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);

which_model_for_traces_loo(which_model_for_traces_loo>trial_no)=trial_no;
%
% %Do S+
% at_end=0;
% this_ii=0;
% ii_post=0;
% ii_pre=0;
% ii=0;
% trial_no=0;
%
% while (at_end==0)
%     next_ii=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
%     if ~isempty(next_ii)
%         if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
%             measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
%             trial_no=trial_no+1;
%             which_model_for_traces_loo(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=trial_no;
%             ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
%             epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
%             measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
%             epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
%             this_ii=this_ii+next_ii+no_points_post;
%             ii_post=ii_post+no_points_post;
%             ii_pre=ii_pre+no_points_pre;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
%     end
% end
%
% training_decisions_post=ones(1,size(measurements_post,1));
% ii_sp_post=size(measurements_post,1);
%
%
% %Do S-
% at_end=0;
% this_ii=0;
% ii=0;
%
%
% while (at_end==0)
%     next_ii=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
%     if ~isempty(next_ii)
%         if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
%             measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
%             trial_no=trial_no+1;
%             which_model_for_traces_loo(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=trial_no;
%             ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
%             epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
%             measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
%             epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
%             this_ii=this_ii+next_ii+no_points_post;
%             ii_post=ii_post+no_points_post;
%             ii_pre=ii_pre+no_points_pre;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
%     end
% end
%
% ii_sm_post=ii_post-ii_sp_post;
%
% training_decisions_post=[training_decisions_post zeros(1,size(measurements_post,1)-ii_sp_post)];

%Now let's limit the ROIs to those below p_threshold

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
noROIs_before_trimming=size(measurements_post,2);
dFF_per_trial_sp=dFF_per_trial_sp(:,p_value_mask,:);
dFF_per_trial_sm=dFF_per_trial_sm(:,p_value_mask,:);
measurements_post=measurements_post(:,p_value_mask);
measurements_pre=measurements_pre(:,p_value_mask);
traces=traces(p_value_mask,:);
no_traces=size(traces,1);

%Save odor times
handles_out.sp_times=[];
handles_out.sp_times_ii=0;
handles_out.sm_times=[];
handles_out.sm_times_ii=0;

%For S+ and S- plot odor on and reinforcement
for epoch=1:handles.dropcData.epochIndex
    %Epoch 2 is odor on, 3 is odor off
    plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
    if plot_epoch
        if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
            handles_out.sp_times_ii=handles_out.sp_times_ii+1;
            handles_out.sp_times(handles_out.sp_times_ii)=handles.dropcData.epochTime(epoch);
            %             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
            %                 '-r','LineWidth',1)
        else
            handles_out.sm_times_ii=handles_out.sm_times_ii+1;
            handles_out.sm_times(handles_out.sm_times_ii)=handles.dropcData.epochTime(epoch);
            %             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
            %                 '-b','LineWidth',1)
        end
    end
end


%Plot the trimmed traces sorted according to p values
to_sort=zeros(no_traces,2);
to_sort(:,1)=1:no_traces;
to_sort(:,2)=p_values;
sorted=sortrows(to_sort,2);
traces_sorted=zeros(no_traces,1);
traces_sorted(:,1)=sorted(:,1);

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
        plot(time,traces(traces_sorted(trNo),:)+y_shift*trNo,'-k','LineWidth',1)
    end

    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    xlabel('time(sec)')
    title(['dFF timecourses after p value sorting ' num2str(size(measurements_post,2)) ' ROIs'])
end





