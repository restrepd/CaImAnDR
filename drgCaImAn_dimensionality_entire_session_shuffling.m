function handles_out=drgCaImAn_dimensionality_entire_session_shuffling(handles_choices)
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
    
%     k_fold=5; %Only used for processing_algorithm=2,
    post_time=30; %Time span
    post_shift=-post_time/2; %This is the time start
%     pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
    p_threshold=1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
%     dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures
    
    handles_choices.pre_perFileName=pre_perFileName;
    handles_choices.pre_perPathName=pre_perPathName;
%     handles_choices.processing_algorithm=processing_algorithm;
    handles_choices.post_time=post_time;
%     handles_choices.k_fold=k_fold;
    handles_choices.post_shift=post_shift;
%     handles_choices.MLalgo_to_use=MLalgo_to_use;
%     handles_choices.pre_time=pre_time;
    handles_choices.p_threshold=p_threshold;
%     handles_choices.dt_p_threshold=dt_p_threshold;
    handles_choices.show_figures=show_figures;
%     handles_choices.ii_cost=ii_cost;
    
else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
   
    post_time=handles_choices.post_time;

    post_shift=handles_choices.post_shift;
  
    p_threshold=handles_choices.p_threshold;
    
    show_figures=handles_choices.show_figures;
    
end



warning('off')

%Restart random seeds
rng('shuffle');
 
convert_z=1; %Convert dFF traces to z
dt_span=40; %Seconds for per trial traces centered on odor on 
moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
no_shuffles=10; %Number of shuffles for per trial shuffling

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_dimensionality_entire_session_shuffling run for ' pre_perFileName '\n']);
fprintf(1, 'post_time = %d, p_threshold= %d, post_shift= %d\n',post_time,p_threshold,post_shift);

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
% handles_out.k_fold=k_fold;
handles_out.post_shift=post_shift;
% handles_out.pre_time=pre_time;
% handles_not_out.MLalgo_to_use=MLalgo_to_use;

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
    title(['All dFF timecourses ' num2str(size(traces,1)) ' ROIs'])
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
% ii_p_threshold=ceil(dt_p_threshold/dt);
no_points_post=floor(post_time/dt);
no_points_post_shift=floor(post_shift/dt);
% no_points_pre=floor(pre_time/dt);
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
ii_span=ceil((dt_span/dt)/2);
dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI=[];
dFFs_sm_per_trial_per_ROI=[];

no_points_pre=-no_points_post_shift;

%training_decisions is 1 if S+ and 2 if S-
%epochs has masks for the following epochs
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR
while (at_end==0)
    next_ii_sp=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    next_ii_sm=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    
    if (isempty(next_ii_sp))&(isempty(next_ii_sm))
        at_end=1;
    else
        
        if isempty(next_ii_sm)
            %This is S+
            next_ii=next_ii_sp;
            if (no_points_post_shift+this_ii+next_ii+ii_span<length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                 ii_sp_post=ii_sp_post+1;
                dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
               
                trial_no=trial_no+1;
                this_ii=this_ii+next_ii+ii_span;

            else
                if (no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
                    at_end=1;
                else
                    this_ii=this_ii+next_ii+no_points_post;
                end
            end
        end
        
        if isempty(next_ii_sp)
            %This is S-
            next_ii=next_ii_sm;
            if (no_points_post_shift+this_ii+next_ii+ii_span<length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                ii_sm_post=ii_sm_post+1;
                dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                
                trial_no=trial_no+1;
                this_ii=this_ii+next_ii+no_points_post;
               
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
                if (no_points_post_shift+this_ii+next_ii+ii_span<length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                    ii_sm_post=ii_sm_post+1;
                    dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                   
                    trial_no=trial_no+1;
                    this_ii=this_ii+next_ii+no_points_post;
                   
                else
                    this_ii=this_ii+next_ii+ii_span;
%                 else
%                     if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||(no_points_post_shift+this_ii+next_ii-ii_span>0)&(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
%                         at_end=1;
%                     else
%                         this_ii=this_ii+next_ii+no_points_post;
%                     end
                end
            else
                %This is S+
                next_ii=next_ii_sp;
                if (no_points_post_shift+this_ii+next_ii+ii_span<length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                    ii_sp_post=ii_sp_post+1;
                    dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                   
                    trial_no=trial_no+1;
                    this_ii=this_ii+next_ii+no_points_post;
                    
                else
%                     if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||(no_points_post_shift+this_ii+next_ii-ii_span>0)&(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs))
%                         at_end=1;
%                     else
                        this_ii=this_ii+next_ii+ii_span;
%                     end
                end
            end
        end
        
        
    end
    
end

which_model_for_traces_loo(which_model_for_traces_loo>trial_no)=trial_no;


%Now let's limit the ROIs to those below p_threshold
noROIs_before_trimming=size(traces,1);
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
dFFs_sm_per_trial_per_ROI=dFFs_sm_per_trial_per_ROI(:,p_value_mask,:);
dFFs_sp_per_trial_per_ROI=dFFs_sp_per_trial_per_ROI(:,p_value_mask,:);
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
    title(['dFF timecourses after p value trimming ' num2str(size(measurements_post,2)) ' ROIs'])
end

pffft=1;




handles_out.dt=dt;






% this_cost=[0 ii_cost;ii_cost 0];
% labels=[];
% timepoint_processed=[];
% correct_predict=[];
% correct_predict_shuffled=[];



N=size(dFFs_sm_per_trial_per_ROI,1)+size(dFFs_sp_per_trial_per_ROI,1);
time_span=[dt:dt:dt*size(dFFs_sm_per_trial_per_ROI,3)]+post_shift+dt;
time_span=time_span-(dt_span/2);

dimensionality=zeros(1,length(time_span));
dimensionalitysp=zeros(1,length(time_span));
dimensionalitysm=zeros(1,length(time_span));



for time_point=1:length(time_span)
    
    %Dimensionality for all trials
    %dFF
    measurements=zeros(N,size(dFFs_sm_per_trial_per_ROI,2));
    measurements(1:size(dFFs_sm_per_trial_per_ROI,1),:)=dFFs_sm_per_trial_per_ROI(:,:,time_point);
    measurements(size(dFFs_sm_per_trial_per_ROI,1)+1:end,:)=dFFs_sp_per_trial_per_ROI(:,:,time_point);
    
    
    %Rows: trials, Columns: electrodes
    Signal=measurements;
    dimensionality(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);

    %Dimensionality for sm
    %dFF
    measurements=zeros(size(dFFs_sm_per_trial_per_ROI,1),size(dFFs_sm_per_trial_per_ROI,2));
    measurements(1:size(dFFs_sm_per_trial_per_ROI,1),:)=dFFs_sm_per_trial_per_ROI(:,:,time_point);
    
    
    %Rows: trials, Columns: electrodes
    Signal=measurements;
    dimensionalitysm(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);

      %Dimensionality for sp
    %dFF
    measurements=zeros(size(dFFs_sp_per_trial_per_ROI,1),size(dFFs_sp_per_trial_per_ROI,2));
    measurements(1:size(dFFs_sp_per_trial_per_ROI,1),:)=dFFs_sp_per_trial_per_ROI(:,:,time_point);
    
    %Rows: trials, Columns: electrodes
    Signal=measurements;
    dimensionalitysp(time_point) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);

end



%Plot the dimensionality
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    
    hold on
    
    plot(time_span,dimensionality,'-b')
    
    
    title(['Dimensionality for no ROIs ' num2str(size(dFFs_sm_per_trial_per_ROI,2))])
    xlabel('Time(sec)')
    ylabel('Dimensionality')

      figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    
    hold on
    
    plot(time_span,dimensionalitysp,'-b')
    
    
    title(['Dimensionality S+ for no ROIs ' num2str(size(dFFs_sm_per_trial_per_ROI,2))])
    xlabel('Time(sec)')
    ylabel('Dimensionality')

      figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    
    hold on
    
    plot(time_span,dimensionalitysm,'-b')
    
    
    title(['Dimensionality S- for no ROIs ' num2str(size(dFFs_sm_per_trial_per_ROI,2))])
    xlabel('Time(sec)')
    ylabel('Dimensionality')
end

handles_out.time_span=time_span;
handles_out.dimensionality=dimensionality;
handles_out.dimensionalitysm=dimensionalitysm;
handles_out.dimensionalitysp=dimensionalitysp;
handles_out.noROIs_before_trimming=noROIs_before_trimming;
handles_out.noROIs_after_trimming=size(dFFs_sm_per_trial_per_ROI,2);


fprintf(1, ['Dimensionality with %d ROIs (original no ROIs %d) = %d\n\n'],size(dFFs_sm_per_trial_per_ROI,2),noROIs_before_trimming,mean(dimensionality(ceil(length(time_span)/2):end)));
fprintf(1, ['Dimensionality S+ with %d ROIs (original no ROIs %d) = %d\n\n'],size(dFFs_sm_per_trial_per_ROI,2),noROIs_before_trimming,mean(dimensionalitysp(ceil(length(time_span)/2):end)));
fprintf(1, ['Dimensionality S- with %d ROIs (original no ROIs %d) = %d\n\n'],size(dFFs_sm_per_trial_per_ROI,2),noROIs_before_trimming,mean(dimensionalitysm(ceil(length(time_span)/2):end)));

pffft=1;