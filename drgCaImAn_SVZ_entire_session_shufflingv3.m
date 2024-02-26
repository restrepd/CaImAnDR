function handles_out=drgCaImAn_SVZ_entire_session_shufflingv3(handles_choices)
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
    post_shift=0.5; %Set to 0 if you want to train with odor on points
%     pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
    MLalgo_to_use=[6]; %Vector with the decoding algorithms you want to use
    ii_cost=3;
    p_threshold=1.1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
    dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures
    perm_before=0; %Permute the labels before running the decoder
    
    handles_choices.pre_perFileName=pre_perFileName;
    handles_choices.pre_perPathName=pre_perPathName;
    handles_choices.processing_algorithm=processing_algorithm;
    handles_choices.post_time=post_time;
    handles_choices.k_fold=k_fold;
    handles_choices.post_shift=post_shift;
    handles_choices.MLalgo_to_use=MLalgo_to_use;
%     handles_choices.pre_time=pre_time;
    handles_choices.p_threshold=p_threshold;
    handles_choices.dt_p_threshold=dt_p_threshold;
    handles_choices.show_figures=show_figures;
    handles_choices.ii_cost=ii_cost;
    handles_choices.perm_before=perm_before;
    
else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
    processing_algorithm=handles_choices.processing_algorithm;
    post_time=handles_choices.post_time;
    k_fold=handles_choices.k_fold;
    post_shift=handles_choices.post_shift;
    MLalgo_to_use=handles_choices.MLalgo_to_use;
%     pre_time=handles_choices.pre_time;
    p_threshold=handles_choices.p_threshold;
    dt_p_threshold=handles_choices.dt_p_threshold;
    show_figures=handles_choices.show_figures;
    ii_cost=handles_choices.ii_cost;
    if isfield(handles_choices,'perm_before')
        perm_before=handles_choices.perm_before;
    else
        perm_before=0;
    end
end

tic

warning('off')

%Restart random seeds
rng('shuffle');

convert_z=0; %Convert dFF traces to z
dt_span=16; %Seconds for per trial traces centered on odor on, this used to be 40, too long
moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
no_shuffles=10; %Number of shuffles for per trial shuffling
window_no=2;
time_windows=[-1 0;
    2 4.1];
  
load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
fprintf(1, 'post_time = %d, p_threshold= %d, post_shift= %d, cost %d\n',post_time,p_threshold,post_shift, ii_cost);

if convert_z==1
    for trace_no=1:size(traces,1)
        traces(trace_no,:)=traces(trace_no,:)/std(traces(trace_no,:));
    end
end

classifier_names{1}='Linear';
classifier_names{2}='SVM';
classifier_names{3}='Bayes';
classifier_names{4}='ANN';
classifier_names{5}='Tree';
classifier_names{6}='GLM';


delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

handles_out.handles_choices=handles_choices;
handles_out.pre_perFileName=pre_perFileName;
handles_out.pre_perPathName=pre_perPathName;
handles_out.post_time=post_time;
handles_out.k_fold=k_fold;
handles_out.post_shift=post_shift;
% handles_out.pre_time=pre_time;
handles_not_out.MLalgo_to_use=MLalgo_to_use;

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
                    'Color',[80/255 194/255 255/255],'LineWidth',1)
            else
                plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                    'Color',[238/255 111/255 179/255],'LineWidth',1)
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
dt=median(time(2:end)-time(1:end-1));
ii_p_threshold=ceil(dt_p_threshold/dt);
no_points_post=floor(post_time/dt);
no_points_post_shift=floor(post_shift/dt);
% no_points_pre=floor(pre_time/dt);
measurements_post=[];
% measurements_pre=[];
epochs_sp_post=zeros(1,length(time));
epochs_sp_pre=zeros(1,length(time));
epochs_sm_post=zeros(1,length(time));
epochs_sm_pre=zeros(1,length(time));
which_model_for_traces_loo=no_odor_trials*ones(1,size(traces,2));

%Do both S+ and S-
at_end=0;
this_ii=0;
ii_post=0;
ii=0;
trial_no=0;

dt_post_which_model=floor(20/dt); %Points that model will be used beyond the training period
ii_span=ceil(dt_span/dt);
dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI_post_shifted=[];
dFFs_sm_per_trial_per_ROI_post_shifted=[];
hit_per_trial=[];
cr_per_trial=[];
miss_per_trial=[];
fa_per_trial=[];

[hit_per_trial,cr_per_trial,dFFs_sp_per_trial_per_ROI_post_shifted,...
    dFFs_sm_per_trial_per_ROI_post_shifted,dFF_per_trial_sm,dFF_per_trial_sp,training_decisions_post,...
    which_model_for_traces_loo,decisions_per_trial,...
    ii_pointer_to_td,epochs_sp_post,measurements_post,...
    ii_post,trial_no,epochs_sm_post,miss_per_trial,fa_per_trial] = ...
    drgCaImAn_parse_out_trialsv2(dt, dt_span,epochs,no_points_post_shift,no_points_post,traces,ii_p_threshold,no_odor_trials);
 

%Calculate percent correct
handles_out.percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);

fprintf(1, 'percent correct behavior = %d\n',handles_out.percent_correct);

which_model_for_traces_loo(which_model_for_traces_loo>trial_no)=trial_no;


%Now let's limit the ROIs to those below p_threshold

p_values=ones(1,size(dFFs_sm_per_trial_per_ROI_post_shifted,2));
for iiROI=1:size(dFFs_sm_per_trial_per_ROI_post_shifted,2)
    dFF_sm=zeros(size(dFFs_sm_per_trial_per_ROI_post_shifted,1),size(dFFs_sm_per_trial_per_ROI_post_shifted,3));
    dFF_sm(:,:)=dFFs_sm_per_trial_per_ROI_post_shifted(:,iiROI,:);
    dFF_sp=zeros(size(dFFs_sp_per_trial_per_ROI_post_shifted,1),size(dFFs_sp_per_trial_per_ROI_post_shifted,3));
    dFF_sp(:,:)=dFFs_sp_per_trial_per_ROI_post_shifted(:,iiROI,:);
    
    [h,p_values(iiROI)]=ttest2(mean(dFF_sp,2),mean(dFF_sm,2));
end

p_value_mask=logical(p_values<=p_threshold);

%Trim the number of ROIs in all matrices
noROIs_before_trimming=size(measurements_post,2);
dFF_per_trial_sp=dFF_per_trial_sp(:,p_value_mask,:);
dFF_per_trial_sm=dFF_per_trial_sm(:,p_value_mask,:);
measurements_post=measurements_post(:,p_value_mask);
% measurements_pre=measurements_pre(:,p_value_mask);
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
% 
% %Euclidean distance
% %For each time point
% for ii_t=1:size(dFF_per_trial_sp,3)
%     jj=0;
%     d=[];
%     for ii_sp=1:size(dFF_per_trial_sp,1)
%         for ii_sm=1:size(dFF_per_trial_sm,1)
%             sum_of_sq=0;
%             for iiROI=1:size(dFF_per_trial_sm,2)
%                 sum_of_sq=sum_of_sq+(dFF_per_trial_sp(ii_sp,iiROI,ii_t)-dFF_per_trial_sm(ii_sm,iiROI,ii_t))^2;
%             end
%             jj=jj+1;
%             d(jj)=sqrt(sum_of_sq);
%         end
%     end
%     dist_euclid(ii_t)=mean(d);
% end
% 
% handles_out.dist_euclid=dist_euclid;

% time_span=time_span(1:end-1)+post_shift;
% handles_out.time_span=time_span;
% handles_out.dist_euclid_zero=mean(dist_euclid((time_span>-20)&(time_span<=0)));
% 
% %Kullback-Leibler divergence
% KLdivergence=zeros(1,size(dFF_per_trial_sp,3));
% for ii_t=1:size(dFF_per_trial_sp,3)
%     
%     z_sp_mean=zeros(1,size(dFF_per_trial_sp,2));
%     these_z_sp=zeros(size(dFF_per_trial_sp,1),size(dFF_per_trial_sp,2));
%     these_z_sp(:,:)=dFF_per_trial_sp(:,:,ii_t);
%     z_sp_mean(1,:)=mean(these_z_sp,1);
%     
%     z_sm_mean=zeros(1,size(dFF_per_trial_sm,2));
%     these_z_sm=zeros(size(dFF_per_trial_sm,1),size(dFF_per_trial_sm,2));
%     these_z_sm(:,:)=dFF_per_trial_sm(:,:,ii_t);
%     z_sm_mean(1,:)=mean(these_z_sm,1);
%     
%     %Use pdist to find all distances for odor 1
%     ii_sp=0;
%     ii_sm=0;
%     distances_sp=zeros(1,size(dFF_per_trial_sp,1));
%     distances_sm=zeros(1,size(dFF_per_trial_sm,1));
%     
%     for ii_sp=1:size(dFF_per_trial_sp,1)
%         odor_queary=zeros(1,size(dFF_per_trial_sp,2));
%         odor_queary(1,:)=these_z_sp(ii_sp,:);
%         all_points=[z_sp_mean; z_sm_mean; odor_queary];
%         all_distances=pdist(all_points);
%         d12=all_distances(1);
%         dq1=all_distances(2);
%         dq2=all_distances(3);
%         distances_sp(1,ii_sp)=(d12^2 + dq1^2 -dq2^2)/(2*d12);
%     end
%     
%     for ii_sm=1:size(dFF_per_trial_sm,1)
%         odor_queary=zeros(1,size(dFF_per_trial_sm,2));
%         odor_queary(1,:)=these_z_sm(ii_sm,:);
%         all_points=[z_sp_mean; z_sm_mean; odor_queary];
%         all_distances=pdist(all_points);
%         d12=all_distances(1);
%         dq1=all_distances(2);
%         dq2=all_distances(3);
%         distances_sm(1,ii_sm)=(d12^2 + dq1^2 -dq2^2)/(2*d12);
%     end
%     
%     
%     %KL divergence
%     num_bins=50;
%     max_d=max([max(distances_sp) max(distances_sm)]);
%     min_d=min([min(distances_sp) min(distances_sm)]);
%     X=[min_d:(max_d-min_d)/(num_bins-1):max_d];
%     
%     p1=zeros(1,length(X));
%     p2=zeros(1,length(X));
%     
%     for ii=1:length(distances_sp)
%         [min_d min_ii]=min(abs(distances_sp(ii)-X));
%         p1(min_ii)=p1(min_ii)+1;
%     end
%     p1=p1/sum(p1);
%     p1(p1==0)=eps; 
%     
%     for ii=1:length(distances_sm)
%         [min_d min_ii]=min(abs(distances_sm(ii)-X));
%         p2(min_ii)=p2(min_ii)+1;
%     end
%     p2=p2/sum(p2);
%     p2(p2==0)=eps;
%     
%     KLdivergence(ii_t) = kldiv(X,p1,p2);
% end
% handles_out.KLdivergence=KLdivergence;

meandFF_per_trial_sp=zeros(size(dFF_per_trial_sp,1),size(dFF_per_trial_sp,3));
meandFF_per_trial_sp(:,:)=mean(dFF_per_trial_sp,2);

meandFF_per_trial_sm=zeros(size(dFF_per_trial_sm,1),size(dFF_per_trial_sm,3));
meandFF_per_trial_sm(:,:)=mean(dFF_per_trial_sm,2);

time_span=[0:dt:dt*(size(dFF_per_trial_sp,3)-1)]-dt_span;

%Plot the mean timecourse for S+ and S-
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .6])
    
    subplot(3,1,1)
    hold on
    
    
    %S+ trials
    try
        CIsp = bootci(1000, @mean, meandFF_per_trial_sp);
        meansp=mean(dFF_per_trial_sp,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        [hlsp, hpsp] = boundedline(time_span',mean(meandFF_per_trial_sp,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
    catch
    end
    
    
    
    
    %S-
    try
        CIsp = bootci(1000, @mean, meandFF_per_trial_sm);
        meansp=mean(dFF_per_trial_sm,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        
        [hlsp, hpsp] = boundedline(time_span',mean(meandFF_per_trial_sm,1)', CIsp', 'cmap',[238/255 111/255 179/255]);
    catch
    end
    
    
    
    
    plot(time_span',mean(meandFF_per_trial_sp,1)', 'Color',[80/255 194/255 255/255]);
    plot(time_span',mean(meandFF_per_trial_sm,1)', 'Color',[238/255 111/255 179/255]);
    
    text(30,0.75,'S-','Color',[80/255 194/255 255/255])
    text(30,0.65,'S+','Color',[0 114/255 178/255])
    
    ylim([0 1])
    
    title(['Mean dFF'])
    xlabel('Time(sec)')
    ylabel('dFF')
    
    subplot(3,1,2)
    hold on
%     
%     %Euclidean distance
%     plot(time_span',(dist_euclid-mean(dist_euclid((time_span>-20)&(time_span<=0))))', 'Color',[238/255 111/255 179/255]);
%     
%     
%     
%     title(['Euclidean distance'])
%     xlabel('Time(sec)')
%     ylabel('Euclidean d')
%     
%     subplot(3,1,3)
%     hold on
%     
%     %KL divergence
%     plot(time_span',KLdivergence', 'Color',[238/255 111/255 179/255]);
%     
%     
%    
%     title(['KL divergence'])
%     xlabel('Time(sec)')
%     ylabel('KL divergence')
%     
    
end


handles_out.meandFFsp=mean(meandFF_per_trial_sp,1);
handles_out.meandFFsm=mean(meandFF_per_trial_sm,1);

training_decisions_post_sh=zeros(1,length(training_decisions_post));
training_decisions_post_sh(1,:)=training_decisions_post(randperm(length(training_decisions_post)));

for ii=1:no_shuffles
    these_shuffled_decisions_per_trial=decisions_per_trial(randperm(length(decisions_per_trial)));
    ww=0;
    for jj=1:length(decisions_per_trial)
        training_decisions_post_sh2(ii,ww+1:ww+no_points_post)=these_shuffled_decisions_per_trial(jj)*ones(1,no_points_post);
        ww=ww+no_points_post;
    end
end

% handles_out.Nall=Nall;
handles_out.dt=dt;
% handles_out.no_points_post=no_points_post;
% handles_out.no_points_post_shift=no_points_post_shift;
% handles_out.no_points_pre=no_points_pre;
% handles_out.measurements_post=measurements_post;
% handles_out.measurements_pre=measurements_pre;
% handles_out.training_decisions_post=training_decisions_post;
% handles_out.epochs_sp_post=epochs_sp_post;
% handles_out.epochs_sp_pre=epochs_sp_pre;
% handles_out.epochs_sm_post=epochs_sm_post;
% handles_out.epochs_sm_pre=epochs_sm_pre;


fprintf(1, ['Training post with %d ROIs (original no ROIs %d)...\n'],size(measurements_post,2),noROIs_before_trimming);



for MLalgo=MLalgo_to_use
    
    this_cost=[0 ii_cost;ii_cost 0];
    labels=[];
    timepoint_processed=[];
    correct_predict=[];
    correct_predict_shuffled=[];
    
    
    
    Nall_post=size(measurements_post,1);
    
    
    
    switch processing_algorithm
        case 1
            %Train with all the post data
            %Store the training data in a table.
            tblTrn=[];
            tblTrn = array2table(measurements_post);
            
            %Store the decisions in Y
            Y=training_decisions_post;
            
            switch MLalgo
                case 1
                    Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                case 2
                    Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                case 3
                    Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                case 4
                    Mdl = fitcnet(tblTrn,Y);
                case 5
                    Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
            end
            
            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
            %Note: for some reason tisdid not work for net when I did this:
            % [label_traces,score] = predict(Mdl,traces');
            % [label_post,score] = predict(Mdl,measurements_post);
            %I had to resolrt to the for loop:
            switch MLalgo
                case 3,4
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label_traces(ii),score] = predict(Mdl,this_time_point);
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label_post(ii),score] = predict(Mdl,this_time_point);
                    end
                otherwise
                    [label_traces,score] = predict(Mdl,traces');
                    [label_post,score] = predict(Mdl,measurements_post);
            end
            
        case 2
            %k-fold cross validation evaluating performance with left-out data
            %Store the training data in a table.
             
            handles_not_out.MLalgo(MLalgo).models=[];
            handles_not_out.MLalgo(MLalgo).processed_succesfully=1;
            points_masked=floor(no_points_post/k_fold);
            which_model=ones(1,size(measurements_post,1));
            for kk=1:k_fold
                
                
                training_mask=ones(size(measurements_post,1),1);
                at_end=0;
                ii=0;
                
                while at_end==0
                    if ii+no_points_post<=size(measurements_post,1)
                        %This is what I originally used, these samples are adjacent
                        %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                        %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                        %
                        %I then changed it to this, spaced by delta_ii
                        delta_ii=floor(no_points_post/points_masked);
                        for ii_pm=1:points_masked
                            training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                            which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                        end
                        
                        ii=ii+no_points_post;
                    else
                        at_end=1;
                    end
                end
                
                these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                these_training_decisions=zeros(1,sum(training_mask));
                
                
                jj=0;
                for ii=1:size(measurements_post,1)
                    if training_mask(ii)==1
                        jj=jj+1;
                        these_training_measurements(jj,:)=measurements_post(ii,:);
                        these_training_decisions(jj)=training_decisions_post(ii);
                    end
                end
                
                tblTrn=[];
                tblTrn = array2table(these_training_measurements);
                
                %Store the decisions in Y
                Y=these_training_decisions;
                
                switch MLalgo
                    case 1
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                    case 2
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                    case 3
                        %The try catch was entered here because of this
                        %error
                        % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                        % A normal distribution cannot be fit for the combination of class 1 and predictor
                        % these_training_measurements107. The data has zero variance.
                        try
                            handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                        catch
                            handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                        end
                    case 4
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                    case 5
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                    case 6
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                end
            end
            
            if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                %Note: for some reason tisdid not work for net when I did this:
                % [label_traces,score] = predict(Mdl,traces');
                % [label_post,score] = predict(Mdl,measurements_post);
                %I had to resolrt to the for loop:
                
                which_model_for_traces=randi(k_fold,1,size(traces,2));
                
                for ii=1:length(which_model)
                    which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                end
                
                if MLalgo==6
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces(ii)).Mdl,this_time_point);
                        if label>0.5
                            label_traces(ii)=1;
                        else
                            label_traces(ii)=0;
                        end
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                        if label>0.5
                            label_post(ii)=1;
                        else
                            label_post(ii)=0;
                        end
                    end
                else
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label_traces(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces(ii)).Mdl,this_time_point);
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label_post(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                    end
                end
                
                
                handles_not_out.MLalgo(MLalgo).label_post=label_post;
                handles_not_out.MLalgo(MLalgo).label_traces=label_traces;
                
                %Now do prediction with shuffled training decisions
                %k-fold cross validation evaluating performance with left-out data
                %Store the training data in a table.
                
                handles_not_out.MLalgo(MLalgo).sh_models=[];
                points_masked=floor(no_points_post/k_fold);
                which_model=ones(1,size(measurements_post,1));
                for kk=1:k_fold
                    
                    
                    training_mask=ones(size(measurements_post,1),1);
                    at_end=0;
                    ii=0;
                    
                    while at_end==0
                        if ii+no_points_post<=size(measurements_post,1)
                            %This is what I originally used, these samples are adjacent
                            %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                            %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                            %
                            %I then changed it to this, spaced by delta_ii
                            delta_ii=floor(no_points_post/points_masked);
                            for ii_pm=1:points_masked
                                training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                            end
                            
                            ii=ii+no_points_post;
                        else
                            at_end=1;
                        end
                    end
                    
                    these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                    these_training_decisions=zeros(1,sum(training_mask));
                    
                    
                    jj=0;
                    for ii=1:size(measurements_post,1)
                        if training_mask(ii)==1
                            jj=jj+1;
                            these_training_measurements(jj,:)=measurements_post(ii,:);
                            these_training_decisions(jj)=training_decisions_post_sh(ii);
                        end
                    end
                    
                    tblTrn=[];
                    tblTrn = array2table(these_training_measurements);
                    
                    %Store the decisions in Y
                    Y=these_training_decisions;
                    
                    handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;
                    
                    switch MLalgo
                        case 1
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        case 2
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            %The try catch was entered here because of this
                            %error
                            % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                            % A normal distribution cannot be fit for the combination of class 1 and predictor
                            % these_training_measurements107. The data has zero variance.
                            try
                                handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                            catch
                                handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                            end
                        case 4
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                        case 5
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                        case 6
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                    end
                end
                
                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                %Note: for some reason this did not work for net when I did this:
                % [label_traces,score] = predict(Mdl,traces');
                % [label_post,score] = predict(Mdl,measurements_post);
                %I had to resort to the for loop:
                
                which_model_for_traces=randi(k_fold,1,size(traces,2));
                
                for ii=1:length(which_model)
                    which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                end
                
                if MLalgo==6
                    
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces(ii)).Mdl,this_time_point);
                        if label>0.5
                            label_traces_sh(ii)=1;
                        else
                            label_traces_sh(ii)=0;
                        end
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                        if label>0.5
                            label_post_sh(ii)=1;
                        else
                            label_post_sh(ii)=0;
                        end
                    end
                else
                    label_traces_sh=zeros(1,size(traces,2));
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label_traces_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces(ii)).Mdl,this_time_point);
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label_post_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                    end
                end
                
                handles_not_out.MLalgo(MLalgo).label_post=label_post;
                handles_not_out.MLalgo(MLalgo).label_traces=label_traces;
                
                handles_not_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
                handles_not_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
                
                %Now do shuffling on a per trial basis sh2
                for ii_sh=1:no_shuffles
                    
                    
                    %Now do prediction with shuffled training decisions
                    %k-fold cross validation evaluating performance with left-out data
                    %Store the training data in a table.
                    which_model
                    handles_not_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                    points_masked=floor(no_points_post/k_fold);
                    which_model=ones(1,size(measurements_post,1));
                    for kk=1:k_fold
                        
                        
                        training_mask=ones(size(measurements_post,1),1);
                        at_end=0;
                        ii=0;
                        
                        while at_end==0
                            if ii+no_points_post<=size(measurements_post,1)
                                %This is what I originally used, these samples are adjacent
                                %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                                %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                                %
                                %I then changed it to this, spaced by delta_ii
                                delta_ii=floor(no_points_post/points_masked);
                                for ii_pm=1:points_masked
                                    training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                    which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                                end
                                
                                ii=ii+no_points_post;
                            else
                                at_end=1;
                            end
                        end
                        
                        these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                        these_training_decisions=zeros(1,sum(training_mask));
                        
                        
                        jj=0;
                        for ii=1:size(measurements_post,1)
                            if training_mask(ii)==1
                                jj=jj+1;
                                these_training_measurements(jj,:)=measurements_post(ii,:);
                                these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                            end
                        end
                        
                        tblTrn=[];
                        tblTrn = array2table(these_training_measurements);
                        
                        %Store the decisions in Y
                        Y=these_training_decisions;
                        
                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;
                        
                        switch MLalgo
                            case 1
                                try
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                catch
                                    % Error using ClassificationDiscriminant (line 380)
                                    % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                    % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                end
                            case 2
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                            case 3
                                %The try catch was entered here because of this
                                %error
                                % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                % A normal distribution cannot be fit for the combination of class 1 and predictor
                                % these_training_measurements107. The data has zero variance.
                                try
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                catch
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                                end
                            case 4
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                            case 5
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                            case 6
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                        end
                    end
                    
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason this did not work for net when I did this:
                    % [label_traces,score] = predict(Mdl,traces');
                    % [label_post,score] = predict(Mdl,measurements_post);
                    %I had to resort to the for loop:
                    
                    which_model_for_traces=randi(k_fold,1,size(traces,2));
                    
                    for ii=1:length(which_model)
                        which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                    end
                    
                    if MLalgo==6
                        
                        
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_post_sh2(ii_sh,ii)=1;
                            else
                                label_post_sh2(ii_sh,ii)=0;
                            end
                        end
                    else
                        
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label_post_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                        end
                    end
                    
                    handles_not_out.MLalgo(MLalgo).label_post=label_post_sh2;
                    
                end
                
                
            end
        case 3
            %leave one trial out
            %Store the training data in a table.
            
            %If perm_before==1 permutate the labels
            if perm_before==1
                these_perm=randperm(length(training_decisions_post));
                handles_out.these_perm=these_perm;
                training_decisions_post=training_decisions_post(these_perm);
            end

            handles_not_out.MLalgo(MLalgo).models=[];
            handles_not_out.MLalgo(MLalgo).processed_succesfully=1;
            points_masked=floor(no_points_post/k_fold);
            which_model=ones(1,size(measurements_post,1));
            no_trials=ii_post/no_points_post;
            for kk=1:no_trials
                
                
                training_mask=ones(size(measurements_post,1),1);
                training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;
                
                
                these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                these_training_decisions=zeros(1,sum(training_mask));
                
                jj=0;
                for ii=1:size(measurements_post,1)
                    if training_mask(ii)==1
                        jj=jj+1;
                        these_training_measurements(jj,:)=measurements_post(ii,:);
                        these_training_decisions(jj)=training_decisions_post(ii);
                    end
                end
                
 
                tblTrn=[];
                tblTrn = array2table(these_training_measurements);
                
                %Store the decisions in Y
                Y=these_training_decisions;
                
                switch MLalgo
                    case 1
                        try
                            handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        catch
                            % Error using ClassificationDiscriminant (line 380)
                            % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                            % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                            handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                        end
                    case 2
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                    case 3
                        %The try catch was entered here because of this
                        %error
                        % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                        % A normal distribution cannot be fit for the combination of class 1 and predictor
                        % these_training_measurements107. The data has zero variance.
                        try
                            handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                        catch
                            handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                        end
                    case 4
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                    case 5
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                    case 6
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                end
            end
            
            if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                %Note: for some reason this did not work for net when I did this:
                % [label_traces,score] = predict(Mdl,traces');
                % [label_post,score] = predict(Mdl,measurements_post);
                %I had to resort to the for loop:
                
                
                if MLalgo==6
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        %which_model_for_traces_loo uses a leave one out
                        %strategy
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        scores(ii,:)=score;
                        %Please note this is a winner take all
                        if label>0.5
                            label_traces(ii)=1;
                        else
                            label_traces(ii)=0;
                        end
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                        scores_post(ii,:)=score;
                        if label>0.5
                            label_post(ii)=1;
                        else
                            label_post(ii)=0;
                        end
                    end
                else
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        try
                            [label_traces(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            scores(ii,:)=score;
                        catch
                            scores(ii,:)=rand(1,2);
                            if rand(1)>0.5
                                label_traces(ii)=1;
                            else
                                label_traces(ii)=0;
                            end
                        end
                        
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        try
                            [label_post(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                            scores_post(ii,:)=score;
                        catch
                            scores_post(ii,:)=rand(1,2);
                            if rand(1)>0.5
                                label_post(ii)=1;
                            else
                                label_post(ii)=0;
                            end
                        end
                        
                    end
                end
                
                
                handles_not_out.MLalgo(MLalgo).label_post=label_post;
                handles_out.MLalgo(MLalgo).label_traces=label_traces;
                handles_out.MLalgo(MLalgo).scores=scores;
                handles_not_out.MLalgo(MLalgo).scores_post=scores_post;
                
                %Now do prediction with shuffled training decisions
                %k-fold cross validation evaluating performance with left-out data
                %Store the training data in a table.
                
                handles_not_out.MLalgo(MLalgo).sh_models=[];
                points_masked=floor(no_points_post/k_fold);
                which_model=ones(1,size(measurements_post,1));
                for kk=1:no_trials
                    
                    
                    training_mask=ones(size(measurements_post,1),1);
                    training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                    which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;
                    
                    
                    these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                    these_training_decisions=zeros(1,sum(training_mask));
                    
                    
                    jj=0;
                    for ii=1:size(measurements_post,1)
                        if training_mask(ii)==1
                            jj=jj+1;
                            these_training_measurements(jj,:)=measurements_post(ii,:);
                            these_training_decisions(jj)=training_decisions_post_sh(ii);
                        end
                    end
                    
                    tblTrn=[];
                    tblTrn = array2table(these_training_measurements);
                    
                    %Store the decisions in Y
                    Y=these_training_decisions;
                    
                    handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;
                    
                    switch MLalgo
                        case 1
                            try
                                handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                            catch
                                % Error using ClassificationDiscriminant (line 380)
                                % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                            end
                        case 2
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            
                            %The try catch was entered here because of this
                            %error
                            % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                            % A normal distribution cannot be fit for the combination of class 1 and predictor
                            % these_training_measurements107. The data has zero variance.
                            try
                                handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                            catch
                                handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                            end
                        case 4
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                        case 5
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                        case 6
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                    end
                end
                
                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                %Note: for some reason this did not work for net when I did this:
                % [label_traces,score] = predict(Mdl,traces');
                % [label_post,score] = predict(Mdl,measurements_post);
                %I had to resort to the for loop:
                
                
                
                if MLalgo==6
                    
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        scores_sh(ii,:)=score;
                        if label>0.5
                            label_traces_sh(ii)=1;
                        else
                            label_traces_sh(ii)=0;
                        end
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                        scores_post_sh(ii,:)=score;
                        if label>0.5
                            label_post_sh(ii)=1;
                        else
                            label_post_sh(ii)=0;
                        end
                    end
                else
                    label_traces_sh=zeros(1,size(traces,2));
                    for ii=1:size(traces,2)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=traces(:,ii);
                        [label_traces_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        scores_sh(ii,:)=score;
                    end
                    
                    for ii=1:size(measurements_post,1)
                        this_time_point=zeros(1,size(traces,1));
                        this_time_point(1,:)=measurements_post(ii,:);
                        [label_post_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                        scores_post_sh(ii,:)=score;
                    end
                end
                
                
                
                handles_not_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
                handles_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
                handles_out.MLalgo(MLalgo).scores_sh=scores_sh;
                handles_not_out.MLalgo(MLalgo).scores_post_sh=scores_post_sh;
                
                %Now do shuffling on a per trial basis
                for ii_sh=1:no_shuffles
                    
                    
                    %Now do prediction with shuffled training decisions
                    %k-fold cross validation evaluating performance with left-out data
                    %Store the training data in a table.
                    
                    handles_not_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                    points_masked=floor(no_points_post/k_fold);
                    which_model=ones(1,size(measurements_post,1));
                    for kk=1:no_trials
                        
                        
                        training_mask=ones(size(measurements_post,1),1);
                        training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                        which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;
                        
                        these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                        these_training_decisions=zeros(1,sum(training_mask));
                        
                        
                        jj=0;
                        for ii=1:size(measurements_post,1)
                            if training_mask(ii)==1
                                jj=jj+1;
                                these_training_measurements(jj,:)=measurements_post(ii,:);
                                these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                            end
                        end
                        
                        tblTrn=[];
                        tblTrn = array2table(these_training_measurements);
                        
                        %Store the decisions in Y
                        Y=these_training_decisions;
                        
                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;
                        
                        switch MLalgo
                            case 1
                                try
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                catch
                                    % Error using ClassificationDiscriminant (line 380)
                                    % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                    % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl  = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                end
                            case 2
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                            case 3
                                %The try catch was entered here because of this
                                %error
                                % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                % A normal distribution cannot be fit for the combination of class 1 and predictor
                                % these_training_measurements107. The data has zero variance.
                                try
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                catch
                                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                                end
                            case 4
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                            case 5
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                            case 6
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                        end
                    end
                    
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason this did not work for net when I did this:
                    % [label_traces,score] = predict(Mdl,traces');
                    % [label_post,score] = predict(Mdl,measurements_post);
                    %I had to resort to the for loop:
                    
                    
                    if MLalgo==6
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            scores_sh2(ii_sh,ii,:)=score;
                            if label>0.5
                                label_traces_sh2(ii_sh,ii)=1;
                            else
                                label_traces_sh2(ii_sh,ii)=0;
                            end
                        end
                        
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                            scores_post_sh2(ii_sh,ii,:)=score;
                            if label>0.5
                                label_post_sh2(ii_sh,ii)=1;
                            else
                                label_post_sh2(ii_sh,ii)=0;
                            end
                        end
                    else
                        
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label_traces_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            scores_sh2(ii_sh,ii,:)=score;
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label_post_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                            scores_post_sh2(ii_sh,ii,:)=score;
                        end
                    end
                    
                    handles_out.MLalgo(MLalgo).label_traces_sh2=label_traces_sh2;
                    handles_not_out.MLalgo(MLalgo).label_post_sh2=label_post_sh2;
                    handles_out.MLalgo(MLalgo).scores_sh2=scores_sh2;
                    handles_not_out.MLalgo(MLalgo).scores_post_sh2=scores_post_sh2;
                    
                end
                
                
            end
    end
    
    if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
        %label is the predicted label, and score is the predicted class
        %posterior probability
        for ii=1:length(training_decisions_post)
            if label_post(ii)==training_decisions_post(ii)
                correct_predict_tr(ii)=1;
            else
                correct_predict_tr(ii)=0;
            end
        end
        
        
        %Calculate wta for windows of no_points_post
        correct_predict_tr_wta=zeros(1,length(training_decisions_post));
        for ii=1:length(training_decisions_post)-no_points_post
            this_correct=zeros(1,no_points_post);
            for jj=1:no_points_post
                if label_post(ii+jj-1)==training_decisions_post(ii+jj-1)
                    this_correct(jj)=1;
                end
            end
            if sum(this_correct)>(no_points_post/2)
                correct_predict_tr_wta(ii+floor(no_points_post/2))=1;
            else
                correct_predict_tr_wta(ii+floor(no_points_post/2))=0;
            end
        end
        
        %Now do shuffled
        %label is the predicted label, and score is the predicted class
        %posterior probability
        for ii=1:length(training_decisions_post_sh)
            if label_post_sh(ii)==training_decisions_post_sh(ii)
                correct_predict_tr_sh(ii)=1;
            else
                correct_predict_tr_sh(ii)=0;
            end
        end
        
        %Calculate wta for windows of no_points_post
        correct_predict_tr_wta_sh=zeros(1,length(training_decisions_post_sh));
        for ii=1:length(training_decisions_post_sh)-no_points_post
            this_correct=zeros(1,no_points_post);
            for jj=1:no_points_post
                if label_post_sh(ii+jj-1)==training_decisions_post_sh(ii+jj-1)
                    this_correct(jj)=1;
                end
            end
            if sum(this_correct)>(no_points_post/2)
                correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=1;
            else
                correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=0;
            end
        end
        
        %Now do shuffled per trial
        correct_predict_tr_wta_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
        correct_predict_tr_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
        for ii_sh=1:no_shuffles
            %label is the predicted label, and score is the predicted class
            %posterior probability
            for ii=1:size(training_decisions_post_sh2,2)
                if label_post_sh2(ii_sh,ii)==training_decisions_post_sh2(ii_sh,ii)
                    correct_predict_tr_sh2(ii_sh,ii)=1;
                else
                    correct_predict_tr_sh2(ii_sh,ii)=0;
                end
            end
            
            %Calculate wta for windows of no_points_post
            
            for ii=1:size(training_decisions_post_sh2,2)-no_points_post
                this_correct=zeros(1,no_points_post);
                for jj=1:no_points_post
                    if label_post_sh2(ii_sh,ii+jj-1)==training_decisions_post_sh2(ii_sh,ii+jj-1)
                        this_correct(jj)=1;
                    end
                end
                if sum(this_correct)>(no_points_post/2)
                    correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=1;
                else
                    correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=0;
                end
            end
            
        end
        
        
        fprintf(1, ['Training accuracy for ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
            ,ii_cost, sum(correct_predict_tr)/length(correct_predict_tr),sum(correct_predict_tr_wta)/length(correct_predict_tr_wta));
        fprintf(1, ['Shuffled training accuracy for ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
            ,ii_cost, sum(correct_predict_tr_sh)/length(correct_predict_tr_sh),sum(correct_predict_tr_wta_sh)/length(correct_predict_tr_wta_sh));
        fprintf(1, ['Training accuracy for shuffled per trial ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
            ,ii_cost, sum(correct_predict_tr_sh2(:))/length(correct_predict_tr_sh2(:)),sum(correct_predict_tr_wta_sh2(:))/length(correct_predict_tr_wta_sh2(:)));
        fprintf(1, ['Mean label trace %d, variance %d\n'],mean(label_traces),var(label_traces));
        
        
        handles_not_out.MLalgo(MLalgo).correct_predict_tr=correct_predict_tr;
        handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta=correct_predict_tr_wta;
        handles_out.MLalgo(MLalgo).accuracy_tr=sum(correct_predict_tr)/length(correct_predict_tr);
        handles_out.MLalgo(MLalgo).accuracy_tr_wta=sum(correct_predict_tr_wta)/length(correct_predict_tr_wta);
        
        handles_not_out.MLalgo(MLalgo).correct_predict_tr_sh=correct_predict_tr_sh;
        handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta_sh=correct_predict_tr_wta_sh;
        handles_out.MLalgo(MLalgo).accuracy_tr_sh=sum(correct_predict_tr_sh)/length(correct_predict_tr_sh);
        handles_out.MLalgo(MLalgo).accuracy_tr_wta_sh=sum(correct_predict_tr_wta_sh)/length(correct_predict_tr_wta_sh);
        
        handles_not_out.MLalgo(MLalgo).correct_predict_tr_sh2=correct_predict_tr_sh2;
        handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta_sh2=correct_predict_tr_wta_sh2;
        handles_out.MLalgo(MLalgo).accuracy_tr_sh2=sum(correct_predict_tr_sh2(:))/length(correct_predict_tr_sh2(:));
        handles_out.MLalgo(MLalgo).accuracy_tr_wta_sh2=sum(correct_predict_tr_wta_sh2(:))/length(correct_predict_tr_wta_sh2(:));
        
        handles_not_out.MLalgo(MLalgo).mean_label_traces=mean(label_traces);
        handles_not_out.MLalgo(MLalgo).var_label_traces=var(label_traces);
        
        moving_mean_label_traces = movmean(label_traces,moving_mean_n);
        handles_not_out.MLalgo(MLalgo).label_traces=label_traces;
        
        %Now let's do the carpentry
        %             moving_mean_label_traces_sh = movmean(mean(label_traces_sh2),moving_mean_n);
        
        moving_mean_label_traces_sh2 = movmean(label_traces_sh2,moving_mean_n);
        handles_not_out.MLalgo(MLalgo).label_traces_sh2=label_traces_sh2;
        moving_mean_label_traces_sh = movmean(label_traces_sh,moving_mean_n);
        handles_not_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
        
        if show_figures==1
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig = figure(figNo);
            
            set(hFig, 'units','normalized','position',[.05 .1 .85 .3])
            
            hold on
            
            
            %                 CIsm = bootci(1000, @mean, moving_mean_label_traces_sh2);
            %                 meansm=mean(moving_mean_label_traces_sh2,1);
            %                 CIsm(1,:)=meansm-CIsm(1,:);
            %                 CIsm(2,:)=CIsm(2,:)-meansm;
            %
            %                 %S- Proficient
            %                 [hlsm, hpsm] = boundedline(time',mean(moving_mean_label_traces_sh2,1)', CIsm', 'cmap',[80/255 194/255 255/255]);
            %
            
            per95=prctile(moving_mean_label_traces_sh2(:),95);
            per5=prctile(moving_mean_label_traces_sh2(:),5);
            CIsh=[mean(moving_mean_label_traces_sh2(:))-per5 per95-mean(moving_mean_label_traces_sh2(:))]';
            [hlCR, hpCR] = boundedline([time(1) time(end)],[mean(moving_mean_label_traces_sh2(:)) mean(moving_mean_label_traces_sh2(:))], CIsh', 'cmap',[80/255 194/255 255/255]);
            
            
            %                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
            plot(time,moving_mean_label_traces,'-k','LineWidth',1)
            plot(time,1.1*((epochs==8)+(epochs==9)),'-b')
            plot(time,1.1*((epochs==6)+(epochs==7)),'-r')
            
            %                 plot(time,moving_mean_label_traces_sh(1,:),'-b')
            
            ylim([-0.2 1.2])
            title(['Label prediction for entire session for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
            
        end
        
        handles_out.time=time;
        
        %Calculate Shanon enthropy
        pone=sum(label_traces==1)/length(label_traces);
        pzero=sum(label_traces==0)/length(label_traces);
        handles_out.MLalgo(MLalgo).shannon_e=-pzero*log2(pzero) - pone*log2(pone);
        
        pone=sum(label_traces_sh==1)/length(label_traces_sh);
        pzero=sum(label_traces_sh==0)/length(label_traces_sh);
        handles_out.MLalgo(MLalgo).shannon_e_sh=-pzero*log2(pzero) - pone*log2(pone);
        
        for ii=1:no_shuffles
            pone=sum(label_traces_sh2(ii,:)==1)/size(label_traces_sh2,2);
            pzero=sum(label_traces_sh2(ii,:)==0)/size(label_traces_sh2,2);
            handles_out.MLalgo(MLalgo).shannon_e_sh2(ii)=-pzero*log2(pzero) - pone*log2(pone);
        end
        fprintf(1, ['Shannon entropy %d, shuffled 1 %d, shuffled 2 %d\n'],handles_out.MLalgo(MLalgo).shannon_e...
            ,handles_out.MLalgo(MLalgo).shannon_e_sh,mean(handles_out.MLalgo(MLalgo).shannon_e_sh2));
        
        %Now let's do accounting and show it in a bar graph
        
        %post Splus
        post_label_sp=label_traces(logical(epochs_sp_post));
        points_per_cut=no_points_post;
        no_cuts=floor(length(post_label_sp)/points_per_cut);
        mean_post_label_sp=zeros(1,no_cuts);
        for ii=1:no_cuts
            mean_post_label_sp(ii)=mean(post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        end
        %     sh_post_label_sp=sum(sh_label_traces,1)/size(sh_label_traces,1);
        %     sh_post_label_sp=sh_post_label_sp(logical(epochs_sp_post));
        %     mean_sh_post_label_sp=zeros(1,no_cuts);
        %     for ii=1:no_cuts
        %         mean_sh_post_label_sp(ii)=mean(sh_post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        %     end
        
        %post Sminus
        post_label_sm=label_traces(logical(epochs_sm_post));
        points_per_cut=no_points_post;
        no_cuts=floor(length(post_label_sm)/points_per_cut);
        mean_post_label_sm=zeros(1,no_cuts);
        for ii=1:no_cuts
            mean_post_label_sm(ii)=mean(post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        end
        %     sh_post_label_sm=sum(sh_label_traces,1)/size(sh_label_traces,1);
        %     sh_post_label_sm=sh_post_label_sm(logical(epochs_sm_post));
        %     mean_sh_post_label_sm=zeros(1,no_cuts);
        %     for ii=1:no_cuts
        %         mean_sh_post_label_sm(ii)=mean(sh_post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        %     end
        
        
%         %pre Sminus
%         pre_label_sm=label_traces(logical(epochs_sm_pre));
%         points_per_cut=no_points_pre;
%         no_cuts=floor(length(pre_label_sm)/points_per_cut);
%         mean_pre_label_sm=zeros(1,no_cuts);
%         for ii=1:no_cuts
%             mean_pre_label_sm(ii)=mean(pre_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
%         end
%         
        
%         %pre Splus
%         pre_label_sp=label_traces(logical(epochs_sp_pre));
%         points_per_cut=no_points_pre;
%         no_cuts=floor(length(pre_label_sp)/points_per_cut);
%         mean_pre_label_sp=zeros(1,no_cuts);
%         for ii=1:no_cuts
%             mean_pre_label_sp(ii)=mean(pre_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
%         end
        
        %     sh_pre_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
        %     sh_pre_label=sh_pre_label(logical(epochs_sm_pre+epochs_sp_pre));
        %     mean_sh_pre_label=zeros(1,no_cuts);
        %     for ii=1:no_cuts
        %         mean_sh_pre_label(ii)=mean(sh_pre_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        %     end
        
        %all
        all_label=label_traces;
        %points_per_cut=no_points_pre;
        points_per_cut=no_points_post;
        no_cuts=floor(length(all_label)/points_per_cut);
        mean_all_label=zeros(1,no_cuts);
        for ii=1:no_cuts
            mean_all_label(ii)=mean(all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        end
        %     sh_all_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
        %     mean_sh_all_label=zeros(1,no_cuts);
        %     for ii=1:no_cuts
        %         mean_sh_all_label(ii)=mean(sh_all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
        %     end
        
        if show_figures==1
            %Note that pre and post are refrenced to the start of the training period
            edges=[0:0.033:1.2];
            rand_offset=0.8;
            
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig = figure(figNo);
            
            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
            
            hold on
            
            bar_offset=1;
            
            %                 %S- pre
            %                 bar_offset=1;
            %
            %                 bar(bar_offset,mean(mean_pre_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            %
            %                 %Violin plot
            %                 [mean_out, CIout]=drgViolinPoint(mean_pre_label_sm...
            %                     ,edges,bar_offset,rand_offset,'k','k',3);
            %
            %                 bar_offset=bar_offset+1;
            %
            %                 %S+ pre
            %                 bar(bar_offset,mean(mean_pre_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            %
            %                 %Violin plot
            %                 [mean_out, CIout]=drgViolinPoint(mean_pre_label_sp...
            %                     ,edges,bar_offset,rand_offset,'k','k',3);
            %
            %                 bar_offset=bar_offset+2;
            
            %S- post
            bar(bar_offset,mean(mean_post_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(mean_post_label_sm...
                ,edges,bar_offset,rand_offset,'k','k',3);
            bar_offset=bar_offset+1;
            
            %S+ post
            bar(bar_offset,mean(mean_post_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(mean_post_label_sp...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            bar_offset=bar_offset+2;
            bar(bar_offset,mean(mean_all_label),'LineWidth', 3,'EdgeColor','none','FaceColor','m')
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(mean_all_label...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            xticks([1 2 4 5 7])
            xticklabels({'S- post', 'S+ post','Entire session'})
            title(['Label prediction for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
            
        end
        
        %Now show average labels for S+, S-, etc for 30 sec
        handles_not_out.MLalgo(MLalgo).mean_all_label=mean_all_label;
%         handles_not_out.MLalgo(MLalgo).mean_pre_label_sm=mean_pre_label_sm;
%         handles_not_out.MLalgo(MLalgo).mean_pre_label_sp=mean_pre_label_sp;
        handles_not_out.MLalgo(MLalgo).mean_post_label_sm=mean_post_label_sm;
        handles_not_out.MLalgo(MLalgo).mean_post_label_sp=mean_post_label_sp;
        handles_not_out.MLalgo(MLalgo).mean_all_label=mean_all_label;

        %Now find all the trials
        %epochs is a vector of the length of time that gives information on
        %behavior
        % 1=Final Valve
        % 6=Hit (on for the duration of odor on)
        % 7=Miss
        % 8=FA
        % 9=CR
        
        at_end=0;
        ii=1;
        tr_ii=0;
        sp_ii=0;
        hit_ii=0;
        miss_ii=0;
        per_trial_sp_timecourse=[];
        per_trial_hit_timecourse=[];
        per_trial_miss_timecourse=[];
        epoch_before_sp=[];
        sm_ii=0;
        cr_ii=0;
        fa_ii=0;
        per_trial_sm_timecourse=[];
        per_trial_cr_timecourse=[];
        per_trial_fa_timecourse=[];
        epoch_before_sm=[];
        per_trial_scores_sp=[];
        per_trial_scores_sm=[];
        these_hits=[];
        these_miss=[];
        these_sp_hits=[];
        these_sp_miss=[];
        these_crs=[];
        these_fas=[];
        these_sm_crs=[];
        these_sm_fas=[];
        these_sps=[];
        theese_sms=[];

        last_sp_sm=-1;


        while at_end==0
            next_ii=[];
            next_ii_sp=find(epochs_sp_post(ii:end)==1,1,'first');
            next_ii_sm=find(epochs_sm_post(ii:end)==1,1,'first');
            if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
                if next_ii_sp<next_ii_sm
                    next_ii=next_ii_sp;
                    this_sp_sm=1;
                else
                    next_ii=next_ii_sm;
                    this_sp_sm=0;
                end
            else
                if ~isempty(next_ii_sp)
                    next_ii=next_ii_sp;
                    this_sp_sm=1;
                end
                if ~isempty(next_ii_sm)
                    next_ii=next_ii_sm;
                    this_sp_sm=0;
                end
            end

            if ~isempty(next_ii)
                if ((ii+next_ii-ii_span)>0)&((ii+next_ii+ii_span<length(label_traces)))
                    if this_sp_sm==1
                        tr_ii=tr_ii+1;
                        sp_ii=sp_ii+1;

                        these_crs(tr_ii)=0;
                        these_fas(tr_ii)=0;
                        these_sps(tr_ii)=1;
                        these_sms(tr_ii)=0;
                        per_trial_sp_timecourse(sp_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                        per_trial_scores_sp(sp_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
                        epoch_before_sp(sp_ii)=last_sp_sm;
                        last_sp_sm=1;

                        if epochs(ii+next_ii)==6
                            hit_ii=hit_ii+1;
                            per_trial_hit_timecourse(hit_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            these_hits(tr_ii)=1;
                            these_miss(tr_ii)=0;
                            these_sp_hits(sp_ii)=1;
                            these_sp_miss(sp_ii)=0;
                        end
                        if epochs(ii+next_ii)==7
                            miss_ii=miss_ii+1;
                            per_trial_miss_timecourse(miss_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            these_hits(tr_ii)=0;
                            these_miss(tr_ii)=1;
                            these_sp_hits(sp_ii)=0;
                            these_sp_miss(sp_ii)=1;
                        end

                        ii_next_post=find(epochs_sp_post(ii+next_ii:end)==0,1,'first');
                        ii=ii+next_ii+ii_next_post;
                    else
                        sm_ii=sm_ii+1;
                        tr_ii=tr_ii+1;
                        per_trial_sm_timecourse(sm_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                        per_trial_scores_sm(sm_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
                        epoch_before_sm(sm_ii)=last_sp_sm;
                        last_sp_sm=0;

                        these_hits(tr_ii)=0;
                        these_miss(tr_ii)=0;
                        these_sps(tr_ii)=0;
                        theese_sms(tr_ii)=1;
                        if epochs(ii+next_ii)==9
                            cr_ii=cr_ii+1;
                            per_trial_cr_timecourse(cr_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            these_crs(tr_ii)=1;
                            these_fas(tr_ii)=0;
                            these_sm_crs(sm_ii)=1;
                            these_sm_fas(sm_ii)=0;
                        end
                        if epochs(ii+next_ii)==8
                            fa_ii=fa_ii+1;
                            per_trial_fa_timecourse(fa_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            these_crs(tr_ii)=0;
                            these_fas(tr_ii)=1;
                            these_sm_crs(sm_ii)=0;
                            these_sm_fas(sm_ii)=1;
                        end
                        ii_next_post=find(epochs_sm_post(ii+next_ii:end)==0,1,'first');
                        ii=ii+next_ii+ii_next_post;
                    end
                else
                    if  ((ii+next_ii+ii_span>length(label_traces)))
                        at_end=1;
                    else
                        ii=ii+next_ii;
                    end
                end
            else
                at_end=1;
            end

        end

        handles_out.MLalgo(MLalgo).these_sm_crs=these_sm_crs;
        handles_out.MLalgo(MLalgo).these_sm_fas=these_sm_fas;
        handles_out.MLalgo(MLalgo).these_sp_hits=these_sp_hits;
        handles_out.MLalgo(MLalgo).these_sp_miss=these_sp_miss;
        handles_out.MLalgo(MLalgo).these_hits=these_hits;
        handles_out.MLalgo(MLalgo).these_miss=these_miss;
        handles_out.MLalgo(MLalgo).these_fas=these_fas;
        handles_out.MLalgo(MLalgo).these_crs=these_crs;

        pfft=1;
%         
%         at_end=0;
%         ii=1;
%         sp_ii=0;
%         per_trial_sp_timecourse=[];
%         epoch_before_sp=[];
%         sm_ii=0;
%         per_trial_sm_timecourse=[];
%         epoch_before_sm=[];
%         per_trial_scores_sp=[];
%         per_trial_scores_sm=[];
%         
%         last_sp_sm=-1;
%         
%         
%         while at_end==0
%             next_ii=[];
%             next_ii_sp=find(epochs_sp_post(ii:end)==1,1,'first');
%             next_ii_sm=find(epochs_sm_post(ii:end)==1,1,'first');
%             if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
%                 if next_ii_sp<next_ii_sm
%                     next_ii=next_ii_sp;
%                     this_sp_sm=1;
%                 else
%                     next_ii=next_ii_sm;
%                     this_sp_sm=0;
%                 end
%             else
%                 if ~isempty(next_ii_sp)
%                     next_ii=next_ii_sp;
%                     this_sp_sm=1;
%                 end
%                 if ~isempty(next_ii_sm)
%                     next_ii=next_ii_sm;
%                     this_sp_sm=0;
%                 end
%             end
%             
%             if ~isempty(next_ii)
%                 if ((ii+next_ii-ii_span)>0)&((ii+next_ii+ii_span<length(label_traces)))
%                     if this_sp_sm==1
%                         sp_ii=sp_ii+1;
%                         per_trial_sp_timecourse(sp_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
%                         per_trial_scores_sp(sp_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
%                         epoch_before_sp(sp_ii)=last_sp_sm;
%                         last_sp_sm=1;
%                         ii_next_post=find(epochs_sp_post(ii+next_ii:end)==0,1,'first');
%                         ii=ii+next_ii+ii_next_post;
%                     else
%                         sm_ii=sm_ii+1;
%                         per_trial_sm_timecourse(sm_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
%                         per_trial_scores_sm(sm_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
%                         epoch_before_sm(sm_ii)=last_sp_sm;
%                         last_sp_sm=0;
%                         ii_next_post=find(epochs_sm_post(ii+next_ii:end)==0,1,'first');
%                         ii=ii+next_ii+ii_next_post;
%                     end
%                 else
%                     if  ((ii+next_ii+ii_span>length(label_traces)))
%                         at_end=1;
%                     else
%                         ii=ii+next_ii;
%                     end
%                 end
%             else
%                 at_end=1;
%             end
%             
%         end
        
        handles_out.MLalgo(MLalgo).per_trial_sp_timecourse=per_trial_sp_timecourse;
        handles_not_out.MLalgo(MLalgo).epoch_before_sp=epoch_before_sp;
        handles_out.MLalgo(MLalgo).per_trial_sm_timecourse=per_trial_sm_timecourse;
        handles_not_out.MLalgo(MLalgo).epoch_before_sm=epoch_before_sm;
        
        handles_out.MLalgo(MLalgo).per_trial_scores_sp=per_trial_scores_sp;
        handles_out.MLalgo(MLalgo).per_trial_scores_sm=per_trial_scores_sm;
        
        this_moving_mean_n=10;
        moving_mean_per_trial_sp_timecourse = movmean(per_trial_sp_timecourse',this_moving_mean_n)';
        moving_mean_per_trial_sm_timecourse = movmean(per_trial_sm_timecourse',this_moving_mean_n)';
        
%         time_span=[0:dt:dt*size(per_trial_sp_timecourse,2)]-dt_span+dt;
%         time_span=time_span(1:end-1)+post_shift;
        
        %Calculate correct predict
        for ii_tr=1:size(per_trial_sp_timecourse,1)
            for ii_time=1:size(per_trial_sp_timecourse,2)
                if per_trial_sp_timecourse(ii_tr,ii_time)==1
                    this_correct_predict(ii_tr,ii_time)=1;
                else
                    this_correct_predict(ii_tr,ii_time)=0;
                end
            end
        end
        
        for ii_tr=1:size(per_trial_sm_timecourse,1)
            for ii_time=1:size(per_trial_sm_timecourse,2)
                if per_trial_sm_timecourse(ii_tr,ii_time)==0
                    this_correct_predict(ii_tr+sp_ii,ii_time)=1;
                else
                    this_correct_predict(ii_tr+sp_ii,ii_time)=0;
                end
            end
        end
        
        
        %Calculate correct predict shuffled
        ii_plus=0;
        for ww=1:10
            for ii_tr=1:size(per_trial_sp_timecourse,1)
                rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                for ii_time=1:size(per_trial_sp_timecourse,2)
                    if per_trial_sp_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                    else
                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                    end
                end
            end
            
            ii_plus=ii_plus+sp_ii;
            
            for ii_tr=1:size(per_trial_sm_timecourse,1)
                rand_stim=randi([0,1],1,size(per_trial_sm_timecourse,2));
                for ii_time=1:size(per_trial_sm_timecourse,2)
                    if per_trial_sm_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                    else
                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                    end
                end
            end
            ii_plus=ii_plus+sm_ii;
        end
        
        handles_out.MLalgo(MLalgo).this_correct_predict=this_correct_predict;
        handles_out.MLalgo(MLalgo).this_correct_predict_sh=this_correct_predict_sh;
        
        if show_figures==1
            %Plot the prediction for S+ and S- and the correct predictions
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig = figure(figNo);
            
            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
            
%             subplot(2,1,1)
            hold on
            
            CIsm = bootci(1000, @mean, moving_mean_per_trial_sm_timecourse);
            meansm=mean(moving_mean_per_trial_sm_timecourse,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;
            
            [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sm_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
            
            CIsp = bootci(1000, @mean, moving_mean_per_trial_sp_timecourse);
            meansp=mean(moving_mean_per_trial_sp_timecourse,1);
            CIsp(1,:)=meansp-CIsp(1,:);
            CIsp(2,:)=CIsp(2,:)-meansp;
            
            
            [hlsp, hpsp] = boundedline(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', CIsp', 'cmap',[0 114/255 178/255]);
            
            plot(time_span',mean(moving_mean_per_trial_sm_timecourse,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
            plot(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
            
            
            xlim([-7 15])
            ylim([0 1])

            this_ylim=ylim;


            %FV
            plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'LineWidth',5, Color=[0.9 0.9 0.9])
            plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

            %Odor on markers
            plot([0 0],this_ylim,'-k')
            odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')


            title(['Label prediction per trial for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
            xlabel('Time(sec)')
            ylabel('Label prediction, S+=1, S-=0')

            %plot accuracy
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);

            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
            hold on

            CIsm = bootci(1000, @mean, this_correct_predict_sh);
            meansm=mean(this_correct_predict_sh,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');


            CIsm = bootci(1000, @mean, this_correct_predict);
            meansm=mean(this_correct_predict,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);

            plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
            plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);

            text(30,0.75,'Shuffled','Color','k')
            text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])

            xlim([-7 15])
            ylim([0.2 1])

            this_ylim=ylim;


            %FV
            plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'LineWidth',5, Color=[0.9 0.9 0.9])
            plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

            %Odor on markers
            plot([0 0],this_ylim,'-k')
            odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')


            title(['Accuracy v3 for ' classifier_names{MLalgo} ' trained from ' num2str(post_shift) ' to ' num2str(post_shift+post_time)])
            xlabel('Time(sec)')
            ylabel('Accuracy')

            %Plot the posterior probabilities for Sp (scores,:,2) and Sm (scores(:,1))
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig = figure(figNo);
            
            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
            
            subplot(2,1,1)
            hold on
            try
                these_scores_sm=zeros(size(per_trial_scores_sp,1),size(per_trial_scores_sp,3));
                these_scores_sm(:,:)=per_trial_scores_sp(:,1,:);
                CIsm = bootci(1000, @mean, these_scores_sm);
                meansm=mean(these_scores_sm,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);

                these_scores_sp=zeros(size(per_trial_scores_sp,1),size(per_trial_scores_sp,3));
                these_scores_sp(:,:)=per_trial_scores_sp(:,2,:);
                CIsp = bootci(1000, @mean, these_scores_sp);
                meansp=mean(these_scores_sp,1);
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;


                [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);


                plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');

            catch
            end
            
            text(30,0.75,'S-','Color',[158/255 31/255 99/255])
            text(30,0.65,'S+','Color',[0 114/255 178/255])
            
            ylim([0 1])
            xlim([-10 20])
            title(['Posterior probability for S+ or S- prediction for S+ trials ' classifier_names{MLalgo} ' ' num2str(p_threshold)])
            xlabel('Time(sec)')
            ylabel('Posterior probability')
            
            subplot(2,1,2)
            hold on
            
            try
                these_scores_sm=zeros(size(per_trial_scores_sm,1),size(per_trial_scores_sm,3));
                these_scores_sm(:,:)=per_trial_scores_sm(:,1,:);
                CIsm = bootci(1000, @mean, these_scores_sm);
                meansm=mean(these_scores_sm,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);

                these_scores_sp=zeros(size(per_trial_scores_sm,1),size(per_trial_scores_sm,3));
                these_scores_sp(:,:)=per_trial_scores_sm(:,2,:);
                CIsp = bootci(1000, @mean, these_scores_sp);
                meansp=mean(these_scores_sp,1);
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;


                [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);

                plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
            catch
            end
            text(30,0.75,'S-','Color',[158/255 31/255 99/255])
            text(30,0.65,'S+','Color',[0 114/255 178/255])
            
            ylim([0 1])
            xlim([-10 20])
            title(['Posterior probability for S+ or S- prediction for S- trials ' classifier_names{MLalgo} ' ' num2str(p_threshold)])
            xlabel('Time(sec)')
            ylabel('Posterior probability')

            %Plot error trials
            %Plot the prediction for S+ and S- and the correct predictions
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig = figure(figNo);
            
            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
            
%             subplot(2,1,1)
            hold on
            
            %Plot CRs
            if sum(these_sm_crs)>=3
                CIsm = bootci(1000, @mean, moving_mean_per_trial_sm_timecourse(logical(these_sm_crs),:));
                meansm=mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_crs),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_crs),:),1)', CIsm', '-b');
            else
                plot(time_span',mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_crs),:),1)','-b')
            end

            %Plot FAs
            if sum(these_sm_fas)>=3
                CIsm = bootci(1000, @mean, moving_mean_per_trial_sm_timecourse(logical(these_sm_fas),:));
                meansm=mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_fas),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_fas),:),1)', CIsm','-m');
            else
                plot(time_span',mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_fas),:),1)','-m')
            end
            
        %Plot Hits
            if sum(these_sp_hits)>=3
                CIsm = bootci(1000, @mean, moving_mean_per_trial_sp_timecourse(logical(these_sp_hits),:));
                meansm=mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_hits),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_hits),:),1)', CIsm', '-r');
            else
                plot(time_span',mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_hits),:),1)','-r')
            end

            %Plot Miss
            if sum(these_sp_miss)>=3
                CIsm = bootci(1000, @mean, moving_mean_per_trial_sp_timecourse(logical(these_sp_miss),:));
                meansm=mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_miss),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_miss),:),1)', CIsm','-c');
            else
                plot(time_span',mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_miss),:),1)','-c')
            end

            plot(time_span',mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_crs),:),1)','-b')
            plot(time_span',mean(moving_mean_per_trial_sm_timecourse(logical(these_sm_fas),:),1)','-m')
            plot(time_span',mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_hits),:),1)','-r')
            plot(time_span',mean(moving_mean_per_trial_sp_timecourse(logical(these_sp_miss),:),1)','-c')
            
            
            xlim([-7 15])
            ylim([0 1])

            this_ylim=ylim;


            %FV
            plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'LineWidth',5, Color=[0.9 0.9 0.9])
            plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

            %Odor on markers
            plot([0 0],this_ylim,'-k')
            odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')


            title(['Prediction, hit=red, miss=cyan, fa=magenta, cr=blue'])
            xlabel('Time(sec)')
            ylabel('Prediction, S+=1, S-=0')

            %plot accuracy
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);

            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
            hold on

            CIsm = bootci(1000, @mean, this_correct_predict_sh);
            meansm=mean(this_correct_predict_sh,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');

            if sum(these_hits)>=3
                CIsm = bootci(1000, @mean, this_correct_predict(logical(these_hits),:));
                meansm=mean(this_correct_predict(logical(these_hits),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict(logical(these_hits),:),1)', CIsm', '-r');
            else
                plot(time_span',mean(this_correct_predict(logical(these_hits),:),1)', '-r');
            end

            if sum(these_miss)>=3
                CIsm = bootci(1000, @mean, this_correct_predict(logical(these_miss),:));
                meansm=mean(this_correct_predict(logical(these_miss),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict(logical(these_miss),:),1)', CIsm', '-c');
            else
                plot(time_span',mean(this_correct_predict(logical(these_miss),:),1)', '-c');
            end

            if sum(these_crs)>=3
                CIsm = bootci(1000, @mean, this_correct_predict(logical(these_crs),:));
                meansm=mean(this_correct_predict(logical(these_crs),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict(logical(these_crs),:),1)', CIsm', '-b');
            else
                plot(time_span',mean(this_correct_predict(logical(these_crs),:),1)', '-b');
            end

            if sum(these_fas)>=3
                CIsm = bootci(1000, @mean, this_correct_predict(logical(these_fas),:));
                meansm=mean(this_correct_predict(logical(these_fas),:),1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;

                [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict(logical(these_fas),:),1)', CIsm', '-m');
            else
                plot(time_span',mean(this_correct_predict(logical(these_fas),:),1)', '-m');
            end



            plot(time_span',mean(this_correct_predict(logical(these_hits),:),1)', '-r');
            plot(time_span',mean(this_correct_predict(logical(these_miss),:),1)', '-c');
            plot(time_span',mean(this_correct_predict(logical(these_crs),:),1)', '-b');
            plot(time_span',mean(this_correct_predict(logical(these_fas),:),1)', '-m');

            plot(time_span',mean(this_correct_predict_sh,1)', '-k');

            xlim([-7 15])
            ylim([0.2 1])

            this_ylim=ylim;


            %FV
            plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'LineWidth',5, Color=[0.9 0.9 0.9])
            plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

            %Odor on markers
            plot([0 0],this_ylim,'-k')
            odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')


            title(['Accuracy v3 hit=red, miss=cyan, cr=blue, fa=magenta'])
            xlabel('Time(sec)')
            ylabel('Accuracy')
            
        end
    else
        fprintf(1, [classifier_names{MLalgo} ' was not processed succesfully\n']);
    end
    
    fprintf(1,['Accuracy for ' classifier_names{MLalgo} ' = %d\n\n'],mean(mean(handles_out.MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),1)))

end

handles_out.time_span=time_span;

% this_correct_predict=handles_out.MLalgo(MLalgo).this_correct_predict;
%  meansm=mean(this_correct_predict,1);
%  odor_window=[2 4.1]; %Note: This is the time window I use in drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse_v2
%     odor_acc=mean(meansm((time_span>=odor_window(1))&(time_span<=odor_window(2))));
%     fprintf(1, ['Accuracy for odor window  is %d\n'],...
%     odor_acc);

fprintf(1,'Elapsed time (hr) %d\n\n',toc/(60*60))

pffft=1;