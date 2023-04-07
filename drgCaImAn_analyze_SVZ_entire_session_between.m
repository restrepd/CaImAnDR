function handles_out3=drgCaImAn_analyze_SVZ_entire_session_between(handles_choices)
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

    %Load pre_per_dec2 file
    [pre_per_dec2_FileName,pre_per_dec2_PathName] = uigetfile({'*pre_per_dec2.mat'},'Select the .m file with all the choices for analysis');

    %Load pre_per file
    [pre_per_FileName,pre_per_PathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');

    
    iiMLalgo=6; %Vector with the decoding algorithms you want to use
  
    p_threshold=1.1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
  
    show_figures=1; %Show the figures
else
    close all
    pre_per_PathName=handles_choices.pre_per_PathName;
    pre_per_FileName=handles_choices.pre_per_FileName;
    pre_per_dec2_PathName=handles_choices.PathName_Out;
    pre_per_dec2_FileName=[pre_per_FileName(1:end-4) '_dec2.mat'];
%     processing_algorithm=handles_choices.processing_algorithm;
%     post_time=handles_choices.post_time;
%     k_fold=handles_choices.k_fold;
%     post_shift=handles_choices.post_shift;
    iiMLalgo=handles_choices.MLalgo_to_use;
%     pre_time=handles_choices.pre_time;
    p_threshold=handles_choices.p_threshold;
%     dt_p_threshold=handles_choices.dt_p_threshold;
    show_figures=handles_choices.show_figures;
%     show_figures=1; %I am hard coding this to troubleshoot
%     ii_cost=handles_choices.ii_cost;
end

warning('off')

figNo=0;


%Load pre_per and show traces

load([pre_per_PathName pre_per_FileName])

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
 
load([pre_per_dec2_PathName pre_per_dec2_FileName])
fprintf(1, ['\ndrgCaImAn_analyze_SVZ_entire_session run for ' pre_per_dec2_FileName '\n\n']);


classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';
classifier_names{6}='Binomial glm';


delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

ii_thr=find(handles_out.handles.p_threshold==p_threshold);
time_span=handles_out.ii_out(ii_thr).handles_out.time_span;
handles_out3.time_span=time_span;


accuracy_tr=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).accuracy_tr;
accuracy_tr_wta=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).accuracy_tr_wta;
accuracy_tr_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh;
accuracy_tr_wta_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).accuracy_tr_wta_sh;
accuracy_tr_sh2=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh2;
accuracy_tr_wta_sh2=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).accuracy_tr_wta_sh2;

fprintf(1, ['Training accuracy for ' classifier_names{iiMLalgo} ' is %d, wta accuracy is %d\n']...
    ,accuracy_tr,accuracy_tr_wta);
fprintf(1, ['Shuffled training accuracy for ' classifier_names{iiMLalgo} ' is %d, wta accuracy is %d\n']...
    ,accuracy_tr_sh,accuracy_tr_wta_sh);
fprintf(1, ['Training accuracy for shuffled per trial for ' classifier_names{iiMLalgo} ' is %d, wta accuracy is %d\n']...
    ,accuracy_tr_sh2,accuracy_tr_wta_sh2);
fprintf(1, ['Percent correct for behavior  is %d\n'],...
    handles_out.ii_out(1).handles_out.percent_correct);

handles_out3.percent_correct=handles_out.ii_out(1).handles_out.percent_correct;

%Plot the timecourses for accuracy and label

per_trial_sp_timecourse=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).per_trial_sp_timecourse;
per_trial_sm_timecourse=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).per_trial_sm_timecourse;
this_correct_predict=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).this_correct_predict;
this_correct_predict_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh;
per_trial_scores_sp=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).per_trial_scores_sp;
per_trial_scores_sm=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).per_trial_scores_sm;

this_moving_mean_n=10;
moving_mean_n=this_moving_mean_n;
moving_mean_per_trial_sp_timecourse = movmean(per_trial_sp_timecourse',this_moving_mean_n)';
moving_mean_per_trial_sm_timecourse = movmean(per_trial_sm_timecourse',this_moving_mean_n)';


if show_figures==1
    %Plot the prediction for S+ and S- and the accuracy
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .4 .4])

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

    text(30,0.75,'S-','Color',[158/255 31/255 99/255])
    text(30,0.85,'S+','Color',[0 114/255 178/255])

    ylim([0 1.1])
    xlim([-7 15])
    title(['Label prediction per trial for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Label prediction, S+=1, S-=0')
end

handles_out3.moving_mean_per_trial_sm_timecourse=moving_mean_per_trial_sm_timecourse;
handles_out3.moving_mean_per_trial_sp_timecourse=moving_mean_per_trial_sp_timecourse;
handles_out3.moving_time_span=time_span;

if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .4 .7])
    hold on

    y_shift=0;
    for ii=1:size(moving_mean_per_trial_sm_timecourse,1)
        plot(time_span',moving_mean_per_trial_sm_timecourse(ii,:)+y_shift,'Color',[158/255 31/255 99/255]);
        y_shift=y_shift+1.5;
    end

    for ii=1:size(moving_mean_per_trial_sp_timecourse,1)
        plot(time_span',moving_mean_per_trial_sp_timecourse(ii,:)+y_shift,'Color',[0 114/255 178/255]);
        y_shift=y_shift+1.5;
    end

    plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
    plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);

    text(30,0.75,'Shuffled','Color','k')
    text(30,0.85,'S+ vs S-','Color',[0 114/255 178/255])


    xlim([-7 15])
    title(['Prediction for each trial for ' classifier_names{iiMLalgo}])
    xlabel('Time(sec)')
    ylabel('Prediction')

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .4 .4])
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
    text(30,0.85,'S+ vs S-','Color',[0 114/255 178/255])

    ylim([0 1.1])
    this_ylim=ylim;
    %Place epoch markers

    %FV
    rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor
    rectangle(Position=[0,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(delta_odor),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0 0 0], EdgeColor=[0 0 0])
    plot([0 0],[this_ylim],'-k')
    plot([mean(delta_odor) mean(delta_odor)],[this_ylim],'-k')

    %Reinforcement
    rectangle(Position=[mean(delta_odor_on_reinf_on),this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(delta_reinf),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[1 0 0], EdgeColor=[1 0 0])
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[this_ylim],'-r')
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[this_ylim],'-r')



    xlim([-7 15])
    title(['Accuracy for ' classifier_names{iiMLalgo}])
    xlabel('Time(sec)')
    ylabel('Accuracy')
end


odor_window=[2 4.1]; %Note: This is the time window I use in drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse_v2
meansm=mean(this_correct_predict,1);
odor_acc=mean(meansm((time_span>=odor_window(1))&(time_span<=odor_window(2))));
fprintf(1, ['Accuracy for odor window  is %d\n'],...
    odor_acc);

handles_out3.odor_acc=odor_acc;



%Now let's take a look at the entire session, is the animal musing about odorants?!
label_traces=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).label_traces;
label_traces_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).label_traces_sh;
label_traces_sh2=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).label_traces_sh2;


% label_traces_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).label_traces_sh;

moving_mean_label_traces = movmean(label_traces,moving_mean_n);
moving_mean_label_traces_sh2 = movmean(label_traces_sh2,moving_mean_n);
moving_mean_label_traces_sh = movmean(label_traces_sh,moving_mean_n);

time=handles_out.ii_out(ii_thr).handles_out.time;
sp_times=handles_out.ii_out(ii_thr).handles_out.sp_times;
sm_times=handles_out.ii_out(ii_thr).handles_out.sm_times;

if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .85 .6])

    subplot(3,1,1)
    hold on




    %                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
    plot(time,moving_mean_label_traces,'-k','LineWidth',1)



    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end

    for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
    end



    ylim([-0.2 1.2])
    title('original data')

    subplot(3,1,2)
    hold on




    plot(time,moving_mean_label_traces_sh,'-m','LineWidth',1)

    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end

    for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
    end



    ylim([-0.2 1.2])
    title('shuffled')

    subplot(3,1,3)
    hold on


    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end

    for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
    end

    ylim([-0.2 1.2])
    title('trials')

    sgtitle(['Label prediction for entire session for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
end

%Now find the changes in prediction that take place between trials
%Note that "between trial" time period is defined as the time period
%after odor off when the prediction goes back down to baseline
%(<percentile 5) up to next odor on -pre_time
per95=prctile(moving_mean_label_traces_sh2(:),95);
per5=prctile(moving_mean_label_traces_sh2(:),5);
all_trial_times=[sp_times sm_times];
all_trial_times=sortrows(all_trial_times');

window_dt=[-7 15];
spontaneous_trial_times=[];
spontaneous_delta_prediction=[];
before_spontaneous_spm=[];
after_spontaneous_spm=[];
time_before=[];
dt_from_before=[];
ii_spontaneous=0;
pre_time=2;
between_start=[];
between_end=[];
ii_between_periods=0;
no_increase_trial_times=[];
no_increase_delta_prediction=[];
ii_no_inc_included=0;
for ii_trial=1:2:length(all_trial_times)-2
    %find the next trial
    this_ii=find(time>all_trial_times(ii_trial)+delta_odor,1,'first');
    if ~isempty(this_ii)
        %Now find the first baseline period
        delta_ii=find(moving_mean_label_traces(this_ii:end) < per5,1,'first');
        if ~isempty(delta_ii)
            if time(this_ii+delta_ii-1)<all_trial_times(ii_trial+2)
                %This is a between time period
                ii_between_periods=ii_between_periods+1;
                between_start(ii_between_periods)=time(this_ii+delta_ii-1);
                between_end(ii_between_periods)=all_trial_times(ii_trial+2)-pre_time;

                %Find the increase
                delta_ii2=find(moving_mean_label_traces(this_ii+delta_ii-1:end) > per95,1,'first');

                try_no_increase=0;
                if ~isempty(delta_ii2)
                    if time(this_ii+delta_ii-1+delta_ii2-2)<all_trial_times(ii_trial+2)-pre_time
                        ii_spontaneous=ii_spontaneous+1;
                        spontaneous_trial_times(ii_spontaneous)=time(this_ii+delta_ii-1+delta_ii2-2);
                        time_before(ii_spontaneous)=all_trial_times(ii_trial);
                        dt_from_before(ii_spontaneous)=spontaneous_trial_times(ii_spontaneous)-time_before(ii_spontaneous);

                        if sum(sp_times==time_before(ii_spontaneous))>0
                            before_spontaneous_spm(ii_spontaneous)=1;
                        else
                            before_spontaneous_spm(ii_spontaneous)=0;
                        end

                        if ii_trial+2<=length(all_trial_times)
                            if sum(sp_times==all_trial_times(ii_trial+2))>0
                                after_spontaneous_spm(ii_spontaneous)=1;
                            else
                                after_spontaneous_spm(ii_spontaneous)=0;
                            end
                        else
                            after_spontaneous_spm(ii_spontaneous)=NaN;
                        end

                        spontaneous_delta_prediction(ii_spontaneous,1:sum((time>=-7+spontaneous_trial_times(ii_spontaneous))&(time<=15+spontaneous_trial_times(ii_spontaneous))))...
                            =moving_mean_label_traces((time>=-7+spontaneous_trial_times(ii_spontaneous))&(time<=15+spontaneous_trial_times(ii_spontaneous)));

                        pffft=1;
                    else
                        try_no_increase=1;
                    end
                else
                    try_no_increase=1;
                end

                %If there is no increase between these two trials
                %extract a control between trials
                if try_no_increase==1
                    this_no_increase_time=(time(this_ii+delta_ii-1)+all_trial_times(ii_trial+2))/2;
                    ii_no_inc_included=ii_no_inc_included+1;
                    no_increase_trial_times(ii_no_inc_included)=this_no_increase_time;

                    no_increase_delta_prediction(ii_no_inc_included,1:sum((time>=-7+no_increase_trial_times(ii_no_inc_included))&(time<=15+no_increase_trial_times(ii_no_inc_included))))...
                        =moving_mean_label_traces((time>=-7+no_increase_trial_times(ii_no_inc_included))&(time<=15+no_increase_trial_times(ii_no_inc_included)));
                end
            end
        end
    end
end

if show_figures==1
    %Plot label predictions alone
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .85 .3])

    hold on




    CIsh=[mean(moving_mean_label_traces_sh2(:))-per5 per95-mean(moving_mean_label_traces_sh2(:))]';
    [hlCR, hpCR] = boundedline([time(1) time(end)],[mean(moving_mean_label_traces_sh2(:)) mean(moving_mean_label_traces_sh2(:))], CIsh', 'cmap',[80/255 194/255 255/255]);


    %                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
    plot(time,moving_mean_label_traces,'-k','LineWidth',1)

    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end

    for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
    end
    %                 plot(time,moving_mean_label_traces_sh(1,:),'-b')

    ylim([-0.2 1.2])
    title(['Label prediction for entire session for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
end

dt=time(2)-time(1);
this_time_span=[-7:dt:15];
this_time_span=this_time_span(1:size(spontaneous_delta_prediction,2));

if show_figures==1



    %Plot the spontaneous between trial dFF
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .4 .4])

    hold on




 


    hold on

    if ~isempty(no_increase_delta_prediction)
        if size(no_increase_delta_prediction,1)>2
            CIsm = bootci(1000, @mean, no_increase_delta_prediction(:,1:length(this_time_span)));
            meansm=mean(no_increase_delta_prediction(:,1:length(this_time_span)),1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(this_time_span',mean(no_increase_delta_prediction(:,1:length(this_time_span)),1)', CIsm', 'cmap',[158/255 31/255 99/255]);
        end
    end

    if ~isempty(spontaneous_delta_prediction)
        if size(spontaneous_delta_prediction,1)>2
            CIsm = bootci(1000, @mean, spontaneous_delta_prediction);
            meansm=mean(spontaneous_delta_prediction,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(this_time_span',mean(spontaneous_delta_prediction,1)', CIsm', 'cmap',[0 114/255 178/255]);
        end
    end

    %     for ii=1:size(spontaneous_delta_prediction,1)
    %         plot(this_time_span',spontaneous_delta_prediction(ii,:),  'Color',[150/255 150/255 150/255])
    %     end
    if ~isempty(no_increase_delta_prediction)
        plot(this_time_span',mean(no_increase_delta_prediction(:,1:length(this_time_span)),1)', 'Color',[158/255 31/255 99/255]);
    end

    if ~isempty(spontaneous_delta_prediction)
        plot(this_time_span',mean(spontaneous_delta_prediction,1)', 'Color',[0 114/255 178/255]);
    end

    text(10,0.75,'No change','Color',[158/255 31/255 99/255])
    text(10,0.85,'Spontaneous','Color',[0 114/255 178/255])

    ylim([0 1.1])
    xlim([-7 15])
    title(['Label prediction spontaneous between trials for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Label prediction')
end

handles_out3.this_time_span_spon=this_time_span;
handles_out3.spontaneous_delta_prediction=spontaneous_delta_prediction;
handles_out3.no_increase_delta_prediction=no_increase_delta_prediction;


%Now calculate the histograms for the fraction of S+ predictions
%in between, S+ and S- pre odor and odor, experimental and shuffled
baseline_t=[-2.5 -1.5];
odor_t=[2 4.1];



%Do calculations within trials
frac_pred_sp_baseline=[];
frac_pred_sp_baseline_sh=[];
frac_pred_sp_odor=[];
frac_pred_sp_odor_sh=[];
frac_pred_sm_baseline=[];
frac_pred_sm_baseline_sh=[];
frac_pred_sm_odor=[];
frac_pred_sm_odor_sh=[];



%S+
ii_trial_inc=0;
for ii_trial=1:2:length(sp_times)

    these_times_base=(time>=sp_times(ii_trial)+baseline_t(1))&(time<=sp_times(ii_trial)+baseline_t(2));
    these_times=(time>=sp_times(ii_trial)+odor_t(1))&(time<=sp_times(ii_trial)+odor_t(2));

    if (sum(these_times_base~=0))&(sum(these_times~=0))
        ii_trial_inc=ii_trial_inc+1;

        %Baseline
        frac_pred_sp_baseline(ii_trial_inc)=sum(label_traces(these_times_base)==1)/sum(these_times_base);
        frac_pred_sp_baseline_sh(ii_trial_inc)=sum(label_traces_sh(these_times_base)==1)/sum(these_times_base);

        %Odor
        frac_pred_sp_odor(ii_trial_inc)=sum(label_traces(these_times)==1)/sum(these_times);
        frac_pred_sp_odor_sh(ii_trial_inc)=sum(label_traces_sh(these_times)==1)/sum(these_times);
    end
end

%S-
ii_trial_inc=0;
for ii_trial=1:2:length(sm_times)
    these_times_base=(time>=sm_times(ii_trial)+baseline_t(1))&(time<=sm_times(ii_trial)+baseline_t(2));
    these_times=(time>=sm_times(ii_trial)+odor_t(1))&(time<=sm_times(ii_trial)+odor_t(2));

    if (sum(these_times_base~=0))&(sum(these_times~=0))
        ii_trial_inc=ii_trial_inc+1;

        %Baseline
        frac_pred_sm_baseline(ii_trial_inc)=sum(label_traces(these_times_base)==1)/sum(these_times_base);
        frac_pred_sm_baseline_sh(ii_trial_inc)=sum(label_traces_sh(these_times_base)==1)/sum(these_times_base);

        %Odor
        frac_pred_sm_odor(ii_trial_inc)=sum(label_traces(these_times)==1)/sum(these_times);
        frac_pred_sm_odor_sh(ii_trial_inc)=sum(label_traces_sh(these_times)==1)/sum(these_times);
    end
end


%and calculation for spontaneous and no change between trials
frac_pred_spon_baseline=[];
frac_pred_spon_baseline_sh=[];
frac_pred_spon_onset=[];
frac_pred_spon_onset_sh=[];
frac_pred_nc_baseline=[];
frac_pred_nc_baseline_sh=[];
frac_pred_nc_onset=[];
frac_pred_nc_onset_sh=[];



%Spontaneous increase in prediction between trials
ii_trial_inc=0;
for ii_trial=1:length(spontaneous_trial_times)

    these_times_base=(time>=spontaneous_trial_times(ii_trial)-2)&(time<=spontaneous_trial_times(ii_trial)-0.5);
    these_times=(time>=spontaneous_trial_times(ii_trial))&(time<=spontaneous_trial_times(ii_trial)+0.2);

    if (sum(these_times_base~=0))&(sum(these_times~=0))
        ii_trial_inc=ii_trial_inc+1;

        %Baseline
        frac_pred_spon_baseline(ii_trial_inc)=sum(label_traces(these_times_base)==1)/sum(these_times_base);
        frac_pred_spon_baseline_sh(ii_trial_inc)=sum(label_traces_sh(these_times_base)==1)/sum(these_times_base);

        %Odor
        frac_pred_spon_onset(ii_trial_inc)=sum(label_traces(these_times)==1)/sum(these_times);
        frac_pred_spon_onset_sh(ii_trial_inc)=sum(label_traces_sh(these_times)==1)/sum(these_times);
    end
end

%No change in prediction between trials
ii_trial_inc=0;
for ii_trial=1:2:length(no_increase_trial_times)
    these_times_base=(time>=no_increase_trial_times(ii_trial)-2)&(time<=no_increase_trial_times(ii_trial)-0.5);
    these_times=(time>=no_increase_trial_times(ii_trial))&(time<=no_increase_trial_times(ii_trial)+0.2);

    if (sum(these_times_base~=0))&(sum(these_times~=0))
        ii_trial_inc=ii_trial_inc+1;

        %Baseline
        frac_pred_nc_baseline(ii_trial_inc)=sum(label_traces(these_times_base)==1)/sum(these_times_base);
        frac_pred_nc_baseline_sh(ii_trial_inc)=sum(label_traces_sh(these_times_base)==1)/sum(these_times_base);

        %Odor
        frac_pred_nc_onset(ii_trial_inc)=sum(label_traces(these_times)==1)/sum(these_times);
        frac_pred_nc_onset_sh(ii_trial_inc)=sum(label_traces_sh(these_times)==1)/sum(these_times);
    end
end

%and calculations for the entire interval between trials
frac_pred_between=[];
frac_pred_between_sh=[];

for ii_between=1:length(between_start)

    these_times=(time>=between_start(ii_between))&(time<=between_end(ii_between));
    
    if (sum(these_times~=0))
        frac_pred_between(ii_between)=sum(label_traces(these_times)==1)/sum(these_times);
        frac_pred_between_sh(ii_between)=sum(label_traces_sh(these_times)==1)/sum(these_times);
    end

end

if show_figures==1
    edges=[0:0.05:1];
    rand_offset=0.8;

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    bar_offset=0;

    %S+
    bar(bar_offset,mean(frac_pred_sp_baseline_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sp_baseline_sh...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_sp_baseline),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 0/255 0/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sp_baseline...
        ,edges,bar_offset,rand_offset,'k','k',3);

    for ii=1:length(frac_pred_sp_baseline_sh)
        plot([bar_offset-1 bar_offset],[frac_pred_sp_baseline_sh(ii) frac_pred_sp_baseline(ii)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+2;

    bar(bar_offset,mean(frac_pred_sp_odor_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sp_odor_sh...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_sp_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 0/255 0/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sp_odor...
        ,edges,bar_offset,rand_offset,'k','k',3);

    for ii=1:length(frac_pred_sp_odor_sh)
        plot([bar_offset-1 bar_offset],[frac_pred_sp_odor_sh(ii) frac_pred_sp_odor(ii)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+3;

    %S-
    bar(bar_offset,mean(frac_pred_sm_baseline_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sm_baseline_sh...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_sm_baseline),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 0/255 0/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sm_baseline...
        ,edges,bar_offset,rand_offset,'k','k',3);

    for ii=1:length(frac_pred_sm_baseline_sh)
        plot([bar_offset-1 bar_offset],[frac_pred_sm_baseline_sh(ii) frac_pred_sm_baseline(ii)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+2;

    bar(bar_offset,mean(frac_pred_sm_odor_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sm_odor_sh...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_sm_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 0/255 0/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sm_odor...
        ,edges,bar_offset,rand_offset,'k','k',3);

    for ii=1:length(frac_pred_sm_odor_sh)
        plot([bar_offset-1 bar_offset],[frac_pred_sm_odor_sh(ii) frac_pred_sm_odor(ii)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+3;

    %between
    bar(bar_offset,mean(frac_pred_between_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_between_sh...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_between),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 0/255 0/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_between...
        ,edges,bar_offset,rand_offset,'k','k',3);

    for ii=1:length(frac_pred_between_sh)
        plot([bar_offset-1 bar_offset],[frac_pred_between_sh(ii) frac_pred_between(ii)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end



    xticks([0.5 3.5 7.5 10.5 14.5])
    xticklabels({'base S+', 'odor S+', 'base S-', 'odor S-','between'})

    ylabel('Fraction S+ prediction')



    title(['Fraction of S+ predictions'])

end

handles_out3.frac_pred_between=frac_pred_between;
handles_out3.frac_pred_between_sh=frac_pred_between_sh;
handles_out3.frac_pred_sm_odor=frac_pred_sm_odor;
handles_out3.frac_pred_sm_odor_sh=frac_pred_sm_odor_sh;
handles_out3.frac_pred_sp_odor=frac_pred_sp_odor;
handles_out3.frac_pred_sp_odor_sh=frac_pred_sp_odor_sh;
handles_out3.frac_pred_sm_baseline=frac_pred_sm_baseline;
handles_out3.frac_pred_sm_baseline_sh=frac_pred_sm_baseline_sh;
handles_out3.frac_pred_sp_baseline=frac_pred_sp_baseline;
handles_out3.frac_pred_sp_baseline_sh=frac_pred_sp_baseline_sh;
handles_out3.before_spontaneous_spm=before_spontaneous_spm;
handles_out3.after_spontaneous_spm=after_spontaneous_spm;
handles_out3.frac_pred_spon_baseline=frac_pred_spon_baseline;
handles_out3.frac_pred_spon_baseline_sh=frac_pred_spon_baseline_sh;
handles_out3.frac_pred_spon_onset=frac_pred_spon_onset;
handles_out3.frac_pred_spon_onset_sh=frac_pred_spon_onset_sh;
handles_out3.frac_pred_nc_baseline=frac_pred_nc_baseline;
handles_out3.frac_pred_nc_baseline_sh=frac_pred_nc_baseline_sh;
handles_out3.frac_pred_nc_onset=frac_pred_nc_onset;
handles_out3.frac_pred_nc_onset_sh=frac_pred_nc_onset_sh;


%Now find the correlation between the spontaneous dFFs and the dFFs per
%trial for S+ and S-
corr_sp=[];
corr_sm=[];
for ii_spont=1:length(spontaneous_trial_times)
    ii_time_spont=find(time>=spontaneous_trial_times(ii_spont),1,'first');

    spont_dFF=zeros(1,size(traces,1));
    spont_dFF(1,:)=traces(:,ii_time_spont);

    %Now find r for each time point for S+
    ii_sp_included=0;
    these_corr_sp=[];
    for sp_ii=1:2:length(sp_times)

        ii_t_first=find(time>=sp_times(sp_ii)-7,1,'first');
        ii_t_last=find(time<=sp_times(sp_ii)+15,1,'last');
        if (~isempty(ii_t_first))&(~isempty(ii_t_last))
            ii_sp_included=ii_sp_included+1;

            for ii_t=ii_t_first:ii_t_last
                this_dFF=zeros(1,size(traces,1));
                this_dFF(1,:)=traces(:,ii_t);
                this_cc=corrcoef(spont_dFF+0.0001,this_dFF+0.0001);
                these_corr_sp(ii_sp_included,ii_t-ii_t_first+1)=this_cc(1,2);
            end

        end
    end
    corr_sp(ii_spont,:)=mean(these_corr_sp);

    %Now find r for each time point for S-
    ii_sm_included=0;
    these_corr_sm=[];
    for sm_ii=1:2:length(sm_times)

        ii_t_first=find(time>=sm_times(sm_ii)-7,1,'first');
        ii_t_last=find(time<=sm_times(sm_ii)+15,1,'last');
        if (~isempty(ii_t_first))&(~isempty(ii_t_last))
            ii_sm_included=ii_sm_included+1;

            for ii_t=ii_t_first:ii_t_last
                this_dFF=zeros(1,size(traces,1));
                this_dFF(1,:)=traces(:,ii_t);
                this_cc=corrcoef(spont_dFF+0.0001,this_dFF+0.0001);
                these_corr_sm(ii_sm_included,ii_t-ii_t_first+1)=this_cc(1,2);
            end

        end
    end
    corr_sm(ii_spont,:)=mean(these_corr_sm);



end

this_time_span=([1:size(corr_sp,2)]-1)*(mean(time(2:end)-time(1:end-1)))-7;

if show_figures==1
    %Plot the correlation for spontaneous S+ prediction vs. S+ and S-
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .4 .4])

    hold on

    CIsm = bootci(1000, @mean, corr_sm);
    meansm=mean(corr_sm,1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(this_time_span',mean(corr_sm,1), CIsm', 'cmap',[158/255 31/255 99/255]);

    CIsp = bootci(1000, @mean, corr_sp);
    meansp=mean(corr_sp,1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;


    [hlsp, hpsp] = boundedline(this_time_span',mean(corr_sp,1), CIsp', 'cmap',[0 114/255 178/255]);

    plot(this_time_span',mean(corr_sm,1),'-','Color',[158/255 31/255 99/255],'DisplayName','S-')
    plot(this_time_span',mean(corr_sp,1), '-','Color',[0 114/255 178/255],'DisplayName','S+');


    ylim([-0.05 0.3])
    this_ylim=ylim;

    text(12,0.85*this_ylim(2) ,'S-','Color',[158/255 31/255 99/255])
    text(12,0.8*this_ylim(2),'S+','Color',[0 114/255 178/255])

    %FV
    rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor
    rectangle(Position=[0,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(delta_odor),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0 0 0], EdgeColor=[0 0 0])
    plot([0 0],[this_ylim],'-k')
    plot([mean(delta_odor) mean(delta_odor)],[this_ylim],'-k')

    %Reinforcement
    rectangle(Position=[mean(delta_odor_on_reinf_on),this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(delta_reinf),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[1 0 0], EdgeColor=[1 0 0])
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[this_ylim],'-r')
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[this_ylim],'-r')



    %     ylim([0 1.1])
    xlim([-7 15])
    title(['Correlation for spontaneous vs. actual dFF ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Rho')
end

handles_out3.this_time_span_corr=this_time_span;
handles_out3.corr_sp=corr_sp;
handles_out3.corr_sm=corr_sm;

%Now do the correlation for the no change traces
nc_corr_sp=[];
nc_corr_sm=[];
for ii_no_inc=1:length(no_increase_trial_times)
    ii_time_spont=find(time>=no_increase_trial_times(ii_no_inc),1,'first');

    if ~isempty(ii_time_spont)
        spont_dFF=zeros(1,size(traces,1));
        spont_dFF(1,:)=traces(:,ii_time_spont);

        %Now find r for each time point for S+
        ii_sp_included=0;
        these_nc_corr_sp=[];
        for sp_ii=1:2:length(sp_times)

            ii_t_first=find(time>=sp_times(sp_ii)-7,1,'first');
            ii_t_last=find(time<=sp_times(sp_ii)+15,1,'last');
            if (~isempty(ii_t_first))&(~isempty(ii_t_last))
                ii_sp_included=ii_sp_included+1;

                for ii_t=ii_t_first:ii_t_last
                    this_dFF=zeros(1,size(traces,1));
                    this_dFF(1,:)=traces(:,ii_t);
                    this_cc=corrcoef(spont_dFF+0.0001,this_dFF+0.0001);
                    these_nc_corr_sp(ii_sp_included,ii_t-ii_t_first+1)=this_cc(1,2);
                end

            end
        end
        nc_corr_sp(ii_no_inc,:)=mean(these_nc_corr_sp);

        %Now find r for each time point for S-
        ii_sm_included=0;
        these_nc_corr_sm=[];
        for sm_ii=1:2:length(sm_times)

            ii_t_first=find(time>=sm_times(sm_ii)-7,1,'first');
            ii_t_last=find(time<=sm_times(sm_ii)+15,1,'last');
            if (~isempty(ii_t_first))&(~isempty(ii_t_last))
                ii_sm_included=ii_sm_included+1;

                for ii_t=ii_t_first:ii_t_last
                    this_dFF=zeros(1,size(traces,1));
                    this_dFF(1,:)=traces(:,ii_t);
                    this_cc=corrcoef(spont_dFF+0.0001,this_dFF+0.0001);
                    these_nc_corr_sm(ii_sm_included,ii_t-ii_t_first+1)=this_cc(1,2);
                end

            end
        end
        nc_corr_sm(ii_no_inc,:)=mean(these_nc_corr_sm);

    end

end

this_time_span=([1:size(nc_corr_sp,2)]-1)*(mean(time(2:end)-time(1:end-1)))-7;

if show_figures==1
    %Plot the correlation for spontaneous S+ prediction vs. S+ and S-
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .4 .4])

    hold on

    if ~isempty(nc_corr_sm)
        if size(nc_corr_sm,1)>2
            CIsm = bootci(1000, @mean, nc_corr_sm);
            meansm=mean(nc_corr_sm,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(this_time_span',mean(nc_corr_sm,1), CIsm', 'cmap',[158/255 31/255 99/255]);
        end
    end

    if ~isempty(nc_corr_sp)
        if size(nc_corr_sp,1)>2
            CIsp = bootci(1000, @mean, nc_corr_sp);
            meansp=mean(nc_corr_sp,1);
            CIsp(1,:)=meansp-CIsp(1,:);
            CIsp(2,:)=CIsp(2,:)-meansp;


            [hlsp, hpsp] = boundedline(this_time_span',mean(nc_corr_sp,1), CIsp', 'cmap',[0 114/255 178/255]);
        end
    end

    if ~isempty(nc_corr_sm)
        plot(this_time_span',mean(nc_corr_sm,1),'-','Color',[158/255 31/255 99/255],'DisplayName','S-')
    end

    if ~isempty(nc_corr_sp)
        plot(this_time_span',mean(nc_corr_sp,1), '-','Color',[0 114/255 178/255],'DisplayName','S+');
    end


    ylim([-0.05 0.3])
    this_ylim=ylim;

    text(12,0.85*this_ylim(2) ,'S-','Color',[158/255 31/255 99/255])
    text(12,0.8*this_ylim(2),'S+','Color',[0 114/255 178/255])

    %FV
    rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor
    rectangle(Position=[0,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(delta_odor),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0 0 0], EdgeColor=[0 0 0])
    plot([0 0],[this_ylim],'-k')
    plot([mean(delta_odor) mean(delta_odor)],[this_ylim],'-k')

    %Reinforcement
    rectangle(Position=[mean(delta_odor_on_reinf_on),this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(delta_reinf),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[1 0 0], EdgeColor=[1 0 0])
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[this_ylim],'-r')
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[this_ylim],'-r')



    %     ylim([0 1.1])
    xlim([-7 15])
    title(['Correlation for no change between trials vs. actual dFF ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Rho')
end


handles_out3.this_time_span_nc_corr=this_time_span;
handles_out3.nc_corr_sp=nc_corr_sp;
handles_out3.nc_corr_sm=nc_corr_sm;



%Plot the traces showing the spontaneous times
%time has the time for the dF/F traces(ROI,time)
if show_figures==1

    dt_p_threshold=4.2; %Time to be used after the odor on for the p_threshold t_test
     post_time=5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on
    post_shift=0; %Set to 0 if you want to train with odor on points
    pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
%     MLalgo_to_use=[1]; %Vector with the decoding algorithms you want to use
%     ii_cost=3;
    p_threshold=1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold

    %Post points
    Nall=size(traces,1);
    dt=median(time(2:end)-time(1:end-1));
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
    dt_span=40;
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

                next_ii=next_ii_sp;
                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)&(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post,:)=1;
                    ii_sp_post=ii_sp_post+1;
                    dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                    dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;

                    %Save trial classification
                    if epochs(this_ii+1+next_ii_sp)==6
                        hit_per_trial=[hit_per_trial 1];
                        miss_per_trial=[miss_per_trial 0];
                    else
                        hit_per_trial=[hit_per_trial 0];
                        miss_per_trial=[miss_per_trial 1];
                    end
                    cr_per_trial=[cr_per_trial 0];
                    fa_per_trial=[fa_per_trial 0];

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


                next_ii=next_ii_sm;
                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post,:)=0;
                    ii_sm_post=ii_sm_post+1;
                    dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                    dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;

                    %Save trial classification
                    if epochs(this_ii+1+next_ii_sm)==9
                        cr_per_trial=[cr_per_trial 1];
                        fa_per_trial=[fa_per_trial 0];
                    else
                        cr_per_trial=[cr_per_trial 0];
                        fa_per_trial=[fa_per_trial 1];
                    end
                    hit_per_trial=[hit_per_trial 0];
                    miss_per_trial=[miss_per_trial 0];

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



                    if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)...
                            &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                        measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                        training_decisions_post(ii_post+1:ii_post+no_points_post,:)=0;
                        ii_sm_post=ii_sm_post+1;
                        dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                        dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                        dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;

                        %Save trial classification
                        if epochs(this_ii+1+next_ii_sm)==9
                            cr_per_trial=[cr_per_trial 1];
                            fa_per_trial=[fa_per_trial 0];
                        else
                            cr_per_trial=[cr_per_trial 0];
                            fa_per_trial=[fa_per_trial 1];
                        end
                        hit_per_trial=[hit_per_trial 0];
                        miss_per_trial=[miss_per_trial 0];

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



                    if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)...
                            &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)
                        measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                        training_decisions_post(ii_post+1:ii_post+no_points_post,:)=1;
                        ii_sp_post=ii_sp_post+1;
                        dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                        dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                        dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;

                        %Save trial classification
                        if epochs(this_ii+1+next_ii_sp)==6
                            hit_per_trial=[hit_per_trial 1];
                            miss_per_trial=[miss_per_trial 0];
                        else
                            hit_per_trial=[hit_per_trial 0];
                            miss_per_trial=[miss_per_trial 1];
                        end
                        cr_per_trial=[cr_per_trial 0];
                        fa_per_trial=[fa_per_trial 0];

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
    p_values=ones(1,size(dFFs_sm_per_trial_per_ROI,2));
    for iiROI=1:size(dFFs_sm_per_trial_per_ROI,2)
        dFF_sm=zeros(size(dFFs_sm_per_trial_per_ROI,1),size(dFFs_sm_per_trial_per_ROI,3));
        dFF_sm(:,:)=dFFs_sm_per_trial_per_ROI(:,iiROI,:);
        dFF_sp=zeros(size(dFFs_sp_per_trial_per_ROI,1),size(dFFs_sp_per_trial_per_ROI,3));
        dFF_sp(:,:)=dFFs_sp_per_trial_per_ROI(:,iiROI,:);

        [h,p_values(iiROI)]=ttest2(mean(dFF_sp,2),mean(dFF_sm,2));
    end

    p_value_mask=logical(p_values<=p_threshold);

    pFDR=drsFDRpval(p_values);
    fprintf(1, ['pFDR =%d, among %d ROIs, %d ROIs significant\n'],pFDR,length(p_values),sum(p_values<=pFDR));
    %Plot the trimmed traces sorted according to p values
    to_sort=zeros(no_traces,2);
    to_sort(:,1)=1:no_traces;
    to_sort(:,2)=p_values;
    sorted=sortrows(to_sort,2);
    traces_sorted=zeros(no_traces,1);
    traces_sorted(:,1)=sorted(:,1);

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
                    'Color',[80/255 194/255 255/255],'LineWidth',1.5)
            else
                plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                    'Color',[238/255 111/255 179/255],'LineWidth',1.5)
            end
        end
    end

    %Plot spontaneous times
    for ii_spont=1:length(spontaneous_trial_times)
        plot([spontaneous_trial_times(ii_spont) spontaneous_trial_times(ii_spont)], [0 (no_traces+2)*y_shift],...
            'Color',[150/255 150/255 150/255],'LineWidth',1.5)

    end


    for trNo=1:no_traces
        % for trNo=1:20
        plot(decimate(time,5),decimate(traces(traces_sorted(trNo),:),5)+y_shift*trNo,'-k','LineWidth',1)
    end

    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    xlabel('time(sec)')
    title(['dFF timecourses after p value sorting ' num2str(size(measurements_post,2)) ' ROIs'])
end

pffft=1;