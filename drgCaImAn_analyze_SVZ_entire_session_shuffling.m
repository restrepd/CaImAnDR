function handles_out2=drgCaImAn_analyze_SVZ_entire_session_shuffling(handles_choices)
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
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per_dec2.mat'},'Select the .m file with all the choices for analysis');
    
%     processing_algorithm=3; %Use 3
%     k_fold=5; %Only used for processing_algorithm=2,
%     post_time=5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on
%     post_shift=0; %Set to 0 if you want to train with odor on points
%     pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
    iiMLalgo=1; %Vector with the decoding algorithms you want to use
%     ii_cost=3;
    p_threshold=1.1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
%     dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures
    

    
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

figNo=0;

% convert_z=1; %Convert dFF traces to z
% dt_span=40; %Seconds shown before and after odor on in the figures
% moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
% no_shuffles=10; %Number of shuffles for per trial shuffling

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_analyze_SVZ_entire_session run for ' pre_perFileName '\n\n']);
% fprintf(1, 'post_time = %d, p_threshold= %d, post_shift= %d, cost %d\n',post_time,p_threshold,post_shift, ii_cost);

% if convert_z==1
%     for trace_no=1:size(traces,1)
%         traces(trace_no,:)=traces(trace_no,:)/std(traces(trace_no,:));
%     end
% end

classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';
classifier_names{6}='Binomial glm';

%Plot the mean timecourse for S+ and S-, Euclidean distance and KL
%divergence
ii_thr=find(handles_out.handles.p_threshold==p_threshold);
meandFFsp=handles_out.ii_out(ii_thr).handles_out.meandFFsp;
meandFFsm=handles_out.ii_out(ii_thr).handles_out.meandFFsm;
time_span=handles_out.ii_out(ii_thr).handles_out.time_span;
dist_euclid=handles_out.ii_out(ii_thr).handles_out.dist_euclid-handles_out.ii_out(ii_thr).handles_out.dist_euclid_zero;
KLdivergence=handles_out.ii_out(ii_thr).handles_out.KLdivergence;


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

    %S+ mean dFF
    plot(time_span,meandFFsp, '-' ,'Color',[80/255 194/255 255/255]);
    
    %S- mean dFF
    plot(time_span,meandFFsm, '-' ,'Color',[238/255 111/255 179/255]);
    
    
    text(30,0.75,'S-','Color',[80/255 194/255 255/255])
    text(30,0.65,'S+','Color',[0 114/255 178/255])
    
    
    title(['Mean dFF'])
    xlabel('Time(sec)')
    ylabel('dFF')
    
    subplot(3,1,2)
    hold on
    
    %Euclidean distance
    plot(time_span,dist_euclid, 'Color',[238/255 111/255 179/255]);
    
    
    
    title(['Euclidean distance'])
    xlabel('Time(sec)')
    ylabel('Distance')
    
    subplot(3,1,3)
    hold on
    
    %KL divergence
    plot(time_span',KLdivergence', 'Color',[238/255 111/255 179/255]);
    
    
    
    title(['KL divergence'])
    xlabel('Time(sec)')
    ylabel('KL divergence')
    
    
end

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
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    subplot(2,1,1)
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
    
    plot(time_span',mean(moving_mean_per_trial_sm_timecourse,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-','LineWidth',1.5)
    plot(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+','LineWidth',1.5);
    
    text(30,0.75,'S-','Color',[158/255 31/255 99/255])
    text(30,0.85,'S+','Color',[0 114/255 178/255])
    
    ylim([0 1])
    title(['Label prediction per trial for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Label prediction, S+=1, S-=0')
    
    subplot(2,1,2)
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
    
    ylim([0 1])
    title(['Accuracy per trial for ' classifier_names{iiMLalgo}])
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
    
    text(30,0.75,'S-','Color',[158/255 31/255 99/255])
    text(30,0.65,'S+','Color',[0 114/255 178/255])
    
    ylim([0 1])
    title(['Posterior probability for S+ or S- prediction for S+ trials ' classifier_names{iiMLalgo} ' ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Posterior probability')
    
    subplot(2,1,2)
    hold on
    
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
    
    text(30,0.75,'S-','Color',[158/255 31/255 99/255])
    text(30,0.65,'S+','Color',[0 114/255 178/255])
    
    ylim([0 1])
    title(['Posterior probability for S+ or S- prediction for S- trials ' classifier_names{iiMLalgo} ' ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Posterior probability')
    
end

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
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r','LineWidth',2)
    end
    
     for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b','LineWidth',2)
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
    

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

     ax=gca;ax.LineWidth=3;
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
    plot(time,moving_mean_label_traces,'-k','LineWidth',0.5)

    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[-0.1 1.1],'Color',[0 114/255 178/255],'LineWidth',2)
    end

    for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[-0.1 1.1],'Color',[158/255 31/255 99/255],'LineWidth',2)
    end
    %                 plot(time,moving_mean_label_traces_sh(1,:),'-b')

    ylim([-0.2 1.2])
    text(50,1.1,'S+','Color',[0 114/255 178/255])
    text(150,1.1,'S-','Color',[158/255 31/255 99/255])
    title(['Label prediction for entire session for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    hold on
    
    bar_offset=1;
    bar(bar_offset,sum(label_traces==0),'k')
    
     bar_offset=bar_offset+1;
    bar(bar_offset,sum(label_traces_sh==0),'m')
    
       bar_offset=bar_offset+2;
    bar(bar_offset,sum(label_traces==1),'k')
    
     bar_offset=bar_offset+1;
    bar(bar_offset,sum(label_traces_sh==1),'m')
    
     xticks([1 2  4 5])
    xticklabels({'S-', 'S- sh', 'S+','S+ sh'})
    
    ylabel('Number of timepoints')
    
    %Now do a Chi Squared test
    [p, Q]= chi2test([sum(label_traces==1), sum(label_traces==0); sum(label_traces_sh==1), sum(label_traces_sh==0)]);
    
    title(['Histogram for labels, Chi squared p value = ' num2str(p)])
    
    
end

prior_p_sm=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).scores(:,1);
prior_p_sp=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).scores(:,2);

prior_p_sm_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).scores_sh(:,1);
prior_p_sp_sh=handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).scores_sh(:,2);

prior_p_sm_sh2=zeros(1,length(time));
prior_p_sp_sh2=zeros(1,length(time));
prior_p_sm_sh2(1,:)=mean(handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).scores_sh2(:,:,1),1);
prior_p_sp_sh2(1,:)=mean(handles_out.ii_out(ii_thr).handles_out.MLalgo(iiMLalgo).scores_sh2(:,:,2),1);

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
    plot(time,prior_p_sp,'-r','LineWidth',1)
    plot(time,prior_p_sm,'-b','LineWidth',1)
  
    
    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end
    
     for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
     end
    

    
    ylim([-0.2 1.2])
    xlabel('Time(sec)')
    ylabel('p')
    title('original data')
    
    subplot(3,1,2)
    hold on
    
    
 
    
    %                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
    plot(time,prior_p_sp_sh,'-r','LineWidth',1)
    plot(time,prior_p_sm_sh,'-b','LineWidth',1)
  
    
    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end
    
     for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
     end
    

    
    ylim([-0.2 1.2])
    xlabel('Time(sec)')
    ylabel('p')
    title('shuffled')
    
    subplot(3,1,3)
    hold on
    
    
 
    
    %                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
    plot(time,prior_p_sp_sh2,'-r','LineWidth',1)
    plot(time,prior_p_sm_sh2,'-b','LineWidth',1)
  
    
    for ii=1:length(sp_times)
        plot([sp_times(ii) sp_times(ii)],[0 1],'-r')
    end
    
     for ii=1:length(sm_times)
        plot([sm_times(ii) sm_times(ii)],[0 1],'-b')
     end
    

    
    ylim([-0.2 1.2])
    xlabel('Time(sec)')
    ylabel('p')
    title('shuffled 2')
    
    sgtitle(['Prior probability for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    
    
      
    [h,p]=vartest2(prior_p_sp,prior_p_sp_sh);
    fprintf(1, ['p value for variance test for prior p sp vs. prior p sp sh ' num2str(p) '\n']);
    [h,p]=vartest2(prior_p_sm,prior_p_sm_sh);
    fprintf(1, ['p value for variance test for prior p sm vs. prior p sm sh ' num2str(p) '\n']);
end

pffft=1;