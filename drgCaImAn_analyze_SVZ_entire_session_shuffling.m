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
    iiMLalgo=6; %Vector with the decoding algorithms you want to use
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
fprintf(1, ['\ndrgCaImAn_analyze_SVZ_entire_session_shuffling run for ' pre_perFileName '\n\n']);
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


delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

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
fprintf(1, ['Percent correct for behavior  is %d\n'],...
    handles_out.ii_out(1).handles_out.percent_correct);


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

    

    odor_window=[2 4.1]; %Note: This is the time window I use in drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse_v2
    odor_acc=mean(meansm((time_span>=odor_window(1))&(time_span<=odor_window(2))));
    fprintf(1, ['Accuracy for odor window  is %d\n'],...
    odor_acc);

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
    
    %Now find the changes in prediction that take place between trials
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

    for ii_trial=1:2:length(all_trial_times)-2
        %find the next trial
        this_ii=find(time>all_trial_times(ii_trial)+delta_odor,1,'first');
        if ~isempty(this_ii)
            %Now find the first baseline period
            delta_ii=find(moving_mean_label_traces(this_ii:end) < per5,1,'first');
            if ~isempty(delta_ii)
                if time(this_ii+delta_ii-1)<all_trial_times(ii_trial+2)
                    %Find the increase
                    delta_ii2=find(moving_mean_label_traces(this_ii+delta_ii-1:end) > per95,1,'first');
                    
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

                        spontaneous_delta_prediction(ii_spontaneous,:)=moving_mean_label_traces((time>=-7+spontaneous_trial_times(ii_spontaneous))&(time<=15+spontaneous_trial_times(ii_spontaneous)));

                        pffft=1;
                    end
                end
            end
        end
    end

    %Plot label predictions alone
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

    %Plot the spontaneous between trial dFF
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

    dt=time(2)-time(1);
    this_time_span=[-7:dt:15];
    this_time_span=this_time_span(1:size(spontaneous_delta_prediction,2));
    
   
    hold on
    
    CIsm = bootci(1000, @mean, spontaneous_delta_prediction);
    meansm=mean(spontaneous_delta_prediction,1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;
    
    [hlsm, hpsm] = boundedline(this_time_span',mean(spontaneous_delta_prediction,1)', CIsm', 'cmap',[0/255 0/255 0/255]);
    
    for ii=1:size(spontaneous_delta_prediction,1)
        plot(this_time_span',spontaneous_delta_prediction(ii,:),  'Color',[150/255 150/255 150/255])
    end
    
    
    text(30,0.75,'S-','Color',[158/255 31/255 99/255])
    text(30,0.85,'S+','Color',[0 114/255 178/255])
    
    ylim([0 1.1])
    xlim([-7 15])
    title(['Label prediction spontaneous between trials for ' classifier_names{iiMLalgo} ' and p threshold ' num2str(p_threshold)])
    xlabel('Time(sec)')
    ylabel('Label prediction')

 

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

    %Histogram for all points
     figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    hold on
    
    edges=[0:0.1:1];

    histogram(moving_mean_label_traces,edges)
    hold on 
    histogram(moving_mean_label_traces_sh,edges)

    
    title(['Histogram for labels ' ])
    
    
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