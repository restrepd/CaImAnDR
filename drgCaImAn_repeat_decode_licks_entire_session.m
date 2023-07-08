function handles_out=drgCaImAn_repeat_decode_licks_entire_session(handles_choices)
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
    post_time=0.5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on
    post_shift=3; %Set to 0 if you want to train with odor on points
    pre_time=3; %Used to calculate the decoding accuracy pre_time sec before post_shift
    MLalgo_to_use=[4]; %Vector with the decoding algorithms you want to use
    ii_cost=3;
    p_threshold=1.1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
    dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures
    perm_before=0; %Permute the labels before running the decoder
    no_repeats=5;
    
    handles_choices.pre_per_FileName=pre_perFileName;
    handles_choices.pre_per_PathName=pre_perPathName;
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
    handles_choices.perm_before=perm_before;
    handles_choices.no_repeats=no_repeats;
    
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
    no_repeats=handles_choices.no_repeats;
    ii_cost=handles_choices.ii_cost;
    if isfield(handles_choices,'perm_before')
        perm_before=handles_choices.perm_before;
    else
        perm_before=0;
    end
end

classifier_names{1}='LDA';
classifier_names{2}='SVM';
classifier_names{3}='Bayes';
classifier_names{4}='ANN';
classifier_names{5}='Tree';
classifier_names{6}='GLM';


delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

handles_out=[];
these_mean_correct_predict=[];
these_mean_correct_predict_sh=[];
these_per_trial_sp_timecourse=[];
these_per_trial_sm_timecourse=[];

for ii_repeat=1:no_repeats
    tic
    this_handles_out=[];
    this_handles_out=drgCaImAn_decode_licks_entire_session(handles_choices);
    if ii_repeat==1
        handles_out=this_handles_out;
    end

    this_mean_correct_predict=[];
    this_mean_correct_predict=mean(this_handles_out.MLalgo(MLalgo_to_use).this_correct_predict);
    these_mean_correct_predict=[these_mean_correct_predict; this_mean_correct_predict];

    this_mean_correct_predict_sh=[];
    this_mean_correct_predict_sh=mean(this_handles_out.MLalgo(MLalgo_to_use).this_correct_predict_sh);
    these_mean_correct_predict_sh=[these_mean_correct_predict_sh; this_mean_correct_predict_sh];

    this_per_trial_sp_timecourse=[];
    this_per_trial_sp_timecourse=mean(this_handles_out.MLalgo(MLalgo_to_use).per_trial_sp_timecourse);
    these_per_trial_sp_timecourse=[these_per_trial_sp_timecourse; this_per_trial_sp_timecourse];

    this_per_trial_sm_timecourse=[];
    this_per_trial_sm_timecourse=mean(this_handles_out.MLalgo(MLalgo_to_use).per_trial_sm_timecourse);
    these_per_trial_sm_timecourse=[these_per_trial_sm_timecourse; this_per_trial_sm_timecourse];

    fprintf(1,['Repeat %d done in %d seconds\n'],ii_repeat,toc)
    pffft=1;

end

handles_out.MLalgo(MLalgo_to_use).this_correct_predict=these_mean_correct_predict;
handles_out.MLalgo(MLalgo_to_use).this_correct_predict_sh=these_mean_correct_predict_sh;
handles_out.MLalgo(MLalgo_to_use).per_trial_sp_timecourse=these_per_trial_sp_timecourse;
handles_out.MLalgo(MLalgo_to_use).per_trial_sm_timecourse=these_per_trial_sm_timecourse;
 
if show_figures==1
    this_correct_predict=these_mean_correct_predict;
    this_correct_predict_sh=these_mean_correct_predict_sh;
    time_span=handles_out.time_span;
    figNo=0;
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


    title(['Accuracy v3 for ' classifier_names{MLalgo_to_use} ' trained from ' num2str(post_shift) ' to ' num2str(post_shift+post_time)])
    xlabel('Time(sec)')
    ylabel('Accuracy')
end

pffft=1;
