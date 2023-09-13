function handles_out=drgCaImAn_analyze_one_decode_licks_entire_sessionv2(handles_choices)
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
    [pre_perFileName,pre_perPathName] = uigetfile({'*.mat'},'Select the .mat file for analysis');


else

end

load([pre_perPathName pre_perFileName])

 
classifier_names{1}='LDA';
classifier_names{2}='SVM';
classifier_names{3}='Bayes';
classifier_names{4}='ANN';
classifier_names{5}='Tree';
classifier_names{6}='GLM';

moving_mean_n=10;

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

figNo=0;
  
time_span=handles_out.ii_out(1).handles_out.time_span;
training_time_from=handles_out.handles.training_time_from;
training_delta_time=handles_out.handles.training_delta_time;
length_trimmed_time_span=sum((time_span>=training_time_from(1))&(time_span<=training_time_from(end)+training_delta_time));
trimmed_time_span=time_span((time_span>=training_time_from(1))&(time_span<=training_time_from(end)+training_delta_time));



ii_cp=0;
ii_step=size(handles_out.ii_out(1).handles_out.correct_predict_post,2)-2;
n_shuffled=size(handles_out.ii_out(1).handles_out.correct_predict_sh_post,1)/size(handles_out.ii_out(1).handles_out.correct_predict_post,1);
no_trials=size(handles_out.ii_out(1).handles_out.correct_predict_post,1);


%Note that not all trials are included in each of these runs, we have to
%use the smallest number of trials
no_trials_per_time_from=zeros(1,length(training_time_from));
for ii=1:length(training_time_from)
    no_trials_per_time_from(ii)=size(handles_out.ii_out(ii).handles_out.correct_predict_post,1);
end

u_no_trials=unique(no_trials_per_time_from);


ii_standard=find(no_trials_per_time_from==min(u_no_trials),1,'first');
st_no_trials=no_trials_per_time_from(ii_standard);
st_hit_per_trial=handles_out.ii_out(ii_standard).handles_out.hit_per_trial;
st_miss_per_trial=handles_out.ii_out(ii_standard).handles_out.miss_per_trial;
st_cr_per_trial=handles_out.ii_out(ii_standard).handles_out.cr_per_trial;
st_fa_per_trial=handles_out.ii_out(ii_standard).handles_out.fa_per_trial;
st_sp_no_trials=sum(st_hit_per_trial)+sum(st_miss_per_trial);
st_sm_no_trials=sum(st_cr_per_trial)+sum(st_fa_per_trial);

trials_from=zeros(1,length(training_time_from));
trials_to=zeros(1,length(training_time_from));
for ii=1:length(training_time_from)
    if no_trials_per_time_from(ii_standard)==no_trials_per_time_from(ii)
        trials_from(ii)=1;
        trials_to(ii)=min(u_no_trials);
    else
        for jj=0:no_trials_per_time_from(ii)-min(u_no_trials)

            this_hit_per_trial=handles_out.ii_out(ii).handles_out.hit_per_trial(1+jj:st_no_trials+jj);
            sumh=0;
            for iih=1:min(u_no_trials)
                if this_hit_per_trial(iih)==st_hit_per_trial(iih)
                   sumh=sumh+1; 
                end
            end
            
            this_miss_per_trial=handles_out.ii_out(ii).handles_out.miss_per_trial(1+jj:st_no_trials+jj);
            summ=0;
            for iim=1:min(u_no_trials)
                if this_miss_per_trial(iim)==st_miss_per_trial(iim)
                   summ=summ+1; 
                end
            end
            
            this_cr_per_trial=handles_out.ii_out(ii).handles_out.cr_per_trial(1+jj:st_no_trials+jj);
             sumc=0;
            for iic=1:min(u_no_trials)
                if this_cr_per_trial(iim)==st_cr_per_trial(iim)
                   sumc=sumc+1; 
                end
            end
            
            this_fa_per_trial=handles_out.ii_out(ii).handles_out.fa_per_trial(1+jj:st_no_trials+jj);
             sumf=0;
            for iif=1:min(u_no_trials)
                if this_fa_per_trial(iif)==st_fa_per_trial(iif)
                   sumf=sumf+1; 
                end
            end
            
            if (sumh==min(u_no_trials))&(summ==min(u_no_trials))&(sumc==min(u_no_trials))&(sumf==min(u_no_trials))
                trials_from(ii)=1+jj;
                trials_to(ii)=st_no_trials+jj;
            end
        end
    end
    
end

%Now assign correct predict, lick fraction, etc
this_correct_predict=zeros(min(u_no_trials),length_trimmed_time_span);
this_correct_predict_sh=zeros(min(u_no_trials),length_trimmed_time_span);



for ii=1:length(training_time_from)
    
    %Do correct predict
    this_correct_predict(:,ii_cp+1:ii_cp+ii_step)=handles_out.ii_out(ii).handles_out.correct_predict_post(trials_from(ii):trials_to(ii),2:end-1);
    
    these_correct_predict_sh=zeros(min(u_no_trials),ii_step);  
    for jj_step=1:ii_step
        jj_cp=0;
        these_cp_sh=zeros(n_shuffled,min(u_no_trials));
        for jj=1:n_shuffled
            these_cp_sh(jj,:)=handles_out.ii_out(ii).handles_out.correct_predict_sh_post(jj_cp+trials_from(ii):jj_cp+trials_to(ii),jj_step+1);
            jj_cp=jj_cp+ii_step;
        end
        these_correct_predict_sh(:,jj_step)=mean(these_cp_sh);
    end
    this_correct_predict_sh(:,ii_cp+1:ii_cp+ii_step)=these_correct_predict_sh;
    
    %Do licks
    trials_from_sp=sum(handles_out.ii_out(ii).handles_out.hit_per_trial(1:trials_from(ii)))+sum(handles_out.ii_out(ii).handles_out.miss_per_trial(1:trials_from(ii)));
    if trials_from_sp==0
        trials_from_sp=1;
    end
    trials_to_sp=sum(handles_out.ii_out(ii).handles_out.hit_per_trial(1:trials_to(ii)))+sum(handles_out.ii_out(ii).handles_out.miss_per_trial(1:trials_to(ii)));
    
    trials_from_sm=sum(handles_out.ii_out(ii).handles_out.cr_per_trial(1:trials_from(ii)))+sum(handles_out.ii_out(ii).handles_out.fa_per_trial(1:trials_from(ii)));
    if trials_from_sm==0
        trials_from_sm=1;
    end
    trials_to_sm=sum(handles_out.ii_out(ii).handles_out.cr_per_trial(1:trials_to(ii)))+sum(handles_out.ii_out(ii).handles_out.fa_per_trial(1:trials_to(ii)));
    
    if ii==1
        this_per_trial_sp_licks_post=zeros(trials_to_sp-trials_from_sp+1,length_trimmed_time_span);
        this_per_trial_sm_licks_post=zeros(trials_to_sm-trials_from_sm+1,length_trimmed_time_span);
        
        this_prediction_per_trial_sp_timecourse=zeros(trials_to_sp-trials_from_sp+1,length_trimmed_time_span);
        this_prediction_per_trial_sm_timecourse=zeros(trials_to_sm-trials_from_sm+1,length_trimmed_time_span);
    end

    this_per_trial_sp_licks_post(:,ii_cp+1:ii_cp+ii_step)=handles_out.ii_out(ii).handles_out.per_trial_sp_licks_post(trials_from_sp:trials_to_sp,2:end-1);
    this_per_trial_sm_licks_post(:,ii_cp+1:ii_cp+ii_step)=handles_out.ii_out(ii).handles_out.per_trial_sm_licks_post(trials_from_sm:trials_to_sm,2:end-1);
    
    %Do predictions
    this_prediction_per_trial_sp_timecourse(:,ii_cp+1:ii_cp+ii_step) = handles_out.ii_out(ii).handles_out.per_trial_sp_timecourse_post(trials_from_sp:trials_to_sp,2:end-1);
    this_prediction_per_trial_sm_timecourse(:,ii_cp+1:ii_cp+ii_step) = handles_out.ii_out(ii).handles_out.per_trial_sm_timecourse_post(trials_from_sm:trials_to_sm,2:end-1);
    
    ii_cp=ii_cp+ii_step;
end

this_per_trial_mean_sp_licks_post = movmean(this_per_trial_sp_licks_post',moving_mean_n)';
this_per_trial_mean_sm_licks_post = movmean(this_per_trial_sm_licks_post',moving_mean_n)';


this_prediction_mean_per_trial_sp_timecourse = movmean(this_prediction_per_trial_sp_timecourse',moving_mean_n)';
this_prediction_mean_per_trial_sm_timecourse = movmean(this_prediction_per_trial_sm_timecourse',moving_mean_n)';

%Plot accuracy
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

[hlsm, hpsm] = boundedline(trimmed_time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');


CIsm = bootci(1000, @mean, this_correct_predict);
meansm=mean(this_correct_predict,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(trimmed_time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);

plot(trimmed_time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
plot(trimmed_time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);

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


title(['Accuracy v3  trained per time point' ])
xlabel('Time(sec)')
ylabel('Accuracy')

%Plor lick fraction
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on


CIsm = bootci(1000, @mean, this_per_trial_mean_sm_licks_post);
meansm=mean(this_per_trial_mean_sm_licks_post,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(trimmed_time_span',mean(this_per_trial_mean_sm_licks_post,1)', CIsm', 'cmap',[158/255 31/255 99/255]);


CIsm = bootci(1000, @mean, this_per_trial_mean_sp_licks_post);
meansm=mean(this_per_trial_mean_sp_licks_post,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(trimmed_time_span',mean(this_per_trial_mean_sp_licks_post,1)', CIsm', 'cmap',[0 114/255 178/255]);

plot(trimmed_time_span',mean(this_per_trial_mean_sm_licks_post,1)','Color',[158/255 31/255 99/255])
plot(trimmed_time_span',mean(this_per_trial_mean_sp_licks_post,1)', '-','Color',[0 114/255 178/255]);



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


title(['Lick fraction' ])
xlabel('Time(sec)')
ylabel('Lick fraction')


%Plor prediction
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on


CIsm = bootci(1000, @mean, this_prediction_mean_per_trial_sm_timecourse);
meansm=mean(this_prediction_mean_per_trial_sm_timecourse,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(trimmed_time_span',mean(this_prediction_mean_per_trial_sm_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);


CIsm = bootci(1000, @mean, this_prediction_mean_per_trial_sp_timecourse);
meansm=mean(this_prediction_mean_per_trial_sp_timecourse,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(trimmed_time_span',mean(this_prediction_mean_per_trial_sp_timecourse,1)', CIsm', 'cmap',[0 114/255 178/255]);

plot(trimmed_time_span',mean(this_prediction_mean_per_trial_sm_timecourse,1)','Color',[158/255 31/255 99/255])
plot(trimmed_time_span',mean(this_prediction_mean_per_trial_sp_timecourse,1)', '-','Color',[0 114/255 178/255]);



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


title(['Prediction' ])
xlabel('Time(sec)')
ylabel('Prediction')

pffft=1;
