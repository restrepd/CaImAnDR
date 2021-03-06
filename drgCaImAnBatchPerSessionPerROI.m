function drgCaImAnBatchPerSessionPerROI

% This function calculates the per ROI dFF changes for go-no go sessions
%
% The input is a series of CalmAn_batch_pre_per.mat files with the CaImAn
% data for dFF for each ROI, the data on epochs and licks. 
% Each training session includes several of these files
% The name and location of these files and some choice parameters are 
% entered in a drgCaImAnChoices file
% caimanhandles.caimandr_choices.start_reversal is the file number for the 
% start of a reversal
% Processing of the data is different if there is a reversal
%
%
% Needs a choices file such as drgCaImAnChoices_20180515_mmPVG02_Cerebellum.m
% Needs the CalmAn_batch_pre_per.mat files output files from drgCaImAn_batch_dropc.m

%
warning('off')

close all
clear all

tic
 
[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSessionPerROI run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

caimanhandles=handles;

%Read the files and calculate the dF/F in each window
num_odor_trials=0;
epochs_per_trial=[];
num_odor_trials_dFF=0;
files_per_trial=[];

all_lda_events=[]; 
all_lda_input_timecourse=[];

per_file_trials_included=[];
time_this_file=[];
 
for filNum=1:caimanhandles.caimandr_choices.no_files
        
    %Read the file
    if iscell(caimanhandles.caimandr_choices.PathName)==0
        load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.FileName{filNum}])
    else
        load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.FileName{filNum}])
    end
     
     
    first_num_odor_trials(filNum)=num_odor_trials+1;
    per_file_trials_included.file(filNum).no_tr_inc=0;
    
    
    for trNo=1:no_odor_trials
        
        %Save epoch
        num_odor_trials=num_odor_trials+1;
        files_per_trial(num_odor_trials)=filNum;
        
        %Save lda
%         all_lda_events{num_odor_trials}=lda_event{trNo};
%         szlit=size(lda_input_timecourse);
%         all_lda_input_timecourse(1:length(time_to_event),1:szlit(2),num_odor_trials)=lda_input_timecourse(:,:,trNo);
%         all_lda_no_comp=szlit(2);
        
        if epoch_per_trial(trNo)==6
            %Hit
            epochs_per_trial(1,num_odor_trials)=1;
            epochs_per_trial(2:4,num_odor_trials)=0;
             
            %Was dF/F calculated?
            if sum(which_trial_Hit==trNo)>0
                %Calculate norm dFF
                
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Hit_traces(which_trial_Hit==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                per_file_trials_included.file(filNum).no_tr_inc=per_file_trials_included.file(filNum).no_tr_inc+1;
                per_file_trials_included.file(filNum).trials(per_file_trials_included.file(filNum).no_tr_inc)=trNo;
                per_file_trials_included.file(filNum).sessiontR(trNo)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(Hit_traces(which_trial_Hit==trNo,win),2);
                    szhit=size(Hit_traces(which_trial_Hit==trNo,win));
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    this_Hitii_lick=which_Hitii_lick(find(which_trial_Hit==trNo,1));
                    these_Hitii_lick_times=[];
                    these_Hitii_lick_times=Hit_lick_times(this_Hitii_lick,1:Hit_no_lick_times(this_Hitii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_Hitii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_Hitii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=1;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
                %Calculate the average snip for this trial
                Hitii=handles_out.Hit_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(Hitii).hit_traces);
                num_traces=handles_out.trialNo(Hitii).trace_numHit;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(Hitii).trace_numHit
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(Hitii).hit_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventHit;
                time_this_file.file(filNum).trial(trNo).time_to_event=handles_out.time_to_eventHit;
            end
        end
         
        if epoch_per_trial(trNo)==7
            %Miss
            epochs_per_trial(2,num_odor_trials)=1;
            epochs_per_trial(1,num_odor_trials)=0;
            epochs_per_trial(3:4,num_odor_trials)=0;
            
                        %Was dF/F calculated?
            if sum(which_trial_Miss==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Miss_traces(which_trial_Miss==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                per_file_trials_included.file(filNum).no_tr_inc=per_file_trials_included.file(filNum).no_tr_inc+1;
                per_file_trials_included.file(filNum).trials(per_file_trials_included.file(filNum).no_tr_inc)=trNo;
                per_file_trials_included.file(filNum).sessiontR(trNo)=num_odor_trials;
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(Miss_traces(which_trial_Miss==trNo,win),2);
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                     
                    %Calculate lick frequency for this window
                    this_Missii_lick=which_Missii_lick(find(which_trial_Miss==trNo,1));
                    these_Missii_lick_times=[];
                    these_Missii_lick_times=Miss_lick_times(this_Missii_lick,1:Miss_no_lick_times(this_Missii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_Missii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_Missii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=2;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                 
                %Calculate the average snip for this trial
                Missii=handles_out.Miss_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(Missii).miss_traces);
                num_traces=handles_out.trialNo(Missii).trace_numMiss;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(Missii).trace_numMiss
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(Missii).miss_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventMiss;
                time_this_file.file(filNum).trial(trNo).time_to_event=handles_out.time_to_eventMiss;
            end
        end
        
        if epoch_per_trial(trNo)==8
            %FA
            epochs_per_trial(3,num_odor_trials)=1;
            epochs_per_trial(1:2,num_odor_trials)=0;
            epochs_per_trial(4,num_odor_trials)=0;
            
            %Was dF/F calculated?
            if sum(which_trial_FA==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(FA_traces(which_trial_FA==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                per_file_trials_included.file(filNum).no_tr_inc=per_file_trials_included.file(filNum).no_tr_inc+1;
                per_file_trials_included.file(filNum).trials(per_file_trials_included.file(filNum).no_tr_inc)=trNo;
                per_file_trials_included.file(filNum).sessiontR(trNo)=num_odor_trials;
                for winNo=1:szwins(1)
                    %Calculate dFFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(FA_traces(which_trial_FA==trNo,win),2);
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    this_FAii_lick=which_FAii_lick(find(which_trial_FA==trNo,1));
                    these_FAii_lick_times=[];
                    these_FAii_lick_times=FA_lick_times(this_FAii_lick,1:FA_no_lick_times(this_FAii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_FAii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_FAii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=3;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
                %Calculate the average snip for this trial
                FAii=handles_out.FA_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(FAii).FA_traces);
                num_traces=handles_out.trialNo(FAii).trace_numFA;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(FAii).trace_numFA
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(FAii).FA_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventFA;
                time_this_file.file(filNum).trial(trNo).time_to_event=handles_out.time_to_eventFA;
            end
        end
        
        if epoch_per_trial(trNo)==9
            %CR
            epochs_per_trial(4,num_odor_trials)=1;
            epochs_per_trial(1:3,num_odor_trials)=0;
            
            %Was dF/F calculated?
            if sum(which_trial_CR==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(CR_traces(which_trial_CR==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                per_file_trials_included.file(filNum).no_tr_inc=per_file_trials_included.file(filNum).no_tr_inc+1;
                per_file_trials_included.file(filNum).trials(per_file_trials_included.file(filNum).no_tr_inc)=trNo;
                per_file_trials_included.file(filNum).sessiontR(trNo)=num_odor_trials;
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(CR_traces(which_trial_CR==trNo,win),2);
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    this_CRii_lick=which_CRii_lick(find(which_trial_CR==trNo,1));
                    these_CRii_lick_times=[];
                    these_CRii_lick_times=CR_lick_times(this_CRii_lick,1:CR_no_lick_times(this_CRii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_CRii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_CRii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                 
                 %Calculate the average snip for this trial
                CRii=handles_out.CR_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(CRii).CR_traces);
                num_traces=handles_out.trialNo(CRii).trace_numCR;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(CRii).trace_numCR
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(CRii).CR_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventCR;
                time_this_file.file(filNum).trial(trNo).time_to_event=handles_out.time_to_eventCR;
                
                epochs_per_trial_dFF(num_odor_trials_dFF)=4;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
            end
        end
  
    end
    noROIs(filNum)=szhit(1);
    
end

%Calculate percent correct
sliding_window=20; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;

perCorr=[];

%Note: I am moving the window for calculation of perCorr to the right by nine points
for ii=1:num_odor_trials-sliding_window+1  
    no_Hits=sum(epochs_per_trial(1,ii:ii+sliding_window-1)==1);
    no_CRs=sum(epochs_per_trial(4,ii:ii+sliding_window-1)==1);
    perCorr(ii+sliding_window-1)=100*(no_Hits+no_CRs)/sliding_window;
end


perCorr(1:sliding_window)=perCorr(2*sliding_window+1);

%Note, this is here so that perCorr=0 is included in the 0-10 % bin.
perCorr(perCorr==0)=0.00001;




%Plot percent correct vs trial
figNo=1;

try
    close(figNo)
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])

jj_low=find(perCorr<max_percent_low_beh);
plot(jj_low,perCorr(jj_low),'ob')
hold on
jj_high=find(perCorr>min_precent_high_beh);
plot(jj_high,perCorr(jj_high),'or')

jj_mid=find((perCorr<=min_precent_high_beh)&(perCorr>=max_percent_low_beh));
plot(jj_mid,perCorr(jj_mid),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot([0 num_odor_trials],[50 50],'-k')
 
%Draw the boundaries of each file
for filNum=2:caimanhandles.caimandr_choices.no_files
%     plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
    if isfield(caimanhandles.caimandr_choices,'start_reversal')
        if caimanhandles.caimandr_choices.start_reversal==filNum
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
        end
    end
    if isfield(caimanhandles.caimandr_choices,'start_gogo')
        if caimanhandles.caimandr_choices.start_gogo==filNum
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,80,'Go-go','Color','k','FontSize',18)
        end
    end
    if isfield(caimanhandles.caimandr_choices,'start_session')
        if sum(caimanhandles.caimandr_choices.start_session==filNum)>0
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
%             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
        end
    end
end
 
title(['Percent correct vs. trial number ' choiceFileName(18:end-2)])
xlabel('Trial number')
ylabel('Percent correct')
ylim([0 100])

%Now calculate the per ROI dFFs
noWins=size(caimanhandles.caimandr_choices.wins,1);
dFFavg_per_ROI=[];
glm_dFF=[];
glm_ii=0;

ii_stats=0;
dFF_stats=[];
this_ii=zeros(1,12);


trialNo=0;

for winNo=1:noWins
    for spmNo=1:2
        for perCorrNo=1:3
            dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).noComps=0;
        end
    end
end


for filNum=1:caimanhandles.caimandr_choices.no_files
    
    %Read the file
    if iscell(caimanhandles.caimandr_choices.PathName)==0
        load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.FileName{filNum}])
    else
        load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.FileName{filNum}])
    end
    

    %for each window
    for winNo=1:noWins
        %for each component
        this_no_comps=length(handles_out.file(filNum).trial_this_file(no_odor_trials).componentNo);
        for compNo=1:this_no_comps
            these_dFFs=[];
            these_spms=[];
            these_pcorrs=[];
            no_dFFs=0;
            for trNo=1:handles_out.file(filNum).no_odor_trials
                if sum(per_file_trials_included.file(filNum).trials==trNo)>0
                    no_dFFs=no_dFFs+1;
                    this_time=time_this_file.file(filNum).trial(trNo).time_to_event;
                    these_dFFs(no_dFFs)=mean(handles_out.file(filNum).trial_this_file(trNo).componentNo(compNo).traces((this_time>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (this_time<=caimanhandles.caimandr_choices.wins(winNo,2))));    
                    if (epoch_per_trial(trNo)==6)||(epoch_per_trial(trNo)==9)
                        these_spms(no_dFFs)=2;
                    else
                        these_spms(no_dFFs)=1;
                    end
                    this_perCorr=perCorr(per_file_trials_included.file(filNum).sessiontR(trNo));
                    if this_perCorr<max_percent_low_beh
                        these_pcorrs(no_dFFs)=1;
                    else
                        
                        if this_perCorr>min_precent_high_beh
                            these_pcorrs(no_dFFs)=2;
                        else
                            these_pcorrs(no_dFFs)=-1;
                        end
                        
                    end
                    
                end  
            end
            
            for perCorrNo=1:2
                for spmNo=1:2
                    dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).compNo(compNo).data_exist=0;
                    if sum((these_pcorrs==perCorrNo)&(these_spms==spmNo))>0
                        dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).noComps=dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).noComps+1;
                        this_mean_dFF=mean(these_dFFs((these_pcorrs==perCorrNo)&(these_spms==spmNo)));
                        dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).noComps)=this_mean_dFF;
                        glm_ii=glm_ii+1;
                        glm_dFF.data(glm_ii)=this_mean_dFF;
                        glm_dFF.spm(glm_ii)=spmNo;
                        glm_dFF.perCorr(glm_ii)=perCorrNo;
                        glm_dFF.winNo(glm_ii)=winNo;
                        
                        ii_stats=(winNo-1)*4+2*(spmNo-1)+perCorrNo-1+1;
                        this_ii(ii_stats)=this_ii(ii_stats)+1;
                        dFF_stats(ii_stats).data(this_ii(ii_stats))=this_mean_dFF;
                            
                    end
                end
            end
            
        end
    end
    
end

dFF_stats(1).description='Pre-odor naive S- ';
dFF_stats(3).description='Pre-odor naive S+ ';
dFF_stats(2).description='Pre-odor proficient S- ';
dFF_stats(4).description='Pre-odor proficient S+ ';

dFF_stats(5).description='Odor naive S- ';
dFF_stats(7).description='Odor naive S+ ';
dFF_stats(6).description='Odor proficient S- ';
dFF_stats(8).description='Odor proficient S+ ';

dFF_stats(9).description='Reinforcement naive S- ';
dFF_stats(11).description='Reinforcement naive S+ ';
dFF_stats(10).description='Reinforcement proficient S- ';
dFF_stats(12).description='Reinforcement proficient S+ ';


%Violin plot the dfFs
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])
hold on
  
x_pos=0;
edges=[-0.3:0.05:3];
rand_offset=0.4;
for perCorrNo=1:2
    for spmNo=1:2
        
        if perCorrNo==1
            for winNo=1:3
                if spmNo==1
                    bar(x_pos,mean(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF),'EdgeColor',[0.7 0.7 1.0],'FaceColor',[0.7 0.7 1.0])
                else
                    bar(x_pos,mean(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF),'EdgeColor',[1 0.7 0.7],'FaceColor',[1 0.7 0.7])
                end
                
                [mean_out, CIout]=drgViolinPoint(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF,edges,x_pos,rand_offset,'k','k',0.2);
                x_pos=x_pos+1;
            end
        else
            for winNo=1:3
                if spmNo==1
                    bar(x_pos,mean(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF),'b')
                else
                    bar(x_pos,mean(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF),'r')
                end
                [mean_out, CIout]=drgViolinPoint(dFFavg_per_ROI.win(winNo).spm(spmNo).perCorr(perCorrNo).dFF,edges,x_pos,rand_offset,'k','k',0.2);
                x_pos=x_pos+1;
            end
        end
        x_pos=x_pos+1;
    end
    x_pos=x_pos+2;
end
 
ylim([-0.3 1])
xticks([0 1 2 4 5 6 10 11 12 14 15 16])
xticklabels({'Pre','Odor','Reinf','Pre','Odor','Reinf','Pre','Odor','Reinf','Pre','Odor','Reinf'})
text(3,0.7,'<65%','FontSize',20)
text(13,0.7,'>80%','FontSize',20)

%Perform the glm for mi
fprintf(1, ['\n\nglm for dFF\n'])
tbl = table(glm_dFF.data',glm_dFF.spm',glm_dFF.perCorr',glm_dFF.winNo',...
    'VariableNames',{'dFF','spm','perCorr','window'});
mdl = fitglm(tbl,'dFF~spm+perCorr+window+spm*perCorr*window'...
    ,'CategoricalVars',[2,3,4])

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for dFF\n'])
try
    [output_data] = drgMutiRanksumorTtest(dFF_stats);
    fprintf(1, '\n\n')
catch
end

%odor vs pre
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .5 .5])
hold on
 
%Splus odor vs pre <65%
plot(dFFavg_per_ROI.win(1).spm(2).perCorr(1).dFF,dFFavg_per_ROI.win(2).spm(2).perCorr(1).dFF,'o','Color',[1.0 0.7 0.7])
plot(dFFavg_per_ROI.win(1).spm(1).perCorr(1).dFF,dFFavg_per_ROI.win(2).spm(1).perCorr(1).dFF,'o','Color',[0.7 0.7 1])
plot([-0.2 1],[-0.2 1],'-k')
xlim([-0.2 1])
ylim([-0.2 2])
title('dFF odor vs pre-odor for percent correct <65%')
xlabel('dFF pre-odor')
ylabel('dFF odor')

%odor vs pre
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .5 .5])
hold on
%Splus odor vs pre >80%
plot(dFFavg_per_ROI.win(1).spm(2).perCorr(2).dFF,dFFavg_per_ROI.win(2).spm(2).perCorr(2).dFF,'or')
plot(dFFavg_per_ROI.win(1).spm(1).perCorr(2).dFF,dFFavg_per_ROI.win(1).spm(2).perCorr(2).dFF,'ob')
plot([-0.2 1],[-0.2 1],'-k')
xlim([-0.2 1])
ylim([-0.2 2])
title('dFF odor vs pre-odor for percent correct >80%')

pffft=1;

