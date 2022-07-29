function drgCaImAnBatchPerSessionEventsPerTrialDimensionality(choiceBatchPathName,choiceFileName)

%This code does dimensionality analysis per trial
% Needs a choices file such as drgCaImAnChoicesDiego20180910_mmPVG04_Cerebellum
% Needs the output files from drgCaImAn_batch_dropc.m
warning('off')

close all
clear all

% handles_outs.lick_slopes=[];
% handles_outs.dFF_slopes=[];
% handles_outs.no_lick_slopes=0;
% handles_outs.no_dFF_slopes=0;
% handles_outs.no_vels=0;
% handles_outs.velocities=[];
% handles_outs.no_accels=0;
% handles_outs.acceleration=[];
% handles_outs.epochs_per_trial=[];
% handles_outs.perCorr=[];

dim_out=[];


min_trials=20;

tic

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAnBatchPerSessionReversalPerTrial run for ' choiceFileName '\n\n']);



addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

caimanhandles=handles;

%Read the files and calculate the dF/F in each window
num_odor_trials=0;
epochs_per_trial=[];
num_odor_trials_dFF=0;

all_lda_events=[];
all_lda_input_timecourse=[];
no_timepoints=2000000;
all_lda_events_miss_FA=[];

lick_times=[];
no_licks=[];
dLickTraces=[];

%Variables for optical flow
skip_ii=19;
baseline_ii=200;

handles_out=[];
these_times=[];

for filNum=1:caimanhandles.caimandr_choices.no_files
    
    %Read the CalmAn_batch file
    if iscell(caimanhandles.caimandr_choices.PathName)==0
        load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.FileName{filNum}])
    else
        load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.FileName{filNum}])
    end
    
%     %Read the Optical flow file
%     if iscell(caimanhandles.caimandr_choices.PathName)==0
%         load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.OpFileName{filNum}])
%     else
%         load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.OpFileName{filNum}])
%     end
%     
%     %Align the optical flow to the start of imaging
%     baseline_end_ii=100*frame_rate;
%     mag_meanlaser99=prctile(mag_meanlaser(1:baseline_end_ii),99);
%     mag_meanlaser1=prctile(mag_meanlaser(1:baseline_end_ii),1);
%     mean_laser_start=mean(mag_meanlaser(skip_ii+1:baseline_ii));
%     if exist('first_image_ii')==0
%         first_image_ii=find(mag_meanlaser(skip_ii+1:end)>0.5*(mag_meanlaser99-mag_meanlaser1)+mean_laser_start,1)+skip_ii;
%     end
%     timeOF=([1:length(mag_meantail)]/(frame_rate))-(first_image_ii/(frame_rate));
    
    first_num_odor_trials(filNum)=num_odor_trials+1;
    
    
    for trNo=1:no_odor_trials
        
        %Save epoch
        num_odor_trials=num_odor_trials+1;
        
        %Save lda
        all_lda_events{num_odor_trials}=lda_event{trNo};
        szlit=size(lda_input_timecourse);
        all_lda_input_timecourse(1:length(time_to_eventLDA),1:szlit(2),num_odor_trials)=lda_input_timecourse(:,:,trNo);
        all_lda_no_comp(num_odor_trials)=szlit(2);
        all_lda_fileNo(num_odor_trials)=filNum;
        no_timepoints=min([no_timepoints length(time_to_eventLDA)]);
        
        if epoch_per_trial(trNo)==6
            %Hit
            epochs_per_trial(1,num_odor_trials)=1;
            epochs_per_trial(2:4,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=1;
            
            %Was dF/F calculated?
            if sum(which_trial_Hit==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Hit_traces(which_trial_Hit==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save the licks for this trial
                this_Hitii_lick=which_Hitii_lick(find(which_trial_Hit==trNo,1));
                these_Hitii_lick_times=[];
                these_Hitii_lick_times=Hit_lick_times(this_Hitii_lick,1:Hit_no_lick_times(this_Hitii_lick));
                if ~isempty(these_Hitii_lick_times)
                    lick_times(num_odor_trials,1:length(these_Hitii_lick_times))=these_Hitii_lick_times;
                    no_licks(num_odor_trials)=length(these_Hitii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
%                 dLickTraces(num_odor_trials,:)=dHit_lick_traces(this_Hitii_lick,:);
                
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
                these_times(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventHit;
                
%                 %Get the movement
%                 snip_maskOF=(timeOF>=epoch_these_times(trNo)-dt_before)...
%                     &(timeOF<=epoch_these_times(trNo)+dt_after);
%                 velocity(num_odor_trials_dFF,1:sum(snip_maskOF))=mag_meantail(snip_maskOF);
%                 time_OFsnip(num_odor_trials_dFF,1:sum(snip_maskOF))=([1:sum(snip_maskOF)]*(1/frame_rate))-dt_before;
%                 no_timepointsOF(num_odor_trials_dFF)=sum(snip_maskOF);
            end
        end
        
        if epoch_per_trial(trNo)==7
            %Miss
            epochs_per_trial(2,num_odor_trials)=1;
            epochs_per_trial(1,num_odor_trials)=0;
            epochs_per_trial(3:4,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=2;
            
            %Was dF/F calculated?
            if sum(which_trial_Miss==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Miss_traces(which_trial_Miss==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save lick times
                this_Missii_lick=which_Missii_lick(find(which_trial_Miss==trNo,1));
                these_Missii_lick_times=[];
                these_Missii_lick_times=Miss_lick_times(this_Missii_lick,1:Miss_no_lick_times(this_Missii_lick));
                if ~isempty(these_Missii_lick_times)
                    lick_times(num_odor_trials,1:length(these_Missii_lick_times))=these_Missii_lick_times;
                    no_licks(num_odor_trials)=length(these_Missii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
%                 dLickTraces(num_odor_trials,:)=dMiss_lick_traces(this_Missii_lick,:);
                
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
                these_times(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventMiss;
                
%                 %Get the movement
%                 snip_maskOF=(timeOF>=epoch_these_times(trNo)-dt_before)...
%                     &(timeOF<=epoch_these_times(trNo)+dt_after);
%                 velocity(num_odor_trials_dFF,1:sum(snip_maskOF))=mag_meantail(snip_maskOF);
%                 time_OFsnip(num_odor_trials_dFF,1:sum(snip_maskOF))=([1:sum(snip_maskOF)]*(1/frame_rate))-dt_before;
%                 no_timepointsOF(num_odor_trials_dFF)=sum(snip_maskOF);
            end
        end
        
        if epoch_per_trial(trNo)==8
            %FA
            epochs_per_trial(3,num_odor_trials)=1;
            epochs_per_trial(1:2,num_odor_trials)=0;
            epochs_per_trial(4,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=4;
            
            %Was dF/F calculated?
            if sum(which_trial_FA==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(FA_traces(which_trial_FA==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save lick times
                this_FAii_lick=which_FAii_lick(find(which_trial_FA==trNo,1));
                these_FAii_lick_times=[];
                these_FAii_lick_times=FA_lick_times(this_FAii_lick,1:FA_no_lick_times(this_FAii_lick));
                if ~isempty(these_FAii_lick_times)
                    lick_times(num_odor_trials,1:length(these_FAii_lick_times))=these_FAii_lick_times;
                    no_licks(num_odor_trials)=length(these_FAii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
%                 dLickTraces(num_odor_trials,:)=dFA_lick_traces(this_FAii_lick,:);
                
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
                these_times(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventFA;
%                 
%                 %Get the movement
%                 snip_maskOF=(timeOF>=epoch_these_times(trNo)-dt_before)...
%                     &(timeOF<=epoch_these_times(trNo)+dt_after);
%                 velocity(num_odor_trials_dFF,1:sum(snip_maskOF))=mag_meantail(snip_maskOF);
%                 time_OFsnip(num_odor_trials_dFF,1:sum(snip_maskOF))=([1:sum(snip_maskOF)]*(1/frame_rate))-dt_before;
%                 no_timepointsOF(num_odor_trials_dFF)=sum(snip_maskOF);
            end
        end
        
        if epoch_per_trial(trNo)==9
            %CR
            epochs_per_trial(4,num_odor_trials)=1;
            epochs_per_trial(1:3,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=3;
            
            %Was dF/F calculated?
            if sum(which_trial_CR==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(CR_traces(which_trial_CR==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save lick times
                this_CRii_lick=which_CRii_lick(find(which_trial_CR==trNo,1));
                these_CRii_lick_times=[];
                these_CRii_lick_times=CR_lick_times(this_CRii_lick,1:CR_no_lick_times(this_CRii_lick));
                if ~isempty(these_CRii_lick_times)
                    lick_times(num_odor_trials,1:length(these_CRii_lick_times))=these_CRii_lick_times;
                    no_licks(num_odor_trials)=length(these_CRii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
%                 dLickTraces(num_odor_trials,:)=dCR_lick_traces(this_CRii_lick,:);
                
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
                these_times(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventCR;
                
                epochs_per_trial_dFF(num_odor_trials_dFF)=4;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
%                 %Get the movement
%                 snip_maskOF=(timeOF>=epoch_these_times(trNo)-dt_before)...
%                     &(timeOF<=epoch_these_times(trNo)+dt_after);
%                 velocity(num_odor_trials_dFF,1:sum(snip_maskOF))=mag_meantail(snip_maskOF);
%                 time_OFsnip(num_odor_trials_dFF,1:sum(snip_maskOF))=([1:sum(snip_maskOF)]*(1/frame_rate))-dt_before;
%                 no_timepointsOF(num_odor_trials_dFF)=sum(snip_maskOF);
            end
        end
        
        if num_odor_trials_dFF==45
            pffft=1;
        end
    end
    noROIs(filNum)=szhit(1);
end


%Trim the time course for LDA
all_lda_input_timecourse=all_lda_input_timecourse(1:no_timepoints,:,:);
time_to_eventLDA=time_to_eventLDA(1,1:no_timepoints);

%Calculate percent correct
sliding_window=min_trials; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;

perCorr=[];

%Note: Because this is a reversal I am moving the window for calculation of perCorr to the right by nine points
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


if ~isfield(caimanhandles.caimandr_choices,'start_reversal')
    caimanhandles.caimandr_choices.start_reversal=200;
end

%For reversals plot violin plot of percent
if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
    %Keep track of the percent correct
    %Plot percent correct vs trial
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])
    hold on
    
    %Parameters for violin plot
    edges=0:5:100;
    rand_offset=0.8;
    
    %Trials before reversal
    x_val=1;
    handles_out2.pctPerWin(1).pct=perCorr(1:first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1);
    bar(x_val,mean(handles_out2.pctPerWin(1).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(1).pct_mean, handles_out2.pctPerWin(1).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(1).pct,edges,x_val,rand_offset,'k','k',2);
    
    %Trials after reversal
    x_val=2;
    this_end=min([first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+65 length(perCorr)]);
    handles_out2.pctPerWin(2).pct=perCorr(first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+15:this_end);
    bar(x_val,mean(handles_out2.pctPerWin(2).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(2).pct_mean, handles_out2.pctPerWin(2).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(2).pct,edges,x_val,rand_offset,'k','k',2);
    
    %Trials at end
    x_val=3;
    handles_out2.pctPerWin(3).pct=perCorr(end-50:end);
    bar(x_val,mean(handles_out2.pctPerWin(3).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(3).pct_mean, handles_out2.pctPerWin(3).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(3).pct,edges,x_val,rand_offset,'k','k',2);
    
    %Draw lines between points
    plot([1 2 3],[mean(handles_out2.pctPerWin(1).pct) mean(handles_out2.pctPerWin(2).pct) mean(handles_out2.pctPerWin(3).pct)],'-k')
    
    ylim([0 120])
    xlim([0.3 3.7])
    ylabel('Percent correct')
    xticks([1 2 3 5 6 7])
    xticklabels({'Before','After','End'})
    
    
    
    window_labels{1}='Before';
    window_labels{2}='After';
    window_labels{3}='End';
    
    fprintf(1, ['\n\nranksum p values for percent correct windows\n\n']);
    
    no_pvals=0;
    for ii=1:3
        for jj=ii+1:3
            no_pvals=no_pvals+1;
            if (adtest(handles_out2.pctPerWin(ii).pct)==1)||(adtest(handles_out2.pctPerWin(ii).pct)==1)
                p_vals_corr(no_pvals)=ranksum(handles_out2.pctPerWin(ii).pct,handles_out2.pctPerWin(jj).pct);
                fprintf(1, ['p values ranksum for ' window_labels{ii} ' vs. ' window_labels{jj} ' =%d\n'],p_vals_corr(no_pvals));
            else
                [h p_vals_corr(no_pvals)]=ttest2(handles_out2.pctPerWin(ii).pct,handles_out2.pctPerWin(jj).pct);
                fprintf(1, ['p values t test for ' window_labels{ii} ' vs. ' window_labels{jj} ' =%d\n'],p_vals_corr(no_pvals));
            end
            
        end
    end
    
    pFDRcorr=drsFDRpval(p_vals_corr);
    fprintf(1, ['pFDR for significant difference percent correct  = %d\n\n'],pFDRcorr);
end

%plot number of ROIs
figNo=figNo+1;
try
    close figNo
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
hold on
plot([1:caimanhandles.caimandr_choices.no_files],noROIs,'-ok')
ylim([0 1.2*max(noROIs)])
title('Number of ROIs per file')
xlabel('File no')
ylabel('No ROIs')

dim_out.noROIs=noROIs;

%Plot norm dFF for each window
if caimanhandles.caimandr_choices.start_reversal>0
    if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
        tr_reversal=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal);
    else
        tr_reversal=200000;
    end
else
    tr_reversal=0;
end


for winNo=1:szwins(1)
    figNo=figNo+1;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .05+0.3*(winNo-1) .5 .25])
    hold on
    
    for dFF_trNo=1:num_odor_trials_dFF
        
        %Plot odor 1 (S+ forward)
        subplot(2,1,1)
        hold on
        if dFF_trNo<tr_reversal
            %If Hit or Miss
            if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
        else
            %If CR or FA
            if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
            
        end
        
        
        %Plot odor 2 (S- forward)
        subplot(2,1,2)
        hold on
        if dFF_trNo>=tr_reversal
            %If Hit or Miss
            if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
        else
            %If CR or FA
            if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
            
        end
    end
    
    for spno=1:2
        subplot(2,1,spno)
        if isfield(caimanhandles.caimandr_choices,'start_reversal')
            filNum=caimanhandles.caimandr_choices.start_reversal;
            if (filNum>0)&(filNum<caimanhandles.caimandr_choices.no_files)
                plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
                text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Reversal','Color','k','FontSize',18)
            end
        end
        
        if isfield(caimanhandles.caimandr_choices,'start_gogo')
            filNum=caimanhandles.caimandr_choices.start_gogo;
            if (filNum>0)&(filNum<caimanhandles.caimandr_choices.no_files)
                plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
                text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Go-go','Color','k','FontSize',18)
            end
        end
        
        if isfield(caimanhandles.caimandr_choices,'start_session')
            if length(caimanhandles.caimandr_choices.start_session)>=2
                for sessionNo=2:length(caimanhandles.caimandr_choices.start_session)
                    filNum=caimanhandles.caimandr_choices.start_session(sessionNo);
                    plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k')
                    %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
                end
            end
        end
        
        plot([1 num_odor_trials_dFF],[0 0], '-k')
        
        if spno==1
            title(['Odor 1 (S+ forward)'])
        else
            title(['Odor 2 (S- forward)'])
        end
        xlabel('Trial number')
        ylabel('Normalized dF/F')
        ylim([prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))])
    end
    sgtitle(['Normalized dF/F for window No ' num2str(winNo) ' Hit(red) Miss(cyan) FA(magenta) CR(blue)'])
end

%Do a linear discriminant analysis
for winNo=1:szwins(1)
    handles_sig.win(winNo).ii_for_sig=0;
end

if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
    %This is a forward run
    total_trial_windows=3;
    supertitle_description{1}='dF/F LDA analysis for trials with percent correct <65';
    supertitle_description{2}='dF/F LDA analysis for trials with percent correct >=65&<80';
    supertitle_description{3}='dF/F LDA analysis for trials with percent correct >=80';
else
    %Forward and reverse
    total_trial_windows=2;
    supertitle_description{1}='dF/F LDA analysis for trials before reversal';
    supertitle_description{2}='dF/F LDA analysis for trials after reversal at end of the session';
end

firstFig=figNo+1;
maxPC1=-2000;
maxPC2=-2000;
minPC1=20000;
minPC2=20000;

trial_window_description{1}='percent correct <65';
trial_window_description{2}='percent correct >=65&<80';
trial_window_description{3}='percent correct >=80';

lick_threshold=20; %This is the threshold to exclude the runs where Ming was adding water manually
t_odor_on=0;
t_odor_off=4;

for no_trial_windows=1:total_trial_windows
    handles_par(no_trial_windows).time_to_eventLDA=time_to_eventLDA;
    dFF_trial_mask=[];
    lick_excluded_trials=[];
    jj=0;
    events_miss_FA=[];
    which_trials_in_PCA=[];
    
    if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
        
        fprintf(1, '\n\nDimensionality processed for dF/F for trials before reversal \n');
        pct_windows=[45 65;65 80;80 100.1];
        
        for ii=1:num_odor_trials_dFF
            if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                lick_excluded_trials(ii)=0;
                if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
                    dFF_trial_mask(ii)=1;
                    jj=jj+1;
                    handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                    events{jj,1}=all_lda_events{ii};
                    events_miss_FA(jj)=all_lda_events_miss_FA(ii);
                    which_trials_in_PCA(jj)=ii;
                else
                    dFF_trial_mask(ii)=0;
                end
            else
                if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
                    lick_excluded_trials(ii)=1;
                else
                    lick_excluded_trials(ii)=0;
                end
                dFF_trial_mask(ii)=0;
            end
        end
    else
        if no_trial_windows==1
            %Forward trials
            fprintf(1, '\n\nLDA processed for dF/F for trials before reversal \n');
            
            for ii=1:num_odor_trials_dFF
                if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                    lick_excluded_trials(ii)=0;
                    if (trial_dFF(ii)>=1)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1)
                        dFF_trial_mask(ii)=1;
                        jj=jj+1;
                        handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                        events{jj,1}=all_lda_events{ii};
                    else
                        dFF_trial_mask(ii)=0;
                    end
                else
                    dFF_trial_mask(ii)=0
                    if (trial_dFF(ii)>=1)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1)
                        lick_excluded_trials(ii)=1;
                    else
                        lick_excluded_trials(ii)=0;
                    end
                end
            end
            
        else
            %Trials at end of reversal
            fprintf(1, '\n\nLDA processed for dF/F for trials after reversal \n');
            for ii=1:num_odor_trials_dFF
                if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                    lick_excluded_trials(ii)=0;
                    if (trial_dFF(ii)>=max(trial_dFF)-100)
                        dFF_trial_mask(ii)=1;
                        jj=jj+1;
                        handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                        events{jj,1}=all_lda_events{ii};
                    else
                        dFF_trial_mask(ii)=0;
                    end
                else
                    if (trial_dFF(ii)>=max(trial_dFF)-100)
                        lick_excluded_trials(ii)=1;
                    else
                        lick_excluded_trials(ii)=0;
                    end
                    dFF_trial_mask(ii)=0;
                end
            end
        end
        
    end
    
    
    
    if sum(dFF_trial_mask)>0
        %Do lick analysis
        %dt_lick=0.3;
        
        %Calculate dimensionality in the different windows
        for winNo=1:szwins(1)
            %Calculate dFF
            win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
            szalit=size(all_lda_input_timecourse);
            no_comps=szalit(2);
            
            %Get the data
            %Columns: cells, Rows: dF/F
            Signal=zeros(sum(win)*sum(dFF_trial_mask),no_comps);
            no_datapoints=0;
            
            for trNo=1:length(dFF_trial_mask)
                
                if dFF_trial_mask(trNo)
                    Signal(no_datapoints+1:no_datapoints+sum(win),:)=all_lda_input_timecourse(win,:,trNo);
                    no_datapoints=no_datapoints+sum(win);
                end
            end
            Dim_out(winNo) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
            dim_out.perCorr(no_trial_windows).window(winNo).dimensionality=Dim_out(winNo);
            
            switch no_trial_windows
                case 1
                    fprintf(1, ['Dimensionality for window %d and percent correct <65  = %d\n'], winNo, Dim_out(winNo));
                case 2
                    fprintf(1, ['Dimensionality for window %d and percent correct >=65&<80  = %d\n'], winNo, Dim_out(winNo));
                case 3
                    fprintf(1, ['Dimensionality for window %d and percent correct >80  = %d\n'], winNo, Dim_out(winNo));
                    
            end
        end
        
        
        
        dt_lick=time_to_eventLDA(2)-time_to_eventLDA(1);
        time_licks=([1:ceil((dt_after+dt_before)/dt_lick)]*dt_lick)-dt_before;
        delta_t_gauss=2; %seconds
        no_conv_points_lick=ceil(delta_t_gauss/(time_licks(2)-time_licks(1)));
        no_conv_points_dFF=ceil(delta_t_gauss/(time_to_eventLDA(2)-time_to_eventLDA(1)));
        
        
        
        %First figure out the lick threshold to exclude trials wheren Ming
        %gave the animal water manually during the odor on
        all_licks_per_dt_per_trial=zeros(num_odor_trials,ceil((dt_after+dt_before)/dt_lick));
        for trial_no=1:num_odor_trials
            for ii_lick=1:no_licks(trial_no)
                all_licks_per_dt_per_trial(trial_no, ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=all_licks_per_dt_per_trial(trial_no, ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
            end
        end
        
        
        hit_odor_lick_freq=zeros(1,sum(all_lda_events_miss_FA==1));
        hit_odor_lick_freq(1,:)=sum(all_licks_per_dt_per_trial(all_lda_events_miss_FA==1,(time_licks>=t_odor_on)&(time_licks<=t_odor_off)),2)/(t_odor_off-t_odor_on);
        
        %         lick_threshold=mean(hit_odor_lick_freq)+2*std(hit_odor_lick_freq);
        %         lick_threshold=20;
        %          lick_threshold=200;
        
        %Plot the lick frequency for S+ and S-
        Splick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        Smlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        
        sp_trno=0;
        sm_trno=0;
        
        for trial_no=1:num_odor_trials
            if dFF_trial_mask(trial_no)==1
                %if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))<lick_threshold
                if strcmp(all_lda_events{trial_no},'S+')
                    %S+
                    for ii_lick=1:no_licks(trial_no)
                        Splick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Splick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                    end
                    sp_trno=sp_trno+1;
                else
                    %S-
                    for ii_lick=1:no_licks(trial_no)
                        Smlick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Smlick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                    end
                    sm_trno=sm_trno+1;
                end
                %end
            end
        end
        
        Splick_freq=(Splick_freq/(sp_trno*dt_lick));
        Smlick_freq=(Smlick_freq/(sm_trno*dt_lick));
        
        %Convolve lick_freq using a window of 0.9 sec
        %         no_conv_points=3;
        %         conv_win=ones(1,no_conv_points);
        
        %Convolve lick_freq using a Gaussian window
        conv_win=gausswin(no_conv_points_lick);
        
        Splick_freq=conv(Splick_freq,conv_win,'same')/sum(conv_win);
        Smlick_freq=conv(Smlick_freq,conv_win,'same')/sum(conv_win);
        
        
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        
        hold on
        
        lfreqmax=20;
        
        p1=plot(time_licks,Smlick_freq,'-b','LineWidth',2);
        p2=plot(time_licks,Splick_freq,'-r','LineWidth',2);
        
        %Odor on markers
        plot([0 0],[-1 lfreqmax],'-k')
        odorhl=plot([0 mean(delta_odor)],[-0.5 -0.5],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[-1 lfreqmax],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[-1 lfreqmax],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-0.5 -0.5],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-1 lfreqmax],'-r')
        
        switch no_trial_windows
            case 1
                title(['Lick frequency for percent correct <65']);
            case 2
                title(['Lick frequency for percent correct >=65&<80']);
            case 3
                title(['Lick frequency for percent correct >80']);
        end
        
        xlabel('Time (sec)')
        ylabel('Lick frequency (Hz)')
        legend([p1 p2],'S+','S-')
        xlim([-10 20])
        ylim([-1 lfreqmax])
        
        %Now plot the p value for the difference in licks between S+ and S-
        %Get lick p values
        try
            
            no_pvals=0;
            p_val_Sp_Sm=[];
            time_p_lick=[];
            
            for ii=1:length(time_licks)-1
                
                sp_trno=0;
                sm_trno=0;
                this_Sp=[];
                this_Sm=[];
                
                for trial_no=1:num_odor_trials
                    %                     if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))<lick_threshold
                    if dFF_trial_mask(trial_no)==1
                        if strcmp(all_lda_events{trial_no},'S+')
                            %S+
                            sp_trno=sp_trno+1;
                            this_Sp(sp_trno,1)=0;
                            for ii_lick=1:no_licks(trial_no)
                                if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                    this_Sp(sp_trno,1)=1;
                                end
                            end
                        else
                            %S-
                            sm_trno=sm_trno+1;
                            this_Sm(sm_trno,1)=0;
                            for ii_lick=1:no_licks(trial_no)
                                if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                    this_Sm(sm_trno,1)=1;
                                end
                            end
                        end
                    end
                    %                     end
                end
                
                
                no_pvals=no_pvals+1;
                
                if (~isempty(this_Sm))&(~isempty(this_Sp))
                    p_val_Sp_Sm(no_pvals)=ranksum(this_Sm,this_Sp);
                else
                    p_val_Sp_Sm(no_pvals)=1;
                end
                
                %ranksum gives NaN if the values are all the same
                if isnan(p_val_Sp_Sm(no_pvals))
                    p_val_Sp_Sm(no_pvals)=1;
                end
                
                time_p_lick(no_pvals)= (time_licks(ii)+time_licks(ii+1))/2;
                
            end
            
            logpmin=-20;
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            figure(figNo)
            plot(time_p_lick,log10(p_val_Sp_Sm),'-k','LineWidth',2)
            hold on
            plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)],'-r','LineWidth',1)
            
            %Odor on markers
            plot([0 0],[logpmin 0.5],'-k')
            odorhl=plot([0 mean(delta_odor)],[logpmin+0.5 logpmin+0.5],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],[logpmin 0.5],'-k')
            
            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[logpmin 0.5],'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin+0.5 logpmin+0.5],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin 0.5],'-r')
            
            
            switch no_trial_windows
                case 1
                    title(['log(p value) for the difference in licks Sp vs Sm for percent correct <65']);
                case 2
                    title(['log(p value) for the difference in licks Sp vs Sm for percent correct >=65&<80']);
                case 3
                    title(['log(p value) for the difference in licks Sp vs Sm for percent correct >80']);
            end
            
            xlabel('Time (sec)')
            ylabel('log10(p value)')
            ylim([logpmin 0.5])
            xlim([-10 20])
            
        catch
        end
%         
%         %Plot the licks
%         figNo=figNo+1;
%         try
%             close(figNo)
%         catch
%         end
%         
%         figure(figNo)
%         hold on
%         
%         
%         per99=prctile(dLickTraces(:),99.9);
%         per1=prctile(dLickTraces(:),1);
%         szalllick=size(dLickTraces);
%         time_licksd=([1:szalllick(2)]/(acq_rate/20))-dt_before;
%         y_shift=0;
%         
%         %Plot Sp lick traces
%         for trial_no=1:num_odor_trials
%             if strcmp(all_lda_events{trial_no},'S+')
%                 if dFF_trial_mask(trial_no)==1
%                     plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-r')
%                     y_shift=y_shift+1.2*(per99-per1);
%                 else
%                     if lick_excluded_trials(trial_no)==1
%                         plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-k')
%                         y_shift=y_shift+1.2*(per99-per1);
%                     end
%                 end
%             end
%             
%             
%             
%         end
%         
%         %Plot Sm lick traces
%         for trial_no=1:num_odor_trials
%             if strcmp(all_lda_events{trial_no},'S-')
%                 if dFF_trial_mask(trial_no)==1
%                     plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-b')
%                     y_shift=y_shift+1.2*(per99-per1);
%                 else
%                     if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))>lick_threshold
%                         plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-k')
%                         y_shift=y_shift+1.2*(per99-per1);
%                     end
%                 end
%             end
%         end
%         
%         %Odor on markers
%         plot([0 0],[0 y_shift],'-k')
%         plot([mean(delta_odor) mean(delta_odor)],[0 y_shift],'-k')
%         
%         %Reinforcement markers
%         plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 y_shift],'-r')
%         plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 y_shift],'-r')
%         
%         switch no_trial_windows
%             case 1
%                 title(['lick traces per trial for percent correct <65']);
%             case 2
%                 title(['lick traces per trial for percent correct >=65&<80']);
%             case 3
%                 title(['lick traces per trial for percent correct >80']);
%         end
%         
%         xlabel('Time (sec)')
%         
%         xlim([-10 20])
%         
        
        
        
        
        %Plot each trial for dF/F, lick frequency, velocity and their
        %time derivatives
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .05 .9 .9])
        
        sgtitle(['Timecourse for dF/F, lick frequency (LF) and velocity and derivatives for ' supertitle_description{no_trial_windows}])
        
        %First dF/F
        maxdFF=-200;
        mindFF=200;
        
        maxdFF_dx=-200;
        mindFF_dx=200;
        
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                conv_win=gausswin(no_conv_points_dFF);
                
                this_conv_dFF=[];
                this_conv_dFF=conv(mean_snip_dFF(trNo,:),conv_win,'same')/sum(conv_win);
                
                this_pct95=prctile(this_conv_dFF(1:132),95);
                maxdFF=max([maxdFF this_pct95]);
                
                this_pct5=prctile(this_conv_dFF(1:132),5);
                mindFF=min([mindFF this_pct5]);
                
                
                %Calculate the derivative of dF/F
                this_conv_dFF_dx=gradient(this_conv_dFF);
                
                this_pct95=prctile(this_conv_dFF_dx(1:132),95);
                maxdFF_dx=max([maxdFF_dx this_pct95]);
                
                this_pct5=prctile(this_conv_dFF_dx(1:132),5);
                mindFF_dx=min([mindFF_dx this_pct5]);
                
            end
        end
        
        ymax=maxdFF+0.1*(maxdFF-mindFF);
        ymin=mindFF-0.1*(maxdFF-mindFF);
        
        ymax_dx=maxdFF_dx+0.1*(maxdFF_dx-mindFF_dx);
        ymin_dx=mindFF_dx-0.1*(maxdFF_dx-mindFF_dx);
        
        t_offset=0;
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                
                evNo=all_lda_events_miss_FA(trNo);
                
                %dF/F
                subplot(6,1,1)
                hold on
                
                conv_win=gausswin(no_conv_points_dFF);
                
                this_conv_dFF=[];
                this_conv_dFF=conv(mean_snip_dFF(trNo,:),conv_win,'same')/sum(conv_win);
                
                
                switch evNo
                    case 1
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','r','LineWidth',2);
                    case 2
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','c','LineWidth',2);
                    case 3
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','b','LineWidth',2);
                    case 4
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','m','LineWidth',2);
                end
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin ymax],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
                
                
                
                %d(dF/F)/dt
                subplot(6,1,2)
                hold on
                
                %Calculate the derivative of dF/F
                this_conv_dFF_dx=gradient(this_conv_dFF);
                
                switch evNo
                    case 1
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','r','LineWidth',2);
                    case 2
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','c','LineWidth',2);
                    case 3
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','b','LineWidth',2);
                    case 4
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','m','LineWidth',2);
                end
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin_dx ymax_dx],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin_dx ymax_dx],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin_dx ymax_dx],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx ymax_dx],'-r')
                
                
%                 handles_outs.no_dFF_slopes=handles_outs.no_dFF_slopes+1;
%                 handles_outs.dFF(handles_outs.no_dFF_slopes,1:132)=this_conv_dFF(1:132);
%                 handles_outs.dFF_derivatives(handles_outs.no_dFF_slopes,1:132)=this_conv_dFF_dx(1:132);
%                 handles_outs.time_to_eventLDA=time_to_eventLDA(1:132);
%                 handles_outs.epochs_per_trial(handles_outs.no_dFF_slopes)=epochs_per_trial_dFF(trNo);
%                 handles_outs.perCorr(handles_outs.no_dFF_slopes)=perCorr(trNo);
%                 
                
                t_offset=t_offset+35;

            end
            
            
        end
        
        %         xlim([-5 t_offset-25])
        subplot(6,1,1)
        xlim([-50 350])
        ylim([ymin ymax])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('dF/F')
%         title(['Timecourse for dF/F for ' supertitle_description{no_trial_windows}])
        
        
        subplot(6,1,2)
        xlim([-50 350])
        ylim([ymin_dx ymax_dx])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('d(dF/F)/dt')
%         title(['Timecourse for derivative dF/F for ' supertitle_description{no_trial_windows}])
        
%         %Now plot the velocity
%         %First get ymin, ymax, etc
%         maxvel=-200;
%         minvel=200;
%         maxacc=-200;
%         minacc=200;
%         for trNo=1:length(dFF_trial_mask)
%             
%             if dFF_trial_mask(trNo)
%                 
%                 
%                 velocity(num_odor_trials_dFF,1:sum(snip_maskOF))=mag_meantail(snip_maskOF);
%                 time_OFsnip(num_odor_trials_dFF,1:sum(snip_maskOF))=([1:sum(snip_maskOF)]*(1/frame_rate))-dt_before;
%                 no_timepointsOF(num_odor_trials_dFF)=sum(snip_maskOF);
%                 
%                 no_conv_points_OF=ceil(delta_t_gauss/(time_OFsnip(trNo,2)-time_OFsnip(trNo,1)));
%                 
%                 conv_win=gausswin(no_conv_points_OF);
%                 
%                 this_conv_velocity=[];
%                 this_conv_velocity=conv(velocity(trNo,1:no_timepointsOF(trNo)),conv_win,'same')/sum(conv_win);
%                 
%                 this_time_OFsnip=[];
%                 this_time_OFsnip=time_OFsnip(trNo,1:no_timepointsOF(trNo));
%                 
%                 this_conv_velocity_decimated=[];
%                 for time_ii=1:length(time_to_eventLDA)
%                     point_before=find(this_time_OFsnip<time_to_eventLDA(time_ii),1,'last');
%                     point_after=find(this_time_OFsnip>time_to_eventLDA(time_ii),1,'first');
%                     if ~isempty(point_after)
%                         this_conv_velocity_decimated(1,time_ii)=(this_conv_velocity(point_after)+this_conv_velocity(point_before))/2;
%                     else
%                         this_conv_velocity_decimated(1,time_ii)=this_conv_velocity(point_before);
%                     end
%                 end
%                 
%                 
%                 this_pct95=prctile(this_conv_velocity_decimated(1:132),95);
%                 maxvel=max([maxvel this_pct95]);
%                 
%                 this_pct5=prctile(this_conv_velocity_decimated(1:132),5);
%                 minvel=min([minvel this_pct5]);
%                 
%                 %Acceleration
%                 %Calculate the derivative of the velocity
%                 this_acceleration=gradient(this_conv_velocity_decimated);
%                 
%                 this_pct95=prctile(this_acceleration,95);
%                 maxacc=max([maxacc this_pct95]);
%                 
%                 this_pct5=prctile(this_acceleration,5);
%                 minacc=min([minacc this_pct5]);
%                 
%             end
%             
%             
%         end
%         
%         ymax=maxvel+0.1*(maxvel-minvel);
%         ymin=minvel-0.1*(maxvel-minvel);
%         
%         t_offset=0;
%         
%         ymax_accel=maxacc+0.1*(maxacc-minacc);
%         ymin_accel=minacc-0.1*(maxacc-minacc);
%         
%         t_offset_accel=0;
%         
%         velocity_decimated=zeros(length(dFF_trial_mask),length(time_to_eventLDA));
%         for trNo=1:length(dFF_trial_mask)
%             
%             if dFF_trial_mask(trNo)
%                 
%                 evNo=all_lda_events_miss_FA(trNo);
%                 
%                 %Velocity
%                 subplot(6,1,5)
%                 hold on
%                 velocity(num_odor_trials_dFF,1:sum(snip_maskOF))=mag_meantail(snip_maskOF);
%                 time_OFsnip(num_odor_trials_dFF,1:sum(snip_maskOF))=([1:sum(snip_maskOF)]*(1/frame_rate))-dt_before;
%                 no_timepointsOF(num_odor_trials_dFF)=sum(snip_maskOF);
%                 
%                 no_conv_points_OF=ceil(delta_t_gauss/(time_OFsnip(trNo,2)-time_OFsnip(trNo,1)));
%                 conv_win=gausswin(no_conv_points_OF);
%                 
%                 this_conv_velocity=[];
%                 this_conv_velocity=conv(velocity(trNo,1:no_timepointsOF(trNo)),conv_win,'same')/sum(conv_win);
%                 
%                 this_time_OFsnip=[];
%                 this_time_OFsnip=time_OFsnip(trNo,1:no_timepointsOF(trNo));
%                 
%                 for time_ii=1:length(time_to_eventLDA)
%                     point_before=find(this_time_OFsnip<time_to_eventLDA(time_ii),1,'last');
%                     point_after=find(this_time_OFsnip>time_to_eventLDA(time_ii),1,'first');
%                     if ~isempty(point_after)
%                         velocity_decimated(trNo,time_ii)=(this_conv_velocity(point_after)+this_conv_velocity(point_before))/2;
%                     else
%                         velocity_decimated(trNo,time_ii)=this_conv_velocity(point_before);
%                     end
%                 end
%                 
%                 this_conv_velocity_decimated=[];
%                 this_conv_velocity_decimated=velocity_decimated(trNo,:);
%                 
%                 switch evNo
%                     case 1
%                         plot(time_to_eventLDA(1:132)'+t_offset,this_conv_velocity_decimated(1:132)','r','LineWidth',2);
%                     case 2
%                         plot(time_to_eventLDA(1:132)'+t_offset,this_conv_velocity_decimated(1:132)','c','LineWidth',2);
%                     case 3
%                         plot(time_to_eventLDA(1:132)'+t_offset,this_conv_velocity_decimated(1:132)','b','LineWidth',2);
%                     case 4
%                         plot(time_to_eventLDA(1:132)'+t_offset,this_conv_velocity_decimated(1:132)','m','LineWidth',2);
%                 end
%                 
%                 
%                 %Odor on markers
%                 plot([0+t_offset 0+t_offset],[ymin ymax],'-k')
%                 odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
%                 plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
%                 
%                 %Reinforcement markers
%                 plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
%                 reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
%                 plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
%                 
%                 
%                
%                 
%                 t_offset=t_offset+35;
%                 
%                 %Acceleration
%                 subplot(6,1,6)
%                 hold on
%                
%                 %Acceleration
%                 %Calculate the derivative of the velocity
%                 this_acceleration=gradient(this_conv_velocity_decimated);
%                 
%                 switch evNo
%                     case 1
%                         plot(time_to_eventLDA(1:132)'+t_offset_accel,this_acceleration(1:132)','r','LineWidth',2);
%                     case 2
%                         plot(time_to_eventLDA(1:132)'+t_offset_accel,this_acceleration(1:132)','c','LineWidth',2);
%                     case 3
%                         plot(time_to_eventLDA(1:132)'+t_offset_accel,this_acceleration(1:132)','b','LineWidth',2);
%                     case 4
%                         plot(time_to_eventLDA(1:132)'+t_offset_accel,this_acceleration(1:132)','m','LineWidth',2);
%                 end
%                 
%                 
%                 %Odor on markers
%                 plot([0+t_offset 0+t_offset],[ymin_accel ymax_accel],'-k')
%                 odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin_accel + 0.1*(ymax_accel-ymin_accel) ymin_accel + 0.1*(ymax_accel-ymin_accel)],'-k','LineWidth',5);
%                 plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin_accel ymax_accel],'-k')
%                 
%                 %Reinforcement markers
%                 plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin_accel ymax_accel],'-r')
%                 reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_accel + 0.1*(ymax_accel-ymin_accel) ymin_accel + 0.1*(ymax_accel-ymin_accel)],'-r','LineWidth',5);
%                 plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_accel ymax_accel],'-r')
%                 
%                 handles_outs.no_vels=handles_outs.no_vels+1;
%                 handles_outs.velocities(handles_outs.no_vels,1:132)=this_conv_velocity_decimated(1:132);
%                 handles_outs.no_accels=handles_outs.no_accels+1;
%                 handles_outs.acceleration(handles_outs.no_vels,1:132)=this_acceleration(1:132);
%                 
%                 t_offset_accel=t_offset_accel+35;
%                 
%             end
%             
%             
%         end
%         
%         
%         subplot(6,1,5)
%         xlim([-50 350])
%         ylim([ymin ymax])
%         %         ylim([-0.4 1])
%         xlabel('Time (sec)')
%         ylabel('dx/dt')
% %         title(['Timecourse for velocity for ' supertitle_description{no_trial_windows}])
%         
%         subplot(6,1,6)
%         xlim([-50 350])
%         ylim([ymin_accel ymax_accel])
%         %         ylim([-0.4 1])
%         xlabel('Time (sec)')
%         ylabel('dv/dt')
% %         title(['Timecourse for acceleration for ' supertitle_description{no_trial_windows}])
%         
%         
        
        %lick frequency
        
        
        maxlick=-200;
        minlick=200;
        
        maxlick_dx=-200;
        minlick_dx=200;
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                
                this_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                for ii_lick=1:no_licks(trNo)
                    this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))=this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))+1;
                end
                
                %Convolve lick_freq using a window of 0.9 sec
                %                 no_conv_points=4;
                %                 conv_win=ones(1,no_conv_points);
                
                %Convolve lick_freq using a Gaussian window
                conv_win=gausswin(no_conv_points_lick);
                
                lick_freq=conv(this_lick_freq,conv_win,'same')/sum(conv_win);
                lick_freq=lick_freq/dt_lick;
                
                this_pct95=prctile(lick_freq,95);
                maxlick=max([maxlick this_pct95]);
                
                this_pct5=prctile(lick_freq,5);
                minlick=min([minlick this_pct5]);
                
                %Calculate the derivative of lick_freq
                lick_freq_dx=gradient(lick_freq);
                
                this_pct95=prctile(lick_freq_dx,95);
                maxlick_dx=max([maxlick_dx this_pct95]);
                
                this_pct5=prctile(lick_freq_dx,5);
                minlick_dx=min([minlick_dx this_pct5]);
                
                
            end
        end
        
        ymax=maxlick+0.1*(maxlick-minlick);
        ymin=minlick-0.1*(maxlick-minlick);
        
        ymax_dx=maxlick_dx+0.1*(maxlick_dx-minlick_dx);
        ymin_dx=minlick_dx-0.1*(maxlick_dx-minlick_dx);
        
        t_offset=0;
        
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                
                evNo=all_lda_events_miss_FA(trNo);
                
                
                this_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                for ii_lick=1:no_licks(trNo)
                    this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))=this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))+1;
                end
                
                
                %Convolve lick_freq using a flat window
                %                 no_conv_points=4;
                %                 conv_win=ones(1,no_conv_points);
                
                %Convolve lick_freq using a Gaussian window
                conv_win=gausswin(no_conv_points_lick);
                
                lick_freq=conv(this_lick_freq,conv_win,'same')/sum(conv_win);
                lick_freq=lick_freq/dt_lick;

                
                subplot(6,1,3)
                hold on
                
                switch evNo
                    case 1
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'r','LineWidth',2);
                    case 2
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'c','LineWidth',2);
                    case 3
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'b','LineWidth',2);
                    case 4
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'m','LineWidth',2);
                end
                
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin ymax],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
                
                %Calculate the derivative of lick_freq
                lick_freq_dx=gradient(lick_freq);
                
                %Plot the derivative
                subplot(6,1,4)
                hold on
                
                switch evNo
                    case 1
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'r','LineWidth',2);
                    case 2
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'c','LineWidth',2);
                    case 3
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'b','LineWidth',2);
                    case 4
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'m','LineWidth',2);
                end
                
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin_dx ymax_dx],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin_dx ymax_dx],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin_dx ymax_dx],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx ymax_dx],'-r')
                
               
%                 handles_outs.no_lick_slopes=handles_outs.no_lick_slopes+1;
%                 handles_outs.lick_freq(handles_outs.no_lick_slopes,1:132)=lick_freq(1:132);
%                 handles_outs.lick_derivatives(handles_outs.no_lick_slopes,1:132)=lick_freq_dx(1:132);
%                 handles_outs.time_licks=time_licks(1:132);
                
                t_offset=t_offset+35;
                
            end
        end
        
        %         xlim([-5 t_offset-25])
        subplot(6,1,3)
        xlim([-50 350])
        ylim([ymin ymax])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('LR')
%         title(['Timecourse for lick rate for ' supertitle_description{no_trial_windows}])
        
        subplot(6,1,4)
        xlim([-50 350])
        ylim([ymin_dx ymax_dx])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('dLR/dt')
%         title(['Timecourse for the drivative of lick rate for ' supertitle_description{no_trial_windows}])
        
        pffft=1;
        
        
        
        
        %Plot the timecourse for dF/F
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        
        hold on
        
        maxdFF=-200;
        mindFF=200;
        
        for evNo=[4 2 3 1]
            %Find out how many trials for this event
            
            these_dFF_trials=logical(dFF_trial_mask)&(all_lda_events_miss_FA==evNo);
            
            
            if sum(these_dFF_trials)>0
                
                
                if sum(these_dFF_trials)>2
                    CI=[];
                    CI = bootci(1000, {@mean, mean_snip_dFF(these_dFF_trials,:)})';
                    CI(:,1)=mean(mean_snip_dFF(these_dFF_trials,:))'-CI(:,1);
                    CI(:,2)=CI(:,2)-mean(mean_snip_dFF(these_dFF_trials,:))';
                    
                    maxdFF=max([maxdFF,max(CI(:,1)+mean(mean_snip_dFF(these_dFF_trials,:))')]);
                    mindFF=min([mindFF,min(mean(mean_snip_dFF(these_dFF_trials,:))'-CI(:,2))]);
                    
                    switch evNo
                        case 1
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'r');
                        case 2
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'c');
                        case 3
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'b');
                        case 4
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'm');
                    end
                else
                    switch evNo
                        case 1
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))','r');
                        case 2
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', 'c');
                        case 3
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', 'b');
                        case 4
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))',  'm');
                    end
                end
            end
            
        end
        
        ymax=maxdFF+0.1*(maxdFF-mindFF);
        ymin=mindFF-0.1*(maxdFF-mindFF);
        
        %Odor on markers
        plot([0 0],[ymin ymax],'-k')
        odorhl=plot([0 mean(delta_odor)],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[ymin ymax],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[ymin ymax],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[ymin ymax],'-r')
        
        xlim([-10 20])
        %         ylim([ymin ymax])
        ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('dF/F')
        title(['Timecourse for dF/F for ' supertitle_description{no_trial_windows}])
        
        
        
        pffft=1;
        
    end
    pffft=1;
end
% 
% %Plot the lick derivative plot
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% set(hFig, 'units','normalized','position',[.1 .25 .75 .25])
% 
% fprintf(1, ['\n\n\nProcessing rho and p value for dFF derivative vs lick derivative\n\n']);
% 
% %For the entire timecourse
% subplot(1,4,1)
% hold on
% 
% lick_derivatives=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,:)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,:)];
%     plot(handles_outs.lick_derivatives(trNo,:),handles_outs.dFF_derivatives(trNo,:),'.b')
% end
% 
% ymax=max(dFF_derivatives)+0.1*(max(dFF_derivatives)-min(dFF_derivatives));
% ymin=min(dFF_derivatives)-0.1*(max(dFF_derivatives)-min(dFF_derivatives));
% 
% xmax=max(lick_derivatives)+0.1*(max(lick_derivatives)-min(lick_derivatives));
% xmin=min(lick_derivatives)-0.1*(max(lick_derivatives)-min(lick_derivatives));
% 
% %Calculate the density for the pseudocolor
% no_bins=50;
% entire_dFFdervsld_density=zeros(no_bins,no_bins);
% entire_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     entire_no_points=entire_no_points+length(lick_derivatives);
%     for ii=1:length(lick_derivatives)
%         entire_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             entire_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(lick_derivatives',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs lick derivative  = %d, %d\n'],rho,pval);
% handles_outs.rho(1)=rho;
% handles_outs.pval(1)=pval;
% xlabel('lick frequency derivative (Hz/sec)')
% ylabel('dF/F derivative (1/sec)')
% title('All times')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Now do the three windows
% time_licks=[];
% time_licks=handles_outs.time_licks;
% 
% %Before odor window
% subplot(1,4,2)
% hold on
% 
% time_mask=(time_licks>=-mean(delta_odor))&(time_licks<=0);
% lick_derivatives=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,time_mask)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
%     plot(handles_outs.lick_derivatives(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
% end
% 
% %Calculate the density for the pseudocolor
% no_bins=50;
% pre_odor_dFFdervsld_density=zeros(no_bins,no_bins);
% pre_odor_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     pre_odor_no_points=pre_odor_no_points+length(lick_derivatives);
%     for ii=1:length(lick_derivatives)
%         pre_odor_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             pre_odor_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% 
% tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(lick_derivatives',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs lick derivative before odor = %d, %d\n'],rho,pval);
% handles_outs.rho(2)=rho;
% handles_outs.pval(2)=pval;
% xlabel('lick frequency derivative (Hz/sec)')
% ylabel('dF/F derivative (1/sec)')
% title('Before odor')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Odor window
% subplot(1,4,3)
% hold on
% 
% 
% time_mask=(time_licks>=0)&(time_licks<=mean(delta_odor));
% lick_derivatives=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,time_mask)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
%     plot(handles_outs.lick_derivatives(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
% end
% 
% %Calculate the density for the pseudocolor
% no_bins=50;
% odor_dFFdervsld_density=zeros(no_bins,no_bins);
% odor_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     odor_no_points=odor_no_points+length(lick_derivatives);
%     for ii=1:length(lick_derivatives)
%         odor_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             odor_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(lick_derivatives',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs lick derivative during odor = %d, %d\n'],rho,pval);
% handles_outs.rho(3)=rho;
% handles_outs.pval(3)=pval;
% xlabel('lick frequency derivative (Hz/sec)')
% ylabel('dF/F derivative (1/sec)')
% title('During odor')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Reinforcement window
% subplot(1,4,4)
% hold on
% 
% time_mask=(time_licks>=mean(delta_odor_on_reinf_on))&(time_licks<=mean(delta_odor_on_reinf_on)+3);
% lick_derivatives=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,time_mask)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
%     plot(handles_outs.lick_derivatives(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
% end
% 
% %Calculate the density for the pseudocolor
% no_bins=50;
% reinforcement_dFFdervsld_density=zeros(no_bins,no_bins);
% reinforcement_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     reinforcement_no_points=reinforcement_no_points+length(lick_derivatives);
%     for ii=1:length(lick_derivatives)
%         reinforcement_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             reinforcement_dFFdervsld_density(fix(no_bins*(lick_derivatives(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(lick_derivatives',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs lick derivative during reinforcement = %d, %d\n'],rho,pval);
% handles_outs.rho(4)=rho;
% handles_outs.pval(4)=pval;
% xlabel('lick frequency derivative (Hz/sec)')
% ylabel('dF/F derivative (1/sec)')
% title('Reinforcement')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Now do a pseudocolor plot
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% set(hFig, 'units','normalized','position',[.1 .25 .75 .25])
% 
% %Entire timecourse
% subplot(1,4,1)
% hold on
% entire_dFFdervsld_density=entire_dFFdervsld_density/entire_no_points;
% max_ent=max(entire_dFFdervsld_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(entire_dFFdervsld_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-4    log10(max_ent)])
% xlabel('Derivative of lick freq (Hz/sec)')
% ylabel('dFF derivative (1/sec)');
% title(['Entire timecourse'])
% xlim([min(lick_derivatives) max(lick_derivatives)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% %Pre-odor timecourse
% subplot(1,4,2)
% hold on
% pre_odor_dFFdervsld_density=pre_odor_dFFdervsld_density/pre_odor_no_points;
% max_ent=max(pre_odor_dFFdervsld_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(pre_odor_dFFdervsld_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-4    log10(max_ent)])
% xlabel('Derivative of lick freq (Hz/sec)')
% ylabel('dFF derivative (1/sec)');
% title(['Pre-odor'])
% xlim([min(lick_derivatives) max(lick_derivatives)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% %Odor timecourse
% subplot(1,4,3)
% hold on
% odor_dFFdervsld_density=odor_dFFdervsld_density/odor_no_points;
% max_ent=max(odor_dFFdervsld_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(odor_dFFdervsld_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-4    log10(max_ent)])
% xlabel('Derivative of lick freq (Hz/sec)')
% ylabel('dFF derivative (1/sec)');
% title(['Odor'])
% xlim([min(lick_derivatives) max(lick_derivatives)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% %Reinforcement timecourse
% subplot(1,4,4)
% hold on
% reinforcement_dFFdervsld_density=reinforcement_dFFdervsld_density/reinforcement_no_points;
% max_ent=max(reinforcement_dFFdervsld_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(reinforcement_dFFdervsld_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-4    log10(max_ent)])
% xlabel('Derivative of lick freq (Hz/sec)')
% ylabel('dFF derivative (1/sec)');
% title(['Reinforcement'])
% xlim([min(lick_derivatives) max(lick_derivatives)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% 
% %Plot the velocity vs dFF derivative plot
% %and save data for the pseudocolor plot
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% set(hFig, 'units','normalized','position',[.1 .25 .75 .25])
% 
% fprintf(1, ['\n\n\nProcessing rho and p value for dFF derivative vs velocity\n\n']);
% 
% %For the entire timecourse
% subplot(1,4,1)
% hold on
% 
% velocities=[];
% dFF_derivatives=[];
% 
% for trNo=1:handles_outs.no_lick_slopes
%     velocities=[velocities handles_outs.velocities(trNo,:)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,:)];
%     plot(handles_outs.velocities(trNo,:),handles_outs.dFF_derivatives(trNo,:),'.b')
% end
% 
% ymax=max(dFF_derivatives)+0.1*(max(dFF_derivatives)-min(dFF_derivatives));
% ymin=min(dFF_derivatives)-0.1*(max(dFF_derivatives)-min(dFF_derivatives));
% 
% xmax=max(velocities)+0.1*(max(velocities)-min(velocities));
% xmin=min(velocities)-0.1*(max(velocities)-min(velocities));
% 
% %Calculate the density for the pseudocolor
% no_bins=50;
% entire_dFFdervsVel_density=zeros(no_bins,no_bins);
% entire_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     entire_no_points=entire_no_points+length(velocities);
%     for ii=1:length(velocities)
%         entire_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             entire_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(velocities',dFF_derivatives','VariableNames',{'velocities','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~velocities');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(velocities',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs velocity  = %d, %d\n'],rho,pval);
% handles_outs.rhovel(1)=rho;
% handles_outs.pvalvel(1)=pval;
% xlabel('Velocity (au)')
% ylabel('dF/F derivative 1/sec')
% title('All times')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Now do the three windows
% time_licks=[];
% time_licks=handles_outs.time_licks;
% 
% %Before odor window
% subplot(1,4,2)
% hold on
% 
% time_mask=(time_licks>=-mean(delta_odor))&(time_licks<=0);
% velocities=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     velocities=[velocities handles_outs.velocities(trNo,time_mask)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
%     plot(handles_outs.velocities(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
% end
% 
% %Calculate the density for the pseudocolor
% pre_odor_dFFdervsVel_density=zeros(no_bins,no_bins);
% pre_odor_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     pre_odor_no_points=pre_odor_no_points+length(velocities);
%     for ii=1:length(velocities)
%         pre_odor_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             pre_odor_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(velocities',dFF_derivatives','VariableNames',{'velocities','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~velocities');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(velocities',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs velocity before odor = %d, %d\n'],rho,pval);
% handles_outs.rhovel(2)=rho;
% handles_outs.pvalvel(2)=pval;
% xlabel('Velocity (au)')
% ylabel('dF/F derivative 1/sec')
% title('Before odor')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Odor window
% subplot(1,4,3)
% hold on
% 
% 
% time_mask=(time_licks>=0)&(time_licks<=mean(delta_odor));
% velocities=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     velocities=[velocities handles_outs.velocities(trNo,time_mask)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
%     plot(handles_outs.velocities(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
% end
% 
% %Calculate the density for the pseudocolor
% odor_dFFdervsVel_density=zeros(no_bins,no_bins);
% odor_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     odor_no_points=odor_no_points+length(velocities);
%     for ii=1:length(velocities)
%         odor_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             odor_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(velocities',dFF_derivatives','VariableNames',{'velocities','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~velocities');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(velocities',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs velocity during odor = %d, %d\n'],rho,pval);
% handles_outs.rhovel(3)=rho;
% handles_outs.pvalvel(3)=pval;
% xlabel('Velocity (au)')
% ylabel('dF/F derivative 1/sec')
% title('During odor')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Reinforcement window
% subplot(1,4,4)
% hold on
% 
% time_mask=(time_licks>=mean(delta_odor_on_reinf_on))&(time_licks<=mean(delta_odor_on_reinf_on)+3);
% velocities=[];
% dFF_derivatives=[];
% for trNo=1:handles_outs.no_lick_slopes
%     velocities=[velocities handles_outs.velocities(trNo,time_mask)];
%     dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
%     plot(handles_outs.velocities(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
% end
% 
% %Calculate the density for the pseudocolor
% reinforcement_dFFdervsVel_density=zeros(no_bins,no_bins);
% reinforcement_no_points=0;
% for trNo=1:handles_outs.no_lick_slopes
%     reinforcement_no_points=reinforcement_no_points+length(velocities);
%     for ii=1:length(velocities)
%         reinforcement_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))=...
%             reinforcement_dFFdervsVel_density(fix(no_bins*(velocities(ii)-xmin)/(xmax-xmin)),fix(no_bins*(dFF_derivatives(ii)-ymin)/(ymax-ymin)))+1;
%     end
% end
% 
% tbl= table(velocities',dFF_derivatives','VariableNames',{'velocities','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~velocities');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(velocities',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs velocity during reinforcement = %d, %d\n'],rho,pval);
% handles_outs.rhovel(4)=rho;
% handles_outs.pvalvel(4)=pval;
% xlabel('Velocity (au)')
% ylabel('dF/F derivative 1/sec')
% title('Reinforcement')
% xlim([xmin xmax])
% ylim([ymin ymax])
% 
% %Now do a pseudocolor plot
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% set(hFig, 'units','normalized','position',[.1 .25 .75 .25])
% 
% %Entire timecourse
% subplot(1,4,1)
% hold on
% entire_dFFdervsVel_density=entire_dFFdervsVel_density/entire_no_points;
% max_ent=max(entire_dFFdervsVel_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(entire_dFFdervsVel_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-6    log10(max_ent)])
% xlabel('Velocity (au)')
% ylabel('dFF derivative (1/sec)');
% title(['Entire timecourse'])
% xlim([0 max(velocities)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% %Pre-odor timecourse
% subplot(1,4,2)
% hold on
% pre_odor_dFFdervsVel_density=pre_odor_dFFdervsVel_density/pre_odor_no_points;
% max_ent=max(pre_odor_dFFdervsVel_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(pre_odor_dFFdervsVel_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-6    log10(max_ent)])
% xlabel('Velocity (au)')
% ylabel('dFF derivative (1/sec)');
% title(['Pre-odor'])
% xlim([0 max(velocities)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% %Odor timecourse
% subplot(1,4,3)
% hold on
% odor_dFFdervsVel_density=odor_dFFdervsVel_density/odor_no_points;
% max_ent=max(odor_dFFdervsVel_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(odor_dFFdervsVel_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-6    log10(max_ent)])
% xlabel('Velocity (au)')
% ylabel('dFF derivative (1/sec)');
% title(['Odor'])
% xlim([0 max(velocities)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% %Pre-odor timecourse
% subplot(1,4,4)
% hold on
% reinforcement_dFFdervsVel_density=reinforcement_dFFdervsVel_density/reinforcement_no_points;
% max_ent=max(reinforcement_dFFdervsVel_density(:));
% min_ent=0.000000000000000000001;
% x_vel=([1:no_bins]*(xmax-xmin)/no_bins)+xmin-0.5*((xmax-xmin)/no_bins);
% y_dFF=([1:no_bins]*(ymax-ymin)/no_bins)+ymin-0.5*((ymax-ymin)/no_bins);
% drg_pcolor(repmat(x_vel,no_bins,1),repmat(y_dFF',1,no_bins),log10(reinforcement_dFFdervsVel_density'))
% colormap jet
% shading flat
% % min_prob=0.0113;
% % max_prob=0.0314;
% caxis([log10(max_ent)-6    log10(max_ent)])
% xlabel('Velocity (au)')
% ylabel('dFF derivative (1/sec)');
% title(['Reinforcement'])
% xlim([0 max(velocities)])
% ylim([min(dFF_derivatives) max(dFF_derivatives)])
% 
% handles_outs.delta_odor=delta_odor;
% handles_outs.delta_reinf=delta_reinf;
% handles_outs.delta_odor_on_reinf_on=delta_odor_on_reinf_on;

save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName(1:end-4) '_dims.mat'],'dim_out')


pffft=1;
