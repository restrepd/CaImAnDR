function drgCaImAn_add_trimmed_licks_to_pre_per(pre_perFileName,pre_perPathName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Load file
if exist('pre_perFileName')==0
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');
end
load([pre_perPathName pre_perFileName])

if exist('trimmed_licks')==0
    %Plot the event lines
    odor_on_times=[];
    ootii=0;

    %For S+ and S- plot odor on and reinforcement
    for epoch=1:handles.dropcData.epochIndex
        %Epoch 2 is odor on, 3 is odor off
        plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
        if plot_epoch
            % if show_figures==1
            %     if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
            %         plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
            %             '-r','LineWidth',1)
            %     else
            %         plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
            %             '-b','LineWidth',1)
            %     end
            % end
            if (handles.dropcData.epochEvent(epoch)==2)
                ootii=ootii+1;
                odor_on_times(ootii)=handles.dropcData.epochTime(epoch);
            end
        end


    end



    %Align the rhd times with the olfactometer

    %Find the FV, odor on and odor off events in digital_in recorded by INTAN
    ii=1;
    at_end=0;
    odor_on_times_rhd=[];
    FV_times_rhd=[];
    odor_off_times_rhd=[];
    iioon=0;
    iiFV=0;
    iiooff=0;
    digital_in=bitand(digital_in,2+4+8+16);
    while at_end==0
        ii_FV=find(digital_in(ii:end)==6,1,'first');
        if isempty(ii_FV)
            at_end=1;
        else
            %FV
            ii=ii+ii_FV-1;
            iiFV=iiFV+1;
            FV_times_rhd(iiFV)=ii/acq_rate;

            %Odor on
            ii_odor_on=find(digital_in(ii:end)==18,1,'first');
            %Odor off
            ii_odor_off=find(digital_in(ii:end)<18,1,'first');

            if (~isempty(ii_odor_on))&(~isempty(ii_odor_off))

                %Odor on
                ii=ii+ii_odor_on-1;
                iioon=iioon+1;
                odor_on_times_rhd(iioon)=ii/acq_rate;

                %Odor off

                ii=ii+ii_odor_off-1;
                iiooff=iiooff+1;
                odor_off_times_rhd(iiooff)=ii/acq_rate;

                ii=ii+1;
                if ii>=length(digital_in)
                    at_end=1;
                end
            else
                at_end=1;
            end
        end
    end

    %Find the alignment of the rhd vs the olfactometer times
    if length(odor_on_times)<length(odor_on_times_rhd)
        sum_delta=[];
        for ii=0:length(odor_on_times_rhd)-length(odor_on_times)
            sum_delta(ii+1)=abs(sum(odor_on_times_rhd(1+ii:ii+length(odor_on_times))-odor_on_times));
        end
        [min_del min_jj]=min(sum_delta);
        odor_on_times_rhd=odor_on_times_rhd(min_jj:min_jj+length(odor_on_times)-1);
    end

    if length(odor_on_times)>length(odor_on_times_rhd)
        sum_delta=[];
        for ii=0:length(odor_on_times)-length(odor_on_times_rhd)
            sum_delta(ii+1)=abs(sum(odor_on_times(1+ii:ii+length(odor_on_times_rhd))-odor_on_times_rhd));
        end
        [min_del min_jj]=min(sum_delta);
        odor_on_times=odor_on_times(min_jj:min_jj+length(odor_on_times_rhd)-1);
    end

    delta_t_rhd=mean(odor_on_times-odor_on_times_rhd);
    time_rhd=([1:length(adc_in)]/acq_rate)+delta_t_rhd;
    adc_in_trimmed=adc_in(time_rhd>0);
    time_rhd_trimmed=time_rhd(time_rhd>0);

    tic
    %Output licks in the time scale of traces
    trimmed_licks=zeros(1,length(time));
    for ii_tdec=1:length(time)
        t_from=time(ii_tdec)-dt/2;
        t_to=time(ii_tdec)+dt/2;
        if mean(adc_in_trimmed((time_rhd_trimmed>=t_from)&(time_rhd_trimmed<=t_to)))>1.5
            trimmed_licks(ii_tdec)=1;
        else
            trimmed_licks(ii_tdec)=0;
        end
    end

    toc



    save([pre_perPathName pre_perFileName],'per_ii','per_input_timecourse','per_targets',...
        'Hitii','Hitii_lick','no_Hit_traces','Hit_traces','which_trace_Hit','which_trial_Hit',...
        'Missii','Missii_lick','no_Miss_traces','Miss_traces','which_trace_Miss','which_trial_Miss',...
        'FAii','FAii_lick','no_FA_traces','FA_traces','which_trace_FA','which_trial_FA',...
        'CRii','CRii_lick','no_CR_traces','CR_traces','which_trace_CR','which_trial_CR',...
        'no_odor_trials','epoch_per_trial','epoch_time','time_to_event','no_traces',...
        'dt_before','delta_odor','delta_odor_on_reinf_on','delta_reinf',...
        'dt_after','dt_odor_onset','splus_traces','sminus_traces','sp_odor_response',...
        'sm_odor_response','CIsp','CIsm','time_to_eventCR','time_to_eventHit',...
        'time_to_eventMiss','time_to_eventFA','CICR','CIHit','CIMiss','CIFA',...
        'which_CRii_lick','which_FAii_lick','which_Hitii_lick','which_Missii_lick',...
        'CR_lick_times','FA_lick_times','Hit_lick_times','Miss_lick_times',...
        'CR_no_lick_times','FA_no_lick_times','Hit_no_lick_times','Miss_no_lick_times',...
        'time_licks','FAlick_freq','CRlick_freq','Hitlick_freq','Misslick_freq',...
        'handles','odor_traces','dt','all_lick_traces','acq_rate',...
        'y_shift','traces','time_rhd','adc_in','no_images','handles_out',...
        'lda_input_timecourse','lda_event','time_to_eventLDA','dHit_lick_traces'...
        ,'dCR_lick_traces','dMiss_lick_traces','dFA_lick_traces','time','epochs','digital_in',...
        'trimmed_licks')
end