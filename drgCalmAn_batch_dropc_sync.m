%% drgCaImAn_batch_dropc_sync.m
%
% Used for the synchronization control and the auditory signal analysis
%
close all
clear all



%Choices
do_warp=0;     %1=warped components from a reference file
um_per_pixel=1.745;

%Other default variables
plot_raw=1; %The default is raw
ref_win=[-5 -1.5];

ROIs=[];
% 1 dropc_nsampler_piriform

% 2 dropcspm_hf before 2/23/2018

% 3 dropcspm_hf after 2/24/2018
%
%  handles.dropcData.epochEvent
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR

dropc_program=3;

% Read choices file

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dropc_choices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles_choice.choiceFileName=choiceFileName;
handles_choice.choiceBatchPathName=choiceBatchPathName;

if ischar(handles_choice.PathName)
    cd(handles_choice.PathName)
end

try_catch_status=[];
for fileNo=handles_choice.first_file:handles_choice.no_files
    try_catch_status(fileNo).processed=1;
    try_catch_status(fileNo).saved=0;
end

for fileNo=handles_choice.first_file:handles_choice.no_files
     try
        figNo=0;

        fprintf(1, ['\nProcessing  ' num2str(fileNo) '\n']);

        dt=handles_choice.dt(fileNo);

        

        %Read tiff file
        imageData=imread_big([handles_choice.PathNamecsv{fileNo} handles_choice.csvFileName{fileNo}]);
        num_imgs=size(imageData,3);
        traces=zeros(num_imgs,1);
        for imNo=1:num_imgs
            this_image=zeros(size(imageData,1),size(imageData,2));
            this_image(:,:)=imageData(:,:,imNo);
            traces(imNo,1)=mean(this_image(:));
        end
        traces=traces';

        fnameca=handles_choice.csvFileName{fileNo};
        %     raw=[];
        %     inferred=[];
        %     try
        %         [raw,inferred]=drgGetCAtraces(Yr,A_or,C_or,b2,f2,Cn,options);
        %     catch
        %         pffft=1
        %     end

        %     % Should we use raw or inferred?
        %     if isfield(handles_choice,'plot_raw')
        %         plot_raw=handles_choice.plot_raw;
        %     end
        %
        %     if plot_raw==1
        %         traces=raw;
        %     else
        %         traces=inferred;
        %     end

        sz_traces=size(traces);
        no_traces=sz_traces(1);
        no_images=sz_traces(2);

        fprintf(1, ['\ndrgCaImAn_dropc run for ' handles_choice.csvFileName{fileNo} '\n\n']);

        %Read the dropc file
        handles=[];
        load([handles_choice.PathName{fileNo} handles_choice.spmFileName{fileNo}])

        %Read the rhd file
        adc_in=[];
        digital_in=[];
        acq_rate=[];
        [adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file([handles_choice.PathName{fileNo} handles_choice.rhdFileName{fileNo}],3);


        % Plot the traces
        these_lines{1}='-b';
        these_lines{2}='-r';
        these_lines{3}='-m';
        these_lines{8}='-g';
        these_lines{5}='-y';
        these_lines{6}='-k';
        these_lines{7}='-c';
        these_lines{4}='-k';

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .85 .25])


        hold on

        % Determine the y spacing of the traces
        y_shift=1.2*prctile(traces(:),98);

        %Plot the event lines
        odor_on_times=[];
        ootii=0;
        switch dropc_program
            case 1
                for event=2:handles.dropcData.eventIndex
                    plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_traces+2)*y_shift],...
                        these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
                end
            case 2
                for event=1:handles.dropcData.allTrialIndex
                    plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_traces+2)*y_shift],...
                        these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                end
            case 3
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
                        if (handles.dropcData.epochEvent(epoch)==2)
                            ootii=ootii+1;
                            odor_on_times(ootii)=handles.dropcData.epochTime(epoch);
                        end
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

        %Plot the licks recorded by the INTAN (adc_in)
        time_rhd=([1:length(adc_in)]/acq_rate)+delta_t_rhd;
        pct998=prctile(adc_in,99.8);
        pct1=prctile(adc_in,1);
        norm_fact=0.8*y_shift/(pct998-pct1);

        plot(time_rhd(time_rhd>0),adc_in(time_rhd>0)*norm_fact)

        %Plot the traces
        time=[1:no_images]*dt;
        for trNo=1:no_traces
            % for trNo=1:20
            plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        end
 
        %This vector will save the epochs (Hit, CR, etc)
        epochs=zeros(1,length(time));


        ylim([-y_shift*0.2 1.2*(prctile(traces(:),98)+y_shift)])
 

        xlabel('time (s)')
        ylabel('deltaF/F')
        title(fnameca(1:end-4), 'Interpreter', 'none')
 
        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig1.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig1.fig'])
        end
        
        

        dt_before=10;
        dt_after=20;

        %Uncomment this to show a subset of traces
        %This code is here to explore the individual traces and licks

        %     first_trace=1;
        %     last_trace=2.5;
        %     trialNo_start=1;
        %     trialNo_end=12;
        %     xlim([odor_on_times(trialNo_start)-dt_before odor_on_times(trialNo_end)+dt_after])
        %     ylim([y_shift*(first_trace-1) y_shift*last_trace])
        %     title(['deltaF/f for traces ' num2str(first_trace) ' to ' num2str(last_trace) ])
%         set(hFig, 'units','normalized','position',[.05 .05 .3 .85])



        %Plot the responses aligned with the onset of the epochs
        switch dropc_program
            case 3
                dt_odor_onset=0.1085;  %This is the time from FV off to odor entering the nose cone

                timesSD=5;
                timesSDodorOn=2.5;




                response_points=1*dt_after;
                odor_response_points=1;
                dt_trace=2000000;
                dt_rhd_trace=1000000000000;


                %Determine the odor on and reinforcement times
                iido=0;
                delta_odor=[];
                iidr=0;
                delta_reinf=[];
                iidro=0;
                delta_odor_on_reinf_on=[];

                for epoch=1:handles.dropcData.epochIndex
                    if ((handles.dropcData.epochTime(epoch)-dt_before)>0)&(handles.dropcData.epochTime(epoch)+dt_after<=max(time))

                        if (handles.dropcData.epochEvent(epoch)==4)
                            %This is a reinforcement on, find the odor on and
                            %the reinforcement off

                            %Find reinforcement off
                            next_epoch=epoch+1;
                            while handles.dropcData.epochEvent(next_epoch)~=5
                                next_epoch=next_epoch+1;
                            end
                            iidr=iidr+1;
                            delta_reinf(iidr)=handles.dropcData.epochTime(next_epoch)-handles.dropcData.epochTime(epoch);

                            %Find odor on
                            next_epoch=epoch-1;
                            while handles.dropcData.epochEvent(next_epoch)~=2
                                next_epoch=next_epoch-1;
                            end
                            iidro=iidro+1;
                            delta_odor_on_reinf_on(iidro)=handles.dropcData.epochTime(epoch)-handles.dropcData.epochTime(next_epoch);
                        end

                        if (handles.dropcData.epochEvent(epoch)==2)
                            %This is an odor on event, find the next odor off
                            next_epoch=epoch+1;
                            while handles.dropcData.epochEvent(next_epoch)~=3
                                next_epoch=next_epoch+1;
                            end
                            iido=iido+1;
                            delta_odor(iido)=handles.dropcData.epochTime(next_epoch)-handles.dropcData.epochTime(epoch);
                        end
                    end
                end


                %Extract all licks and determine the threshold
                all_lick_traces=[];
                allii_lick=0;

                for epoch=1:handles.dropcData.epochIndex
                    if ((handles.dropcData.epochTime(epoch)-dt_before)>0)&(handles.dropcData.epochTime(epoch)+dt_after<=max(time))
                        if (handles.dropcData.epochEvent(epoch)==2)
                            %Odor on
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            allii_lick=allii_lick+1;
                            all_lick_traces(allii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        end
                    end
                end
                threshold_lick=prctile(all_lick_traces(:),1)+((prctile(all_lick_traces(:),99)-prctile(all_lick_traces(:),1))/2);

                %Initialize variables

                odor_traces=[];
                od_ii=0;

                fvii=0;
                fvii_lick=0;
                no_fv_traces=0;
                fv_traces=[];

                sp_odor_response=[];
                splus_traces=[];
                spii=0;
                spii_lick=0;
                smii=0;
                smii_lick=0;
                no_sp_traces=0;
                no_sm_traces=0;
                splus_traces=[];
                sminus_traces=[];
                sm_odor_response=[];
                which_trace_sp=[];
                which_trace_sm=[];
                splus_lick_traces=[];
                sminus_lick_traces=[];


                %lda input
                lda_input_timecourse=[];
                lda_event=[];

                %Perceptron input
                per_ii=0;
                per_input=[];
                per_input_timecourse=[];

                %Diefferent epochs
                Hitii=0;
                Hitii_lick=0;
                no_Hit_traces=0;
                Hit_traces=[];
                which_trace_Hit=[];
                which_trial_Hit=[];
                which_Hitii_lick=[];
                Hit_lick_times=[];
                Hit_no_lick_times=0;
                no_Hit_trials=0;

                Missii=0;
                Missii_lick=0;
                no_Miss_traces=0;
                Miss_traces=[];
                which_trace_Miss=[];
                which_trial_Miss=[];
                which_Missii_lick=[];
                Miss_lick_times=[];
                Miss_no_lick_times=0;

                FAii=0;
                FAii_lick=0;
                no_FA_traces=0;
                FA_traces=[];
                which_trace_FA=[];
                which_trial_FA=[];
                which_FAii_lick=[];
                FA_lick_times=[];
                FA_no_lick_times=0;


                CRii=0;
                CRii_lick=0;
                no_CR_traces=0;
                CR_traces=[];
                which_trace_CR=[];
                which_trial_CR=[];
                which_CRii_lick=[];
                CR_lick_times=[];
                CR_no_lick_times=0;

                per_targets=[];
                per_which_events=[];


                no_odor_trials=0;
                no_spm_odor_trials=0;
                epoch_per_trial=[];
                epoch_time=[];

                %Find Hits, CRs, etc
                for epoch=1:handles.dropcData.epochIndex
                    if ((handles.dropcData.epochTime(epoch)-dt_before)>0)&(handles.dropcData.epochTime(epoch)+dt_after<=max(time))
                        %This event is processed for Ca

                        %Final valve epoch
                        if (handles.dropcData.epochEvent(epoch)==1)

                            if handles.dropcData.epochEvent(epoch+1)==2
                                fv_mask=(time>=handles.dropcData.epochTime(epoch))...
                                    &(time<=handles.dropcData.epochTime(epoch+1));
                                epochs(fv_mask)=1;
                            end

                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after);
                            ref_mask=(time>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(time<=handles.dropcData.epochTime(epoch)+ref_win(2));

                            %FV on
                            %Exclude the first snip if it is too close to the start
                            handles_out.no_components=no_traces;
                            trace_num=0;
                            for trNo=1:no_traces
                                if sum(isinf(traces(trNo,snip_mask)))==0
                                    trace_num=trace_num+1;
                                    fvii=fvii+1;
                                    fv_traces(fvii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo,ref_mask));
                                    handles_out.componentNo(trNo).trialNo(fvii_lick+1).fv_traces=fv_traces(fvii,1:sum(snip_mask));
                                    dt_trace=min([dt_trace dt*sum(snip_mask)]);
                                    no_fv_traces=no_fv_traces+1;
                                end
                            end
                            handles_out.no_components_fv=trace_num;

                            %Get the licks
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after);
                            fvii_lick=fvii_lick+1;
                            fv_lick_traces(fvii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);
                        end
                        handles_out.no_fv_trials=fvii_lick;

                        %Now do S+ and S-
                        if (handles.dropcData.epochEvent(epoch)==2)
                            %Odor on
                            odor_mask=(time>=handles.dropcData.epochTime(epoch))...
                                &(time<=handles.dropcData.epochTime(epoch+1));
                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            ref_mask=(time>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(time<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            per_ii=per_ii+1;
                            row_per=0;
                            tr_ii=0;
                            for trNo=1:no_traces
                                %                         %Save data for perceptron analysis in this script
                                %                         this_trace=[];
                                %                         this_trace=traces(trNo, perceptron_mask);
                                %                         %if sum(this_trace(floor(dt_before/dt)+1:end)>mean(this_trace(1:floor(dt_before/dt)))+timesSD*std(this_trace(1:floor(dt_before/dt))))>=response_points
                                %                         per_input(row_per+1:row_per+length(this_trace),per_ii)=this_trace;
                                %                         row_per=row_per+length(this_trace);

                                %Save data for comprehensive perceptron analysis
                                if sum(isinf(traces(trNo,snip_mask)))==0
                                    tr_ii=tr_ii+1;
                                    per_input_timecourse(1:sum(snip_mask),tr_ii,per_ii)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                end
                            end

                            if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor

                                %S plus
                                %Assign variables for perceptron analysis
                                per_which_events(1,per_ii)=1;
                                per_targets(1,per_ii)=1;
                                per_targets(2,per_ii)=0;

                                no_spm_odor_trials=no_spm_odor_trials+1;
                                handles_out.trialNo(smii_lick+1).trianNo=no_spm_odor_trials;
                                handles_out.no_spm_odor_trials=no_spm_odor_trials;

                                trace_num=0;
                                for trNo=1:no_traces
                                    if sum(isinf(traces(trNo,snip_mask)))==0
                                        trace_num=trace_num+1;
                                        this_trace=[];
                                        this_trace=traces(trNo,snip_mask);
                                        %if sum(this_trace(floor(dt_before/dt)+1:end)>mean(this_trace(1:floor(dt_before/dt)))+timesSD*std(this_trace(1:floor(dt_before/dt))))>=response_points
                                        spii=spii+1;
                                        splus_traces(spii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                        handles_out.componentNo(trace_num).trialNo(spii_lick+1).splus_traces=splus_traces(spii,1:sum(snip_mask));
                                        od_ii=od_ii+1;
                                        odor_traces(od_ii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                        dt_trace=min([dt_trace dt*sum(snip_mask)]);
                                        if sum(this_trace(floor((dt_before-1)/dt)+1:floor(dt_before/dt)+1)...
                                                >mean(this_trace(1:floor((dt_before-1)/dt)))+timesSDodorOn*std(this_trace(1:floor((dt_before-1)/dt))))>=odor_response_points
                                            sp_odor_response(spii)=0;
                                        else
                                            sp_odor_response(spii)=1;
                                        end
                                        no_sp_traces=no_sp_traces+1;
                                        which_trace_sp(spii)=trNo;
                                    end
                                end
                                handles_out.no_components_sp=trace_num;
                                rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                                spii_lick=spii_lick+1;
                                splus_lick_traces(spii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                                dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);

                                handles_out.no_sp_trials=spii_lick;
                            else

                                %S minus
                                %Assign variables for perceptron analysis
                                per_which_events(2,per_ii)=1;
                                per_targets(1,per_ii)=0;
                                per_targets(2,per_ii)=1;

                                no_spm_odor_trials=no_spm_odor_trials+1;
                                handles_out.trialNo(smii_lick+1).trianNo=no_spm_odor_trials;
                                handles_out.no_spm_odor_trials=no_spm_odor_trials;

                                trace_num=0;
                                for trNo=1:no_traces
                                    if sum(isinf(traces(trNo,snip_mask)))==0
                                        trace_num=trace_num+1;
                                        this_trace=[];
                                        this_trace=traces(trNo,snip_mask);
                                        %if sum(this_trace(floor(dt_before/dt)+1:end)>mean(this_trace(1:floor(dt_before/dt)))+timesSD*std(this_trace(1:floor(dt_before/dt))))>=response_points
                                        smii=smii+1;
                                        sminus_traces(smii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                        handles_out.componentNo(trace_num).trialNo(smii_lick+1).splus_traces=sminus_traces(smii,1:sum(snip_mask));
                                        od_ii=od_ii+1;
                                        odor_traces(od_ii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                        dt_trace=min([dt_trace dt*sum(snip_mask)]);

                                        if sum(this_trace(floor((dt_before-1)/dt)+1:floor(dt_before/dt)+1)...
                                                >mean(this_trace(1:floor((dt_before-1)/dt)))+timesSDodorOn*std(this_trace(1:floor(dt_before/dt))))>=odor_response_points
                                            sm_odor_response(smii)=0;
                                        else
                                            sm_odor_response(smii)=1;
                                        end

                                        no_sm_traces=no_sm_traces+1;
                                        which_trace_sm(smii)=trNo;
                                    end
                                end
                                rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                                smii_lick=smii_lick+1;
                                sminus_lick_traces(smii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                                dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);

                                handles_out.no_sm_trials=smii_lick;
                            end
                            %end
                        end

                        %Find Hit, CR, FA and Miss

                        %Hit
                        if (handles.dropcData.epochEvent(epoch)==6)

                            no_odor_trials=no_odor_trials+1;
                            no_Hit_trials=no_Hit_trials+1;
                            valid_trace(no_odor_trials)=1;
                            epoch_per_trial(no_odor_trials)=6;
                            epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S+';

                            handles_out.epoch_per_trial(no_odor_trials)=6;

                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            epochs(odor_mask)=6;
                            ref_mask=(time>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(time<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            %Each trace is from a different neuron
                            trace_num=0;
                            for trNo=1:no_traces
                                if sum(isinf(traces(trNo,snip_mask)))==0
                                    trace_num=trace_num+1;
                                    this_trace=[];
                                    this_trace=traces(trNo,snip_mask);
                                    Hitii=Hitii+1;
                                    Hit_traces(Hitii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                    handles_out.componentNo(trace_num).trialNo(Hitii_lick+1).hit_traces=Hit_traces(Hitii,1:sum(snip_mask));
                                    handles_out.file(fileNo).trial_this_file(no_odor_trials).componentNo(trace_num).traces=Hit_traces(Hitii,1:sum(snip_mask));
                                    handles_out.trialNo(Hitii_lick+1).trace_numHit=trace_num;
                                    no_Hit_traces=no_Hit_traces+1;
                                    Hit_trials(Hitii)=no_Hit_trials;
                                    which_trace_Hit(Hitii)=trNo;
                                    which_trial_Hit(Hitii)=no_odor_trials;
                                    which_Hitii_lick(Hitii)=Hitii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                end
                            end



                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            Hitii_lick=Hitii_lick+1;
                            Hit_lick_traces(Hitii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);

                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            Hit_lick_times(Hitii_lick,1:length(these_lick_times))=these_lick_times;
                            Hit_no_lick_times(Hitii_lick)=length(these_lick_times);


                            handles_out.no_Hit_trials=Hitii_lick;
                            handles_out.Hit_trial_no(no_odor_trials)=Hitii_lick;
                        end

                        %Miss
                        if (handles.dropcData.epochEvent(epoch)==7)

                            no_odor_trials=no_odor_trials+1;
                            valid_trace(no_odor_trials)=1;
                            epoch_per_trial(no_odor_trials)=7;
                            epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S+';

                            handles_out.epoch_per_trial(no_odor_trials)=7;

                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            epochs(odor_mask)=7;
                            ref_mask=(time>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(time<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            trace_num=0;
                            for trNo=1:no_traces
                                if sum(isinf(traces(trNo,snip_mask)))==0
                                    trace_num=trace_num+1;
                                    this_trace=[];
                                    this_trace=traces(trNo,snip_mask);
                                    Missii=Missii+1;
                                    Miss_traces(Missii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                    handles_out.componentNo(trace_num).trialNo(Missii_lick+1).miss_traces=Miss_traces(Missii,1:sum(snip_mask));
                                    handles_out.file(fileNo).trial_this_file(no_odor_trials).componentNo(trace_num).traces=Miss_traces(Missii,1:sum(snip_mask));
                                    handles_out.trialNo(Missii_lick+1).trace_numMiss=trace_num;
                                    no_Miss_traces=no_Miss_traces+1;
                                    which_trace_Miss(Missii)=trNo;
                                    which_trial_Miss(Missii)=no_odor_trials;
                                    which_Missii_lick(Missii)=Missii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                end
                            end
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            Missii_lick=Missii_lick+1;
                            Miss_lick_traces(Missii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);

                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            Miss_lick_times(Missii_lick,1:length(these_lick_times))=these_lick_times;
                            Miss_no_lick_times(Missii_lick)=length(these_lick_times);

                            handles_out.no_Miss_trials=Missii_lick;
                            handles_out.Miss_trial_no(no_odor_trials)=Missii_lick;
                        end

                        %FA
                        if (handles.dropcData.epochEvent(epoch)==8)

                            no_odor_trials=no_odor_trials+1;
                            valid_trace(no_odor_trials)=1;
                            epoch_per_trial(no_odor_trials)=8;
                            epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S-';

                            handles_out.epoch_per_trial(no_odor_trials)=8;

                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            epochs(odor_mask)=8;
                            ref_mask=(time>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(time<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            trace_num=0;
                            for trNo=1:no_traces
                                if sum(isinf(traces(trNo,snip_mask)))==0
                                    trace_num=trace_num+1;
                                    this_trace=[];
                                    this_trace=traces(trNo,snip_mask);
                                    FAii=FAii+1;
                                    FA_traces(FAii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                    handles_out.componentNo(trace_num).trialNo(FAii_lick+1).FA_traces=FA_traces(FAii,1:sum(snip_mask));
                                    handles_out.file(fileNo).trial_this_file(no_odor_trials).componentNo(trace_num).traces=FA_traces(FAii,1:sum(snip_mask));
                                    handles_out.trialNo(FAii_lick+1).trace_numFA=trace_num;
                                    no_FA_traces=no_FA_traces+1;
                                    which_trace_FA(FAii)=trNo;
                                    which_trial_FA(FAii)=no_odor_trials;
                                    which_FAii_lick(FAii)=FAii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                end
                            end
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            FAii_lick=FAii_lick+1;
                            FA_lick_traces(FAii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);

                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            FA_lick_times(FAii_lick,1:length(these_lick_times))=these_lick_times;
                            FA_no_lick_times(FAii_lick)=length(these_lick_times);

                            handles_out.no_FA_trials=FAii_lick;
                            handles_out.FA_trial_no(no_odor_trials)=FAii_lick;
                        end

                        %CR
                        if (handles.dropcData.epochEvent(epoch)==9)

                            no_odor_trials=no_odor_trials+1;
                            valid_trace(no_odor_trials)=1;
                            epoch_per_trial(no_odor_trials)=9;
                            epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S-';

                            handles_out.epoch_per_trial(no_odor_trials)=9;

                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            epochs(odor_mask)=9;
                            ref_mask=(time>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(time<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            trace_num=0;
                            for trNo=1:no_traces
                                if sum(isinf(traces(trNo,snip_mask)))==0
                                    trace_num=trace_num+1;
                                    this_trace=[];
                                    this_trace=traces(trNo,snip_mask);
                                    CRii=CRii+1;
                                    CR_traces(CRii,1:sum(snip_mask))=traces(trNo,snip_mask)-mean(traces(trNo, ref_mask));
                                    handles_out.componentNo(trace_num).trialNo(CRii_lick+1).CR_traces=CR_traces(CRii,1:sum(snip_mask));
                                    handles_out.file(fileNo).trial_this_file(no_odor_trials).componentNo(trace_num).traces=CR_traces(CRii,1:sum(snip_mask));
                                    handles_out.trialNo(CRii_lick+1).trace_numCR=trace_num;
                                    no_CR_traces=no_CR_traces+1;
                                    which_trace_CR(CRii)=trNo;
                                    which_trial_CR(CRii)=no_odor_trials;
                                    which_CRii_lick(CRii)=CRii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                end
                            end
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            CRii_lick=CRii_lick+1;
                            CR_lick_traces(CRii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);

                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            CR_lick_times(CRii_lick,1:length(these_lick_times))=these_lick_times;
                            CR_no_lick_times(CRii_lick)=length(these_lick_times);

                            handles_out.no_CR_trials=CRii_lick;
                            handles_out.CR_trial_no(no_odor_trials)=CRii_lick;


                        end

                        %reinforcement
                        if (handles.dropcData.epochEvent(epoch)==4)
                            reinforcement_mask=(time>=handles.dropcData.epochTime(epoch))...
                                &(time<=handles.dropcData.epochTime(epoch+1));
                            epochs(reinforcement_mask)=4;
                        end


                    end
                end
        end

        handles_out.file(fileNo).no_odor_trials=no_odor_trials;

        %Plot the snips
        szFV=size(fv_traces);
        time_to_event=([1:szFV(2)]*dt-dt_before);
        szFVl=size(fv_lick_traces);
        time_to_event_l=([1:szFVl(2)/20]*(20/acq_rate)-dt_before);
        handles_out.time_to_eventFV=time_to_event;

        %Calculate lick frequency
        dt_lick=0.3;
        lick_t_start=0;
        lick_t_end=3;

        %Hit
        Hitlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        Hitlick_per_trial_timecourse=zeros(Hitii_lick,ceil((dt_after+dt_before)/dt_lick));

        for ii=1:Hitii_lick
            for ii_lick=1:Hit_no_lick_times(ii)
                Hitlick_freq( ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))=Hitlick_freq( ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
                Hitlick_per_trial_timecourse(ii,ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))=Hitlick_per_trial_timecourse(ii,ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
            these_times=zeros(1,Hit_no_lick_times(ii));
            these_times(1,1:Hit_no_lick_times(ii))=Hit_lick_times(ii,1:Hit_no_lick_times(ii));
            Hitlick_per_trial(ii)=sum((these_times>=lick_t_start)&(these_times<=lick_t_end))/(lick_t_end-lick_t_start);
        end

        Hitlick_freq=(Hitlick_freq/(Hitii_lick*dt_lick));
        Hitlick_per_trial_timecourse=Hitlick_per_trial_timecourse/dt_lick;

        %Miss
        Misslick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));


        for ii=1:Missii_lick
            for ii_lick=1:Miss_no_lick_times(ii)
                Misslick_freq( ceil((Miss_lick_times(ii,ii_lick)+dt_before)/dt_lick))=Misslick_freq( ceil((Miss_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
        end

        Misslick_freq=(Misslick_freq/(Missii_lick*dt_lick));


        %CR
        CRlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));


        for ii=1:CRii_lick
            for ii_lick=1:CR_no_lick_times(ii)
                CRlick_freq( ceil((CR_lick_times(ii,ii_lick)+dt_before)/dt_lick))=CRlick_freq( ceil((CR_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
        end

        CRlick_freq=(CRlick_freq/(CRii_lick*dt_lick));

        %FA
        FAlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));


        for ii=1:FAii_lick
            for ii_lick=1:FA_no_lick_times(ii)
                FAlick_freq( ceil((FA_lick_times(ii,ii_lick)+dt_before)/dt_lick))=FAlick_freq( ceil((FA_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
        end

        FAlick_freq=(FAlick_freq/(FAii_lick*dt_lick));

        time_licks=([1:length(Hitlick_freq)]*dt_lick)-dt_before+dt_lick/2;

 
        %FV
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        hold on
        CI = bootci(1000, @mean, fv_traces);
        meanfv=mean(fv_traces,1);
        CI(1,:)=meanfv-CI(1,:);
        CI(2,:)=CI(2,:)-meanfv;
        [hl1, hp1] = boundedline(time_to_event',mean(fv_traces,1)', CI', 'r');
        plot([0 0],[0 max(mean(fv_traces,1)')+max(CI(:))],'-k')
        pct5=prctile(mean(fv_traces,1),5);
        pct95=prctile(mean(fv_traces,1),95);
        ylim([pct5-0.2*(pct95-pct5) pct95+0.2*(pct95-pct5)])
        xlabel('Time (sec)')
        ylabel('dF/F')
        title('Ca changes aligned to final valve diversion')

        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig2.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig2.fig'])
        end


        %S+, S-, all snips
        CIsm = bootci(1000, @mean, sminus_traces);
        meansm=mean(sminus_traces,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;

        CIsp = bootci(1000, @mean, splus_traces);
        meansp=mean(splus_traces,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;


        for ii=1:smii_lick
            dsminus_lick_traces(ii,:)=decimate(sminus_lick_traces(ii,:)',20)';
        end
        pct99=max(dsminus_lick_traces(:));
        pct1=prctile(dsminus_lick_traces(:),1);
        dsminus_lick_traces=(dsminus_lick_traces/(pct99-pct1))+4;
        CIsmlick = bootci(1000, @mean, dsminus_lick_traces);
        meansmlick=mean(dsminus_lick_traces,1);
        CIsmlick(1,:)=meansmlick-CIsmlick(1,:);
        CIsmlick(2,:)=CIsmlick(2,:)-meansmlick;

        for ii=1:spii_lick
            dsplus_lick_traces(ii,:)=decimate(splus_lick_traces(ii,:)',20)';
        end
        pct99=prctile(dsplus_lick_traces(:),99);
        pct1=prctile(dsplus_lick_traces(:),1);
        dsplus_lick_traces=(dsplus_lick_traces/(pct99-pct1))+5.5;
        CIsplick = bootci(1000, @mean, dsplus_lick_traces);
        meansplick=mean(dsplus_lick_traces,1);
        CIsplick(1,:)=meansplick-CIsplick(1,:);
        CIsplick(2,:)=CIsplick(2,:)-meansplick;


        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        hold on

        szSp=size(splus_traces);
        szSm=size(sminus_traces);
        %time_to_event=([1:szSm(2)]*dt-dt_before);
        time_to_eventSm=([1:szSm(2)]*dt-dt_before);
        time_to_eventSp=([1:szSp(2)]*dt-dt_before);
        handles_out.time_to_eventSm=time_to_eventSm;
        handles_out.time_to_eventSp=time_to_eventSp;

        pct1=prctile([mean(sminus_traces,1)'; mean(splus_traces(:,1:szSp(2)),1)'],1);
        pct99=prctile([mean(sminus_traces,1)'; mean(splus_traces(:,1:szSp(2)),1)'],99);

 

        [hlsm, hpsm] = boundedline(time_to_eventSm',mean(sminus_traces,1)', CIsm', 'b');
        [hlsp, hpsp] = boundedline(time_to_eventSp',mean(splus_traces,1)', CIsp', 'r');

        %Odor on markers
        plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
        odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')

        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
 
              
       
        title("Ca changes aligned to odor onset")
        legend([hlsp hlsm odorhl reinfhl],'S+','S-','Odor','Reinforcement')
        xlabel('Time (sec)')
        ylabel('dF/F')
        ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
        xlim([-10 19.8])
   
        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig3.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig3.fig'])
        end

        %Hit, CR, et al
        CIHit=[];
        CIMiss=[];
        CIFA=[];
        CICR=[];
        time_to_eventHit=[];
        time_to_eventCR=[];
        time_to_eventFA=[];
        time_to_eventMiss=[];
        min_no_time_points=20000000;

        %S+, S-, all snips
        if no_Hit_traces>1
            meanHit=mean(Hit_traces,1);
            time_to_eventHit=([1:length(meanHit)]*dt-dt_before);
            handles_out.time_to_eventHit=time_to_eventHit;
            %The snip mask sometimes differ by one point for Hit, Miss, CR, FA
            min_no_time_points=min([min_no_time_points length(time_to_eventHit)]);
        end

        if no_Hit_traces>2
            CIHit = bootci(1000, @mean, Hit_traces);
            CIHit(1,:)=meanHit-CIHit(1,:);
            CIHit(2,:)=CIHit(2,:)-meanHit;
        end


        if no_Miss_traces>1
            meanMiss=mean(Miss_traces,1);
            time_to_eventMiss=([1:length(meanMiss)]*dt-dt_before);
            handles_out.time_to_eventMiss=time_to_eventMiss;
            %The snip mask sometimes differ by one point for Hit, Miss, CR, FA
            min_no_time_points=min([min_no_time_points length(time_to_eventMiss)]);
        end

        if no_Miss_traces>2
            CIMiss = bootci(1000, @mean, Miss_traces);
            CIMiss(1,:)=meanMiss-CIMiss(1,:);
            CIMiss(2,:)=CIMiss(2,:)-meanMiss;
        end


        if no_CR_traces>1
            meanCR=mean(CR_traces,1);
            time_to_eventCR=([1:length(meanCR)]*dt-dt_before);
            handles_out.time_to_eventCR=time_to_eventCR;
            %The snip mask sometimes differ by one point for Hit, Miss, CR, FA
            min_no_time_points=min([min_no_time_points length(time_to_eventCR)]);
        end

        if no_CR_traces>2
            CICR = bootci(1000, @mean, CR_traces);
            CICR(1,:)=meanCR-CICR(1,:);
            CICR(2,:)=CICR(2,:)-meanCR;
        end

        if no_FA_traces>1
            meanFA=mean(FA_traces,1);
            time_to_eventFA=([1:length(meanFA)]*dt-dt_before);
            handles_out.time_to_eventFA=time_to_eventFA;
            %The snip mask sometimes differ by one point for Hit, Miss, CR, FA
            min_no_time_points=min([min_no_time_points length(time_to_eventFA)]);
        end


        if no_FA_traces>2
            CIFA = bootci(1000, @mean, FA_traces);
            CIFA(1,:)=meanFA-CIFA(1,:);
            CIFA(2,:)=CIFA(2,:)-meanFA;
        end

        lda_input_timecourse=lda_input_timecourse(1:min_no_time_points,:,:);
        time_to_eventLDA=([1:min_no_time_points]*dt-dt_before);

        %Decimate the licks
        try
            %CR lick traces
            dCR_lick_traces=[];
            for ii=1:CRii_lick
                dCR_lick_traces(ii,:)=decimate(CR_lick_traces(ii,:)',20)';
            end

            %FA lick traces
            dFA_lick_traces=[];
            for ii=1:FAii_lick
                dFA_lick_traces(ii,:)=decimate(FA_lick_traces(ii,:)',20)';
            end

            %Miss lick traces
            dMiss_lick_traces=[];
            for ii=1:Missii_lick
                dMiss_lick_traces(ii,:)=decimate(Miss_lick_traces(ii,:)',20)';
            end

            %Hit lick traces
            dHit_lick_traces=[];
            for ii=1:Hitii_lick
                dHit_lick_traces(ii,:)=decimate(Hit_lick_traces(ii,:)',20)';
            end

        catch
        end
        %Save the calculated data
        if do_warp==1
            save_name=[handles_choice.PathNamecsv{fileNo} fnameca(1:end-4) '_warpv2_pre_per.mat'];
        else
            save_name=[handles_choice.PathNamecsv{fileNo} fnameca(1:end-4) '_batchv2_pre_per.mat']
        end

        %splus_traces(handles_out.no_sp_trials*no_traces,no_timepoints), order is  for trial 1: ROI1, ROI2, etc.. ROIend, trial 2:

        save(save_name,'per_ii','per_input_timecourse','per_targets',...
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
            ,'dCR_lick_traces','dMiss_lick_traces','dFA_lick_traces','time','epochs','digital_in')

        try_catch_status(fileNo).saved=1;

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        hold on

        pct1=prctile([mean(CR_traces,1)'; mean(Hit_traces,1)';mean(Miss_traces,1)';mean(FA_traces,1)'],1);
        pct99=prctile([mean(CR_traces,1)'; mean(Hit_traces,1)';mean(Miss_traces,1)';mean(FA_traces,1)'],99);

        %Odor on markers
        plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
        odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')

        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0.3 0.3],'-r','LineWidth',5);

        if no_CR_traces>2
            [hlCR, hpCR] = boundedline(time_to_eventCR',mean(CR_traces,1)', CICR', 'b');
        end
        if no_Hit_traces>2
            [hlHit, hpHit] = boundedline(time_to_eventHit',mean(Hit_traces,1)', CIHit', 'r');
        end
        if no_Miss_traces>2
            [hlMiss, hpMiss] = boundedline(time_to_eventMiss',mean(Miss_traces,1)', CIMiss', 'c');
        end
        if no_FA_traces>2
            [hlFA, hpFA] = boundedline(time_to_eventFA',mean(FA_traces,1)', CIFA', 'm');
        end

        title("Ca changes aligned to odor onset")

        if (no_Hit_traces>2)&(no_CR_traces>2)&(no_FA_traces>2)&(no_Miss_traces>2)
            legend([hlHit hlMiss hlCR hlFA],'Hit','Miss','CR','FA')
        end
        if (no_Hit_traces>2)&(no_CR_traces>2)&(no_FA_traces<=2)&(no_Miss_traces<=2)
            legend([hlHit hlCR],'Hit','CR')
        end
        if (no_Hit_traces>2)&(no_CR_traces>2)&(no_FA_traces>2)&(no_Miss_traces<=2)
            legend([hlHit hlCR hlFA],'Hit','CR','FA')
        end
        if (no_Hit_traces>2)&(no_CR_traces>2)&(no_FA_traces<=2)&(no_Miss_traces>2)
            legend([hlHit hlMiss hlCR],'Hit','Miss','CR')
        end

        xlabel('Time (sec)')
        ylabel('dF/F')
        ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
        xlim([-10 19.8])

        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig4.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig4.fig'])
        end
 
        %Plot the licks
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        hold on

        for ii=1:allii_lick
            dall_lick_traces(ii,:)=decimate(all_lick_traces(ii,:)',20)';
        end
        szalllick=size(dall_lick_traces);
        time_licks=([1:szalllick(2)]/(acq_rate/20))-dt_before;
        per99=prctile(dall_lick_traces(:),99.9);
        per1=prctile(dall_lick_traces(:),1);

        mean_Hit_licks=zeros(1,szalllick(2));
        mean_Miss_licks=zeros(1,szalllick(2));
        mean_FA_licks=zeros(1,szalllick(2));
        mean_CR_licks=zeros(1,szalllick(2));

        y_shift=0;

        %Plot CR lick traces
        for ii=1:CRii_lick
            dCR_lick_traces(ii,:)=decimate(CR_lick_traces(ii,:)',20)';
            plot(time_licks,dCR_lick_traces(ii,:)+y_shift,'-b')
            mean_CR_licks(1,:)= mean_CR_licks(1,:) + ((dCR_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end

        %Plot FA lick traces
        for ii=1:FAii_lick
            dFA_lick_traces(ii,:)=decimate(FA_lick_traces(ii,:)',20)';
            plot(time_licks,dFA_lick_traces(ii,:)+y_shift,'-m')
            mean_FA_licks(1,:)= mean_FA_licks(1,:) + ((dFA_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end

        %Plot Miss lick traces
        for ii=1:Missii_lick
            dMiss_lick_traces(ii,:)=decimate(Miss_lick_traces(ii,:)',20)';
            plot(time_licks,dMiss_lick_traces(ii,:)+y_shift,'-c')
            mean_Miss_licks(1,:)= mean_Miss_licks(1,:) + ((dMiss_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end

        %PLot Hit lick traces
        for ii=1:Hitii_lick
            dHit_lick_traces(ii,:)=decimate(Hit_lick_traces(ii,:)',20)';
            plot(time_licks,dHit_lick_traces(ii,:)+y_shift,'-r')
            mean_Hit_licks(1,:)= mean_Hit_licks(1,:) + ((dHit_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end

        %Odor on markers
        plot([0 0],[0 y_shift],'-k')
        plot([mean(delta_odor) mean(delta_odor)],[0 y_shift],'-k')

        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 y_shift],'-r')
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 y_shift],'-r')

        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig5.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig5.fig'])
        end

        %Plot lick frequency
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        hold on
        time_licks=([1:length(Hitlick_freq)]*dt_lick)-dt_before;
        plot(time_licks,Hitlick_freq,'-r')
        plot(time_licks,Misslick_freq,'-c')
        plot(time_licks,CRlick_freq,'-b')
        plot(time_licks,FAlick_freq,'-m')

        title('Lick frequency, r=Hit, c=Miss, b=CR, m=FA')
        xlabel('Time (sec)')
        ylabel('Lick frequency (Hz)')

        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig6.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig6.fig'])
        end


        %Get lick p values
        try
            dt_lick_pval=0.1;
            no_pvals=0;
            p_val_Hit_CR=[];
            for ii=1:dt_lick_pval*acq_rate:sum(rhd_mask)-dt_lick_pval*acq_rate
                %Hit vs CR
                this_CR=zeros(CRii_lick,1);
                for jj=1:CRii_lick
                    if sum(CR_lick_traces(jj,ii:ii+dt_lick_pval*acq_rate)>threshold_lick)>=1
                        this_CR(jj,1)=1;
                    else
                        this_CR(jj,1)=0;
                    end
                end

                this_Hit=zeros(Hitii_lick,1);
                for jj=1:Hitii_lick
                    if sum(Hit_lick_traces(jj,ii:ii+dt_lick_pval*acq_rate)>threshold_lick)>=1
                        this_Hit(jj,1)=1;
                    else
                        this_Hit(jj,1)=0;
                    end
                end

                no_pvals=no_pvals+1;

                if (~isempty(this_CR))&(~isempty(this_CR))
                    p_val_Hit_CR(no_pvals)=ranksum(this_CR,this_Hit);
                else
                    p_val_Hit_CR(no_pvals)=1;
                end

                %ranksum gives NaN if the values are all the same
                if isnan(p_val_Hit_CR(no_pvals))
                    p_val_Hit_CR(no_pvals)=1;
                end

            end

            time_p_lick=[-dt_before+(dt_lick_pval/2):dt_lick_pval:dt_after-(dt_lick_pval)];


            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);

            plot(time_p_lick,log10(p_val_Hit_CR))
            hold on
            plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)])
            title('log(p value) for the difference in licks Hit vs CR')
            xlabel('Time (sec)')
            ylabel('log10(p value)')
            if do_warp==1
                savefig([fnameca(1:end-4) '_dropc_warp_Fig7.fig'])
            else
                savefig([fnameca(1:end-4) '_dropc_batch_Fig7.fig'])
            end
        catch
        end
        save_name_catch=[handles_choice.PathNamecsv{fileNo} fnameca(1:end-4) '_try_catch.mat'];
        save(save_name_catch,'try_catch_status')
    catch ME
        try_catch_status(fileNo).processed=0;
        try_catch_status(fileNo).ME=ME;
        save_name_catch=[handles_choice.PathNamecsv{fileNo} fnameca(1:end-4) '_try_catch.mat'];
        save(save_name_catch,'try_catch_status')
     end


     %Save the rhd times

     save_name=[handles_choice.PathNamecsv{fileNo} fnameca(1:end-4) '_rhd_times.mat'];


     %splus_traces(handles_out.no_sp_trials*no_traces,no_timepoints), order is  for trial 1: ROI1, ROI2, etc.. ROIend, trial 2:

     save(save_name, 'odor_on_times_rhd','FV_times_rhd','odor_off_times_rhd')
end



pffft=1;
















