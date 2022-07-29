%% drgcmCaImAn_batch_dropc_rhd_ext.m
%
% Needs as an input the ouput file from EXTRACT
%
close all
clear all

figNo=0;

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{8}='-g';
these_lines{5}='-y';
these_lines{6}='-k';
these_lines{7}='-c';
these_lines{4}='-k';

these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{4}='c';

%Choices
do_warp=0;     %1=warped components from a reference file
um_per_pixel=0.489;

%Other default variables
plot_raw=1; %Used plot_raw=1 for Ming's data
ref_win=[-5 -1.5];


% Read choices file

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dropc_rhd_choices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles_choice.choiceFileName=choiceFileName;
handles_choice.choiceBatchPathName=choiceBatchPathName;

cd(handles_choice.PathName)

handles_per_trial=[];
handles_per_trial.no_trials=0;
handles_per_file=[];
handles_per_file.no_files=0;

if ~isfield(handles_choice, 'No_images')
    handles_choice.No_images=0;
end

if ~isfield(handles_choice, 'No_rhd')
    handles_choice.No_rhd=0;
end

if isfield(handles_choice, 'read_CaImAn_mat')
    read_CaImAn_mat=handles_choice.read_CaImAn_mat;
else
    read_CaImAn_mat=0; %mat or py file for CaImAn
end


%Set dropc_program
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

if ~isfield(handles_choice, 'dropc_program')
    dropc_program=3;
else
    dropc_program=handles_choice.dropc_program;
end


for fileNo=handles_choice.first_file:handles_choice.no_files
    
    fprintf(1, ['\nProcessing  ' num2str(fileNo) '\n']);
    
    if handles_choice.No_images~=1
        
        switch read_CaImAn_mat
            case 1
                Yr=[];
                A_or=[];
                C_or=[];
                b2=[];
                f2=[];
                Cn=[];
                options=[];
                if do_warp==1
                    [Yr,A_or,C_or,b2,f2,Cn,options] = drgCaImAn_get_warped_components(handles_choice.CaImAn_scriptFileName{fileNo},handles_choice.CaImAn_referenceFileName);
                else
                    load(handles_choice.CaImAn_scriptFileName{fileNo})
                end
                
                
                % In order to understand what these variables are take a look at
                % plot_components_GUI.m under CaImAn-MATLAB/utilities
                % CaImAn calls ROIs "components"
                % Yr hs the 6000 images in a vector of 65536 pixels
                % A_or has the image of each component e.g. 65536x158 double for 256x256
                % images of 85 components
                % C_or these are the dFF traces as a function of time for each component
                % b2 are the background images
                % f2 are the background traces
                % Cn is the correlation image e.g. 256x256
                % options has all the options used by the user
                
                fnameca=handles_choice.CaImAn_scriptFileName{fileNo};
                
                %     close all
                
                %Calculate the center of mass (coms) for each component
                %and draw the components and  their locations
                thr = 0.95;
                d1 = options.d1;
                d2 = options.d2;
                dt=1/options.fr;
                
                %Draw the components for the first image
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                cla
                imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
                Cn1x2=2*Cn;
                title(['Components for file No ' num2str(fileNo)])
                
                
                %Find the number of components
                szA_or=size(A_or);
                coms=[];
                num_coms=0;
                
                for i=1:szA_or(2)
                    A_temp = full(reshape(A_or(:,i),d1,d2));
                    A_temp = medfilt2(A_temp,[3,3]);
                    A_temp = A_temp(:);
                    [temp,ind] = sort(A_temp(:).^2,'ascend');
                    temp =  cumsum(temp);
                    ff = find(temp > (1-thr)*temp(end),1,'first');
                    if ~isempty(ff)
                        [contour_properties,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
                        ww.LineWidth = 2;
                        num_coms=num_coms+1;
                        coms(1,num_coms)=mean(contour_properties(1,2:end));
                        coms(2,num_coms)=mean(contour_properties(2,2:end));
                        plot(coms(1,num_coms),coms(2,num_coms),'.r')
                    end
                end
                
                
                raw=[];
                inferred=[];
                try
                    [raw,inferred]=drgGetCAtraces(Yr,A_or,C_or,b2,f2,Cn,options);
                catch
                    pffft=1
                end
                % Should we use raw or inferred?
                
                traces=raw;
                
                traces_inferred=inferred;
            case 2
                %Python input
                dt=handles_choice.CaImAn_dt;
                raw=readmatrix(handles_choice.CaImAn_rawFileName{fileNo});
                inferred=readmatrix(handles_choice.CaImAn_inferredFileName{fileNo});
                coms=readmatrix(handles_choice.CaImAn_comsFileName{fileNo});
                fnameca=handles_choice.CaImAn_rawFileName{fileNo};
                % Should we use raw or inferred?
                
                traces=raw;
                
                traces_inferred=inferred;
            case 3
                 
                %This reads the extract file
                load([handles_choice.PathNamecsv{fileNo} handles_choice.csvFileName{fileNo}])
                traces=zeros(size(output.temporal_weights,1),size(output.temporal_weights,2));
                for traceNo=1:size(output.temporal_weights,2)
                    traces(:,traceNo)=output.temporal_weights(:,traceNo);
                end
                
                
                traces=traces';
                fnameca=handles_choice.csvFileName{fileNo};
                dt=handles_choice.CaImAn_dt;
        end
  
        
        
        sz_traces=size(traces);
        no_traces=sz_traces(1);
        no_images=sz_traces(2);
        
        fprintf(1, ['\ndrgCaImAn_dropc run for ' fnameca '\n\n']);
    end
    
    %Read the dropc file
    handles=[];
    load(handles_choice.spmFileName{fileNo})
    
    %Read the rhd file and register rhd times to dropc
    if handles_choice.No_rhd~=1
        
        %Read the rhd file
        adc_in=[];
        digital_in=[];
        acq_rate=[];
        amplifier_data=[];
        [adc_in,digital_in,acq_rate,amplifier_data]=drg_read_Intan_RHD2000_file(handles_choice.rhdFileName{fileNo},3);
        
        fnamerhd=handles_choice.rhdFileName{fileNo};
        
        
        %Align the rhd times with the olfactometer
        %Note: The olfactometer is triggered by image acquisition
        %This calculation yields delta_t_rhd that is used to offset the rhd traces
        %to the imaging traces
        
        %Get odor on times for the spm
        odor_on_times=[];
        ootii=0;
        for epoch=1:handles.dropcData.epochIndex
            %Epoch 2 is odor on
            if (handles.dropcData.epochEvent(epoch)==2)
                ootii=ootii+1;
                odor_on_times(ootii)=handles.dropcData.epochTime(epoch);
            end
        end
        
        
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
        
        if length(odor_on_times)==length(odor_on_times_rhd)
            delta_t_rhd=mean(odor_on_times-odor_on_times_rhd);
        else
            if length(odor_on_times)>length(odor_on_times_rhd)
                delta_t_rhd=mean(odor_on_times(1:length(odor_on_times_rhd))-odor_on_times_rhd);
            else
                delta_t_rhd=mean(odor_on_times-odor_on_times_rhd(1:length(odor_on_times)));
            end
        end
    end
    
    % Plot the calcium traces and licks
    handles_per_file.no_files=handles_per_file.no_files+1;
    if handles_choice.No_images~=1
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        
        
        hold on
        
        % Determine the y spacing of the traces
        y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));
        
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
                    try
                        plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_traces+2)*y_shift],...
                            these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                    catch
                    end
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
         
        %Plot the licks recorded by the INTAN (adc_in)
        if handles_choice.No_rhd~=1
            time_rhd=([1:length(digital_in)]/acq_rate)+delta_t_rhd;
            pct998=prctile(adc_in,99.8);
            pct1=prctile(adc_in,1);
            norm_fact=0.8*y_shift/(pct998-pct1);
            
            plot(decimate(time_rhd(time_rhd>0),100),decimate(adc_in(time_rhd>0)*norm_fact,100))
        end
        
        %Plot the traces
        time=[1:no_images]*dt;
        for trNo=1:no_traces
            % for trNo=1:20
            plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        end
        
        
        ylim([-y_shift*0.2 (no_traces+2)*y_shift])
        xlim([0 max(time)])
        
        xlabel('time (s)')
        ylabel('deltaF/F')
        title(['Raw dFF traces for file No ' num2str(fileNo)])
        
        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_Fig1.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_Fig1.fig'])
        end
        
%         %Plot inferred traces
%         figNo=figNo+1;
%         try
%             close(figNo)
%         catch
%         end
%         hFig=figure(figNo);
%         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
%         
%         
%         hold on
%         
%         % Determine the y spacing of the traces
%         y_shift=1.2*(prctile(traces_inferred(:),95)-prctile(traces_inferred(:),5));
%         
%         %Plot the event lines
%         switch dropc_program
%             case 1
%                 for event=2:handles.dropcData.eventIndex
%                     plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_traces+2)*y_shift],...
%                         these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
%                 end
%             case 2
%                 for event=1:handles.dropcData.allTrialIndex
%                     try
%                         plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_traces+2)*y_shift],...
%                             these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
%                     catch
%                     end
%                 end
%             case 3
%                 %For S+ and S- plot odor on and reinforcement
%                 for epoch=1:handles.dropcData.epochIndex
%                     %Epoch 2 is odor on, 3 is odor off
%                     plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
%                     if plot_epoch
%                         if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
%                             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
%                                 '-r','LineWidth',1)
%                         else
%                             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
%                                 '-b','LineWidth',1)
%                         end
%                     end
%                 end
%         end
%         
%         %Plot the licks recorded by the INTAN (adc_in)
%         if handles_choice.No_rhd~=1
%             plot(decimate(time_rhd(time_rhd>0),100),decimate(adc_in(time_rhd>0)*norm_fact,100))
%         end
%         
%         %Plot the inferred traces
%         for trNo=1:no_traces
%             % for trNo=1:20
%             plot(time,traces_inferred(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
%         end
%         
%         
%         ylim([-y_shift*0.2 (no_traces+2)*y_shift])
%         xlim([0 max(time)])
%         
%         xlabel('time (s)')
%         ylabel('deltaF/F')
%         title(['Inferred dFF traces for file No ' num2str(fileNo)])
%         
%         if do_warp==1
%             savefig([fnameca(1:end-4) '_infdropc_warp_Fig1.fig'])
%         else
%             savefig([fnameca(1:end-4) '_infdropc_batch_Fig1.fig'])
%         end
         
        
        handles_per_file.file(handles_per_file.no_files).dFFtraces=traces;
%         handles_per_file.file(handles_per_file.no_files).dFFtraces_inferred=traces_inferred;
        handles_per_file.file(handles_per_file.no_files).dFF_time=time;
        handles_per_file.file(handles_per_file.no_files).no_dFF_components=no_traces;
    end
    
    if length(decimate(time_rhd,100))~=length(decimate(amplifier_data,100))
        handles_choice.No_LFP=1;
    else
        handles_choice.No_LFP=0;
    end
    
    
    %Plot the LFPs
    if handles_choice.No_LFP==0
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        
        hold on
        
        y_shift=1.2*(prctile(amplifier_data(:),99)-prctile(amplifier_data(:),1));
        
        %Plot the event lines
        no_electrodes=size(amplifier_data,1);
        
        switch dropc_program
            case 1
                for event=2:handles.dropcData.eventIndex
                    plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_electrodes+2)*y_shift],...
                        these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
                end
            case 2
                for event=1:handles.dropcData.allTrialIndex
                    plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_electrodes+2)*y_shift],...
                        these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                end
            case 3
                %For S+ and S- plot odor on and reinforcement
                for epoch=1:handles.dropcData.epochIndex
                    %Epoch 2 is odor on, 3 is odor off
                    plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
                    if plot_epoch
                        if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_electrodes+2)*y_shift],...
                                '-r','LineWidth',1)
                        else
                            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_electrodes+2)*y_shift],...
                                '-b','LineWidth',1)
                        end
                        
                    end
                    
                    
                end
        end
        
        %Plot the licks recorded by the INTAN (adc_in)
        time_rhd=([1:length(digital_in)]/acq_rate)+delta_t_rhd;
        pct998=prctile(adc_in,99.8);
        pct1=prctile(adc_in,1);
        norm_fact=0.8*y_shift/(pct998-pct1);
        
        plot(decimate(time_rhd(time_rhd>0),100),decimate(adc_in(time_rhd>0)*norm_fact,100))
        
        
        %Plot the LFP traces
        for elecNo=1:no_electrodes
            % for trNo=1:20
            plot(decimate(time_rhd(time_rhd>0),100),decimate(amplifier_data(elecNo,time_rhd>0)+y_shift*elecNo,100),'-k','LineWidth',1)
        end
         
        handles_per_file.file(handles_per_file.no_files).LFPtraces=amplifier_data;
        handles_per_file.file(handles_per_file.no_files).LFP_time=time_rhd;
        handles_per_file.file(handles_per_file.no_files).no_electrodes=no_electrodes;
        handles_per_file.file(handles_per_file.no_files).acq_rate=acq_rate;
        
        ylim([-y_shift*0.2 (no_electrodes+2)*y_shift])
        xlim([0 max(time_rhd)])
        
        xlabel('time (s)')
        ylabel('LFP')
        title(['LFP traces for file No ' num2str(fileNo)])
        
        if do_warp==1
            savefig([fnamerhd(1:end-4) '_dropc_warp_LFP.fig'])
        else
            savefig([fnamerhd(1:end-4) '_dropc_batch_LFP.fig'])
        end
        
        %Do the wavelet analysis
        dec_n=handles_choice.dec_n; %100
        dectim=decimate(time_rhd(time_rhd>0),dec_n);
        
        %Setup the wavelet scales
        %   scales = helperCWTTimeFreqVector(minfreq,maxfreq,f0,dt,NumVoices)
        %   f0 - center frequency of the wavelet in cycles/unit time
        %   dt - sampling interval
        %   NumVoices - number of voices per octave
         
        NumVoices=5;
        minfreq=handles_choice.minfreq; %1
        maxfreq=handles_choice.maxfreq; %100
%         dt_rhd=1/acq_rate; %Note: this is what we used before
        dt_rhd=1/(acq_rate/dec_n);
        f0=5/(2*pi);
        
        a0 = 2^(1/NumVoices);
        minscale = f0/(maxfreq*dt_rhd);
        maxscale = f0/(minfreq*dt_rhd);
        minscale = floor(NumVoices*log2(minscale));
        maxscale = ceil(NumVoices*log2(maxscale));
        scales = a0.^(minscale:maxscale).*dt_rhd;
        
        %Now do the wavelet transform
        for electNo=1:no_electrodes
            notch60HzFilt = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                'DesignMethod','butter','SampleRate',acq_rate);
            this_LFP=zeros(1,size(amplifier_data(electNo,time_rhd>0),2));
            this_LFP(1,:)=amplifier_data(electNo,time_rhd>0);
            LFP=filtfilt(notch60HzFilt,this_LFP);
            decLFP=decimate(LFP,dec_n);
            decFs=acq_rate/dec_n;
            
            cwtLFP = cwtft({detrend(double(decLFP)),1/decFs},'wavelet','morl','scales',scales);
            Prev=abs(cwtLFP.cfs).^2;
            P=Prev(end:-1:1,:);
            
            frev=cwtLFP.frequencies;
            f=frev(end:-1:1);
            
            if electNo==1
                all_Power_timecourse=zeros(no_electrodes,length(f),length(dectim));
            end
            all_Power_timecourse(electNo,:,:)=P;
            
        end
        
        %Per electrode power plot
        log_P_per_trial_timecourse_sub=zeros(length(f)*no_electrodes,length(dectim));
        log_P_timecourses=zeros(no_electrodes,length(f),length(dectim));
        log_P_timecourses_per_bw=zeros(no_electrodes,length(handles_choice.lowF),length(dectim));
        mean_log_P_timecourse=zeros(length(f),length(dectim));
        y_shift=0;
        sy_shift=0;
        shifted_freq=[];
        for electNo=1:no_electrodes
            this_log_P_timecourse=zeros(length(f),length(dectim));
            this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(electNo,:,:));
            mean_log_P_timecourse(:,:)=mean_log_P_timecourse(:,:)+this_log_P_timecourse(:,:);
            log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(f),:)=this_log_P_timecourse;
            log_P_timecourses(electNo,:,:)=this_log_P_timecourse;
            for bwii=1:length(handles_choice.lowF)
                log_P_timecourses_per_bw(electNo,bwii,:)=mean(this_log_P_timecourse((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),:),1);
            end
            shifted_freq(1,y_shift+1:y_shift+length(f))=f+(electNo-1)*f(end);
            y_shift=y_shift+length(f);
        end
        mean_log_P_timecourse=mean_log_P_timecourse/no_electrodes;
        
        %     handles_per_file.file(handles_per_file.no_files).log_P_timecourses=log_P_timecourses;
        handles_per_file.file(handles_per_file.no_files).log_P_timecourses_per_bw=log_P_timecourses_per_bw;
        %     handles_per_file.file(handles_per_file.no_files).f=f;
        handles_per_file.file(handles_per_file.no_files).log_P_time=dectim;
        
        %Plot the per-electrode timecourse
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        maxLogPper=prctile(log_P_per_trial_timecourse_sub(:),99);
        minLogPper=prctile(log_P_per_trial_timecourse_sub(:),1);
        %Note: Diego added this on purpose to limit the range to 10 dB
        %This results in emphasizing changes in the top 10 dB
        if maxLogPper-minLogPper>12
            minLogPper=maxLogPper-12;
        end
        
        set(hFig, 'units','normalized','position',[.07 .1 .75 .3])
        drg_pcolor(repmat(dectim,length(f)*no_electrodes,1)',repmat(shifted_freq,length(dectim),1),log_P_per_trial_timecourse_sub')
        
        colormap jet
        shading interp
        %     caxis([minLogPper maxLogPper]);
        xlabel('Time (sec)')
        ylabel('Frequency*trialNo');
        title(['Power (dB, wavelet) timecourse per electrode for file No ' num2str(fileNo)])
        
        %Plot the event lines
        hold on
        
        switch dropc_program
            case 1
                for event=2:handles.dropcData.eventIndex
                    plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_electrodes+2)*y_shift],...
                        these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
                end
            case 2
                for event=1:handles.dropcData.allTrialIndex
                    plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_electrodes+2)*y_shift],...
                        these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                end
            case 3
                %For S+ and S- plot odor on and reinforcement
                for epoch=1:handles.dropcData.epochIndex
                    %Epoch 2 is odor on, 3 is odor off
                    plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
                    if plot_epoch
                        if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 max(shifted_freq)],...
                                '-m','LineWidth',2)
                        else
                            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 max(shifted_freq)],...
                                '-b','LineWidth',2)
                        end
                        
                    end
                    
                    
                end
        end
        
        %Plot the mean wavelet power timecourse
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        
        
        set(hFig, 'units','normalized','position',[.07 .1 .75 .3])
        drg_pcolor(repmat(dectim,length(f),1)',repmat(f,length(dectim),1),mean_log_P_timecourse')
        
        colormap jet
        shading interp
        %     caxis([minLogPper maxLogPper]);
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)');
        title(['Mean power (dB, wavelet) timecourse  for file No ' num2str(fileNo)])
        
        %Plot the event lines
        hold on
        
        switch dropc_program
            case 1
                for event=2:handles.dropcData.eventIndex
                    plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_electrodes+2)*y_shift],...
                        these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
                end
            case 2
                for event=1:handles.dropcData.allTrialIndex
                    plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_electrodes+2)*y_shift],...
                        these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                end
            case 3
                %For S+ and S- plot odor on and reinforcement
                for epoch=1:handles.dropcData.epochIndex
                    %Epoch 2 is odor on, 3 is odor off
                    plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
                    if plot_epoch
                        if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 max(shifted_freq)],...
                                '-m','LineWidth',2)
                        else
                            plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 max(shifted_freq)],...
                                '-b','LineWidth',2)
                        end
                        
                    end
                    
                    
                end
        end
        
        %Compute the mean LFP for the different bandwidths
        
        %I will subtract the one percentile
        szmlp=size(mean_log_P_timecourse);
        
        
        for bwii=1:length(handles_choice.lowF)
            ref_logP=prctile(mean(mean_log_P_timecourse((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),:),1),1)*ones(1,szmlp(2));
            handles_out2.logPtrace(bwii,1:szmlp(2))=mean(mean_log_P_timecourse((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),:),1)-ref_logP;
            max_logP(bwii)=prctile(mean(mean_log_P_timecourse((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),:),1)-ref_logP,99);
            min_logP(bwii)=prctile(mean(mean_log_P_timecourse((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),:),1)-ref_logP,1);
        end
        
        
        
        
        for bwii=1:length(handles_choice.lowF)
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            hFig=figure(figNo);
            
            
            set(hFig, 'units','normalized','position',[.07 .1 .75 .3])
            hold on
            plot(dectim,handles_out2.logPtrace(bwii,1:length(dectim)),'-k')
            
            switch dropc_program
                case 1
                    for event=2:handles.dropcData.eventIndex
                        plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_electrodes+2)*y_shift],...
                            these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
                    end
                case 2
                    for event=1:handles.dropcData.allTrialIndex
                        plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_electrodes+2)*y_shift],...
                            these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                    end
                case 3
                    %For S+ and S- plot odor on and reinforcement
                    for epoch=1:handles.dropcData.epochIndex
                        %Epoch 2 is odor on, 3 is odor off
                        plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
                        if plot_epoch
                            if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                                plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [min_logP(bwii)-0.2*(max_logP(bwii)-min_logP(bwii)) max_logP(bwii)+0.2*(max_logP(bwii)-min_logP(bwii))],...
                                    '-m','LineWidth',2)
                            else
                                plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [min_logP(bwii)-0.2*(max_logP(bwii)-min_logP(bwii)) max_logP(bwii)+0.2*(max_logP(bwii)-min_logP(bwii))],...
                                    '-b','LineWidth',2)
                            end
                            
                        end
                        
                        
                    end
            end
            
            xlabel('Time (sec)')
            ylabel('dB')
            xlim([min(dectim) max(dectim)])
            ylim([min_logP(bwii)-0.2*(max_logP(bwii)-min_logP(bwii)) max_logP(bwii)+0.2*(max_logP(bwii)-min_logP(bwii))])
            title(['delta power ' handles_choice.bw_names(bwii)])
        end
         
        
    end
    
    %Plot the event lines
    
    
    %Start extracting trials
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
                if ((handles.dropcData.epochTime(epoch)-dt_before)>0)
                    %&(handles.dropcData.epochTime(epoch)+dt_after<=max(time))
                    
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
                if ((handles.dropcData.epochTime(epoch)-dt_before)>0)
                    %&(handles.dropcData.epochTime(epoch)+dt_after<=max(time))
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
            fvelii=0;
            fv_deltalogP=[];
            
            
            sp_odor_response=[];
            splus_traces=[];
            spii=0;
            spii_lick=0;
            spelii=0;
            sp_deltalogP=[];
            smii=0;
            smii_lick=0;
            smelii=0;
            sm_deltalogP=[];
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
            Hitelii=0;
            Hit_deltalogP=[];
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
            Misselii=0;
            Miss_deltalogP=[];
            which_trace_Miss=[];
            which_trial_Miss=[];
            which_Missii_lick=[];
            Miss_lick_times=[];
            Miss_no_lick_times=0;
            
            FAii=0;
            FAii_lick=0;
            no_FA_traces=0;
            FA_traces=[];
            FAelii=0;
            FA_deltalogP=[];
            which_trace_FA=[];
            which_trial_FA=[];
            which_FAii_lick=[];
            FA_lick_times=[];
            FA_no_lick_times=0;
            
            
            CRii=0;
            CRii_lick=0;
            no_CR_traces=0;
            CR_traces=[];
            CRelii=0;
            CR_deltalogP=[];
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
                if ((handles.dropcData.epochTime(epoch)-dt_before)>0)
                    %&(handles.dropcData.epochTime(epoch)+dt_after<=max(time))
                    %This event is processed for Ca
                    
                    %Final valve epoch
                    if (handles.dropcData.epochEvent(epoch)==1)
                        
                        if handles_choice.No_images~=1
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
                        end
                        
                        %Get the licks
                        rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before)...
                            &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after);
                        fvii_lick=fvii_lick+1;
                        fv_lick_traces(fvii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);
                        
                        %Get the wavelet LFPs for each bandwidth
                        lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                            &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                        lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                            &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                        handles_out.no_electrodes=no_electrodes;
                        
                        for bwii=1:length(handles_choice.lowF)
                            for electNo=1:no_electrodes
                                fvelii=fvelii+1;
                                this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                this_meanlogP(1,:)=mean(these_logP,1);
                                
                                these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                
                                this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                fv_deltalogP(bwii,fvelii,1:sum(lfp_snip_mask))=this_deltalogP;
                                handles_out.electrodeNo(electNo).bwii(bwii).trialNo(fvelii).fv_logP=fv_deltalogP(bwii,fvelii,1:sum(lfp_snip_mask));
                            end
                        end
                        
                    end
                    handles_out.no_fv_trials=fvii_lick;
                    
                    %Now do S+ and S-
                    if (handles.dropcData.epochEvent(epoch)==2)
                        
                        %Odor on calcium
                        if handles_choice.No_images~=1
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
                        end
                        
                        if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                            
                            %S+
                            
                            if (handles_choice.No_images~=1)&(sum(snip_mask)>0)
                                %S plus calcium
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
                            end
                            
                            %Get the wavelet LFPs for each bandwidth
                            lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                                &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                            lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            
                            for bwii=1:length(handles_choice.lowF)
                                for electNo=1:no_electrodes
                                    spelii=spelii+1;
                                    this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                    these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                    these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                    this_meanlogP(1,:)=mean(these_logP,1);
                                    
                                    these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                    these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                    this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                    
                                    this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                    this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                    sp_deltalogP(bwii,spelii,1:sum(lfp_snip_mask))=this_deltalogP;
                                    handles_out.electrodeNo(electNo).bwii(bwii).trialNo(spelii).sp_logP=sp_deltalogP(bwii,spelii,1:sum(lfp_snip_mask));
                                end
                            end
                            
                            
                        else
                            
                            %S minus
                            
                            if (handles_choice.No_images~=1)&(sum(snip_mask)>0)
                                %Calcium
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
                            end
                            
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            smii_lick=smii_lick+1;
                            sminus_lick_traces(smii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);
                            
                            handles_out.no_sm_trials=smii_lick;
                            
                            %Get the wavelet LFPs for each bandwidth
                            lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                                &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                            lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                                &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                            
                            for bwii=1:length(handles_choice.lowF)
                                for electNo=1:no_electrodes
                                    smelii=smelii+1;
                                    this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                    these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                    these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                    this_meanlogP(1,:)=mean(these_logP,1);
                                    
                                    these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                    these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                    this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                    
                                    this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                    this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                    sm_deltalogP(bwii,smelii,1:sum(lfp_snip_mask))=this_deltalogP;
                                    handles_out.electrodeNo(electNo).bwii(bwii).trialNo(smelii).sm_logP=sm_deltalogP(bwii,smelii,1:sum(lfp_snip_mask));
                                end
                            end
                        end
                        %end
                    end
                    
                    %Find Hit, CR, FA and Miss and record data for handles_per_trial
                    
                    %Hit
                    if (handles.dropcData.epochEvent(epoch)==6)
                        
                        handles_per_trial.no_trials=handles_per_trial.no_trials+1;
                        handles_per_trial.trial(handles_per_trial.no_trials).hit=1;
                        handles_per_trial.trial(handles_per_trial.no_trials).miss=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).CR=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).FA=0;
                        
                        no_odor_trials=no_odor_trials+1;
                        no_Hit_trials=no_Hit_trials+1;
                        valid_trace(no_odor_trials)=1;
                        epoch_per_trial(no_odor_trials)=6;
                        epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                        lda_event{no_odor_trials}='S+';
                        
                        handles_out.epoch_per_trial(no_odor_trials)=6;
                        
                        %Do Ca traces
                        if handles_choice.No_images~=1
                            
                            
                            
                            
                            
                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
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
                                    
                                    
                                    
                                    no_Hit_traces=no_Hit_traces+1;
                                    Hit_trials(Hitii)=no_Hit_trials;
                                    which_trace_Hit(Hitii)=trNo;
                                    which_trial_Hit(Hitii)=no_odor_trials;
                                    which_Hitii_lick(Hitii)=Hitii_lick+1;
                                    
                                    handles_out.componentNo(trace_num).trialNo(Hitii_lick+1).hit_traces=Hit_traces(Hitii,1:sum(snip_mask));
                                    handles_out.trialNo(Hitii_lick+1).trace_numHit=trace_num;
%                                     if read_CaImAn_mat==1
%                                         handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(:,trNo);
%                                     else
%                                         handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(trNo,1:2)';
%                                     end
                                    
                                    
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                    
                                    handles_per_trial.trial(handles_per_trial.no_trials).no_traces=trace_num;
                                    
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).fileNo=fileNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).trNo=trNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).traces=Hit_traces(Hitii,1:sum(snip_mask));
                                    
                                end
                            end
                        end
                        
                        %Now do licks
                        rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                            &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                        Hitii_lick=Hitii_lick+1;
                        Hit_lick_traces(Hitii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        
                        
                        these_lick_times=[];
                        these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                        Hit_lick_times(Hitii_lick,1:length(these_lick_times))=these_lick_times;
                        Hit_no_lick_times(Hitii_lick)=length(these_lick_times);
                        
                        handles_per_trial.trial(handles_per_trial.no_trials).lick_trace=adc_in(rhd_mask);
                        handles_per_trial.trial(handles_per_trial.no_trials).threshold_lick=threshold_lick;
                        handles_per_trial.trial(handles_per_trial.no_trials).these_lick_times=these_lick_times;
                        
                        
                        handles_out.no_Hit_trials=Hitii_lick;
                        handles_out.Hit_trial_no(no_odor_trials)=Hitii_lick;
                        
                        %Get the wavelet LFPs for each bandwidth
                        lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                            &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                        lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                            &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                        
                        for bwii=1:length(handles_choice.lowF)
                            for electNo=1:no_electrodes
                                Hitelii=Hitelii+1;
                                this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                this_meanlogP(1,:)=mean(these_logP,1);
                                
                                these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                
                                this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                Hit_deltalogP(bwii,Hitelii,1:sum(lfp_snip_mask))=this_deltalogP;
                                handles_out.electrodeNo(electNo).bwii(bwii).trialNo(Hitelii).Hit_logP=Hit_deltalogP(bwii,Hitelii,1:sum(lfp_snip_mask));
                                handles_per_trial.trial(handles_per_trial.no_trials).bwii(bwii).electrode(electNo).logP=Hit_deltalogP(bwii,Hitelii,1:sum(lfp_snip_mask));
                            end
                        end
                    end
                    
                    %Miss
                    if (handles.dropcData.epochEvent(epoch)==7)
                        
                        
                        handles_per_trial.no_trials=handles_per_trial.no_trials+1;
                        handles_per_trial.trial(handles_per_trial.no_trials).hit=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).miss=1;
                        handles_per_trial.trial(handles_per_trial.no_trials).CR=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).FA=0;
                        
                        no_odor_trials=no_odor_trials+1;
                        
                        valid_trace(no_odor_trials)=1;
                        epoch_per_trial(no_odor_trials)=7;
                        epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                        lda_event{no_odor_trials}='S+';
                        
                        handles_out.epoch_per_trial(no_odor_trials)=7;
                        
                        if handles_choice.No_images~=1
                            
                            
                            
                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
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
                                    handles_out.trialNo(Missii_lick+1).trace_numMiss=trace_num;
                                    no_Miss_traces=no_Miss_traces+1;
                                    which_trace_Miss(Missii)=trNo;
                                    which_trial_Miss(Missii)=no_odor_trials;
                                    which_Missii_lick(Missii)=Missii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                    
                                    handles_per_trial.trial(handles_per_trial.no_trials).no_traces=trace_num;
%                                     if read_CaImAn_mat==1
%                                         handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(:,trNo);
%                                     else
%                                         handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(trNo,1:2)';
%                                     end
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).fileNo=fileNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).trNo=trNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).traces=Miss_traces(Missii,1:sum(snip_mask));
                                    
                                end
                            end
                            
                        end
                        
                        %Now do licks
                        rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                            &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                        Missii_lick=Missii_lick+1;
                        Miss_lick_traces(Missii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        
                        these_lick_times=[];
                        these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                        Miss_lick_times(Missii_lick,1:length(these_lick_times))=these_lick_times;
                        Miss_no_lick_times(Missii_lick)=length(these_lick_times);
                        
                        handles_per_trial.trial(handles_per_trial.no_trials).lick_trace=adc_in(rhd_mask);
                        handles_per_trial.trial(handles_per_trial.no_trials).threshold_lick=threshold_lick;
                        handles_per_trial.trial(handles_per_trial.no_trials).these_lick_times=these_lick_times;
                        
                        handles_out.no_Miss_trials=Missii_lick;
                        handles_out.Miss_trial_no(no_odor_trials)=Missii_lick;
                        
                        %Get the wavelet LFPs for each bandwidth
                        lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                            &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                        lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                            &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                        
                        for bwii=1:length(handles_choice.lowF)
                            for electNo=1:no_electrodes
                                Misselii=Misselii+1;
                                this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                this_meanlogP(1,:)=mean(these_logP,1);
                                
                                these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                
                                this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                Miss_deltalogP(bwii,Misselii,1:sum(lfp_snip_mask))=this_deltalogP;
                                handles_out.electrodeNo(electNo).bwii(bwii).trialNo(Misselii).Miss_logP=Miss_deltalogP(bwii,Misselii,1:sum(lfp_snip_mask));
                                handles_per_trial.trial(handles_per_trial.no_trials).bwii(bwii).electrode(electNo).logP=Miss_deltalogP(bwii,Misselii,1:sum(lfp_snip_mask));
                            end
                        end
                    end
                    
                    %FA
                    if (handles.dropcData.epochEvent(epoch)==8)
                        
                        
                        
                        handles_per_trial.no_trials=handles_per_trial.no_trials+1;
                        handles_per_trial.trial(handles_per_trial.no_trials).hit=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).miss=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).CR=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).FA=1;
                        
                        no_odor_trials=no_odor_trials+1;
                        valid_trace(no_odor_trials)=1;
                        epoch_per_trial(no_odor_trials)=8;
                        epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                        lda_event{no_odor_trials}='S-';
                        
                        handles_out.epoch_per_trial(no_odor_trials)=8;
                        
                        if handles_choice.No_images~=1
                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
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
                                    handles_out.trialNo(FAii_lick+1).trace_numFA=trace_num;
                                    no_FA_traces=no_FA_traces+1;
                                    which_trace_FA(FAii)=trNo;
                                    which_trial_FA(FAii)=no_odor_trials;
                                    which_FAii_lick(FAii)=FAii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                    
                                    handles_per_trial.trial(handles_per_trial.no_trials).no_traces=trace_num;
%                                     handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(:,trNo);
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).fileNo=fileNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).trNo=trNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).traces=FA_traces(FAii,1:sum(snip_mask));
                                    
                                end
                            end
                        end
                        
                        %Now do licks
                        rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                            &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                        FAii_lick=FAii_lick+1;
                        FA_lick_traces(FAii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        
                        these_lick_times=[];
                        these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                        FA_lick_times(FAii_lick,1:length(these_lick_times))=these_lick_times;
                        FA_no_lick_times(FAii_lick)=length(these_lick_times);
                        
                        handles_per_trial.trial(handles_per_trial.no_trials).lick_trace=adc_in(rhd_mask);
                        handles_per_trial.trial(handles_per_trial.no_trials).threshold_lick=threshold_lick;
                        handles_per_trial.trial(handles_per_trial.no_trials).these_lick_times=these_lick_times;
                        
                        handles_out.no_FA_trials=FAii_lick;
                        handles_out.FA_trial_no(no_odor_trials)=FAii_lick;
                        
                        %Get the wavelet LFPs for each bandwidth
                        lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                            &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                        lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                            &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                        
                        for bwii=1:length(handles_choice.lowF)
                            for electNo=1:no_electrodes
                                FAelii=FAelii+1;
                                this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                this_meanlogP(1,:)=mean(these_logP,1);
                                
                                these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                
                                this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                FA_deltalogP(bwii,FAelii,1:sum(lfp_snip_mask))=this_deltalogP;
                                handles_out.electrodeNo(electNo).bwii(bwii).trialNo(FAelii).FA_logP=FA_deltalogP(bwii,FAelii,1:sum(lfp_snip_mask));
                                handles_per_trial.trial(handles_per_trial.no_trials).bwii(bwii).electrode(electNo).logP=FA_deltalogP(bwii,FAelii,1:sum(lfp_snip_mask));
                            end
                        end
                    end
                    
                    %CR
                    if (handles.dropcData.epochEvent(epoch)==9)
                        
                        
                        
                        handles_per_trial.no_trials=handles_per_trial.no_trials+1;
                        handles_per_trial.trial(handles_per_trial.no_trials).hit=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).miss=0;
                        handles_per_trial.trial(handles_per_trial.no_trials).CR=1;
                        handles_per_trial.trial(handles_per_trial.no_trials).FA=0;
                        
                        
                        no_odor_trials=no_odor_trials+1;
                        valid_trace(no_odor_trials)=1;
                        epoch_per_trial(no_odor_trials)=9;
                        epoch_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                        lda_event{no_odor_trials}='S-';
                        
                        handles_out.epoch_per_trial(no_odor_trials)=9;
                        
                        if handles_choice.No_images~=1
                            
                            snip_mask=(time>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
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
                                    handles_out.trialNo(CRii_lick+1).trace_numCR=trace_num;
                                    no_CR_traces=no_CR_traces+1;
                                    which_trace_CR(CRii)=trNo;
                                    which_trial_CR(CRii)=no_odor_trials;
                                    which_CRii_lick(CRii)=CRii_lick+1;
                                    lda_input_timecourse(1:sum(snip_mask),trace_num,no_odor_trials)=traces(trNo, snip_mask)-mean(traces(trNo, ref_mask));
                                    
                                    handles_per_trial.trial(handles_per_trial.no_trials).no_traces=trace_num;
%                                     if read_CaImAn_mat==1
%                                         handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(:,trNo);
%                                     else
%                                         handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).roi_location(1:2)=coms(trNo,1:2)';
%                                     end
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).fileNo=fileNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).trNo=trNo;
                                    handles_per_trial.trial(handles_per_trial.no_trials).trace(trace_num).traces=CR_traces(CRii,1:sum(snip_mask));
                                    
                                end
                            end
                        end
                        
                        %Get the licks
                        
                        rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                            &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                        CRii_lick=CRii_lick+1;
                        CR_lick_traces(CRii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        
                        these_lick_times=[];
                        these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                        CR_lick_times(CRii_lick,1:length(these_lick_times))=these_lick_times;
                        CR_no_lick_times(CRii_lick)=length(these_lick_times);
                        
                        
                        handles_per_trial.trial(handles_per_trial.no_trials).lick_trace=adc_in(rhd_mask);
                        handles_per_trial.trial(handles_per_trial.no_trials).threshold_lick=threshold_lick;
                        handles_per_trial.trial(handles_per_trial.no_trials).these_lick_times=these_lick_times;
                        
                        handles_out.no_CR_trials=CRii_lick;
                        handles_out.CR_trial_no(no_odor_trials)=CRii_lick;
                        
                        %Get the wavelet LFPs for each bandwidth
                        lfp_snip_mask=(dectim>=handles.dropcData.epochTime(epoch)-dt_before)...
                            &(dectim<=handles.dropcData.epochTime(epoch)+dt_after);
                        lfp_ref_mask=(dectim>=handles.dropcData.epochTime(epoch)+ref_win(1))...
                            &(dectim<=handles.dropcData.epochTime(epoch)+ref_win(2));
                        
                        for bwii=1:length(handles_choice.lowF)
                            for electNo=1:no_electrodes
                                CRelii=CRelii+1;
                                this_meanlogP=zeros(1,sum(lfp_snip_mask));
                                these_logP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_snip_mask));
                                these_logP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_snip_mask));
                                this_meanlogP(1,:)=mean(these_logP,1);
                                
                                these_reflogP=zeros(sum((f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii))),sum(lfp_ref_mask));
                                these_reflogP(:,:)=10*log10(all_Power_timecourse(electNo,(f>=handles_choice.lowF(bwii))&(f<handles_choice.highF(bwii)),lfp_ref_mask));
                                this_refLFP=mean(mean(these_reflogP,1))*ones(1,sum(lfp_snip_mask));
                                
                                this_deltalogP=zeros(1,sum(lfp_snip_mask));
                                this_deltalogP(:,:)=this_meanlogP-this_refLFP;
                                CR_deltalogP(bwii,CRelii,1:sum(lfp_snip_mask))=this_deltalogP;
                                handles_out.electrodeNo(electNo).bwii(bwii).trialNo(CRelii).CR_logP=CR_deltalogP(bwii,CRelii,1:sum(lfp_snip_mask));
                                handles_per_trial.trial(handles_per_trial.no_trials).bwii(bwii).electrode(electNo).logP=CR_deltalogP(bwii,CRelii,1:sum(lfp_snip_mask));
                            end
                        end
                    end
                    
                end
            end
    end
    
    if handles_choice.No_images~=1
        %Plot the snips
        szFV=size(fv_traces);
        time_to_event=([1:szFV(2)]*dt-dt_before);
        szFVl=size(fv_lick_traces);
        time_to_event_l=([1:szFVl(2)/20]*(20/acq_rate)-dt_before);
        handles_out.time_to_eventFV=time_to_event;
    end
    
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
    
    if handles_choice.No_images~=1
        %FV dFF
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
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
    end
    
    %FV delta dB power
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    hold on
    maxfv=-20000;
    minfv=+20000;
    lfp_t_event=[-dt_before:(dt_after+dt_before)/sum(lfp_snip_mask):dt_after];
    if length(lfp_t_event)>sum(lfp_snip_mask)
        lfp_t_event=lfp_t_event(1:sum(lfp_snip_mask));
    end
     
    for bwii=1:length(handles_choice.lowF)
        subplot(length(handles_choice.lowF),1,bwii)
        hold on
        these_fv_deltalogPs=zeros(fvelii,sum(lfp_snip_mask));
        these_fv_deltalogPs(:,:)=fv_deltalogP(bwii,:,1:sum(lfp_snip_mask));
        CI = bootci(1000, @mean, these_fv_deltalogPs);
        meanfv=mean(these_fv_deltalogPs,1);
        CI(1,:)=meanfv-CI(1,:);
        CI(2,:)=CI(2,:)-meanfv;
        [hl1, hp1] = boundedline(lfp_t_event',meanfv', CI', these_colors{bwii});
        maxfv=max(meanfv);
        minfv=min(meanfv);
        plot([0 0],[minfv-0.2*(maxfv-minfv) maxfv+0.2*(maxfv-minfv)],'-k')
        ylim([minfv-0.2*(maxfv-minfv) maxfv+0.2*(maxfv-minfv)])
        if bwii==length(handles_choice.lowF)
            xlabel('Time (sec)')
        end
        ylabel('dB')
        title(handles_choice.bw_names(bwii))
    end
    
    sgtitle('delta power aligned to final valve diversion')
    
    if do_warp==1
        savefig([fnamerhd(1:end-4) '_dropc_warp_fvdeltadV.fig'])
    else
        savefig([fnamerhd(1:end-4) '_dropc_batch_fvdeltadV.fig'])
    end
    
    if handles_choice.No_images~=1
        %Plot spm for Ca
        
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
        hFig=figure(figNo);
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
        
        
        
        title("Ca changes aligned to odor onset, S+ and S-")
        legend([hlsp hlsm odorhl reinfhl],'S+','S-','Odor','Reinforcement')
        xlabel('Time (sec)')
        ylabel('dF/F')
        ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
        xlim([-10 19.8])
        
        if do_warp==1
            savefig([fnameca(1:end-4) '_dropc_warp_spmCa.fig'])
        else
            savefig([fnameca(1:end-4) '_dropc_batch_spmCa.fig'])
        end
    end
    
    %Plot spm for delta LFP
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    
    lfp_t_event=[-dt_before:(dt_after+dt_before)/sum(lfp_snip_mask):dt_after];
    if length(lfp_t_event)>sum(lfp_snip_mask)
        lfp_t_event=lfp_t_event(1:sum(lfp_snip_mask));
    end
    
    no_time_points_sm=min([sum(lfp_snip_mask) size(sm_deltalogP,3)]);
    no_time_points_sp=min([sum(lfp_snip_mask) size(sp_deltalogP,3)]);
    for bwii=1:length(handles_choice.lowF)
        subplot(length(handles_choice.lowF),1,bwii)
        hold on
        
        these_sm_deltalogPs=zeros(smelii,no_time_points_sm);
        these_sm_deltalogPs(:,:)=sm_deltalogP(bwii,:,1:no_time_points_sm);
        CIsm = bootci(1000, @mean, these_sm_deltalogPs);
        meansm=mean(these_sm_deltalogPs,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        [hlsm, hpsm] = boundedline(lfp_t_event(1:no_time_points_sm)',meansm', CIsm', 'b');
         
        these_sp_deltalogPs=zeros(spelii,no_time_points_sp);
        these_sp_deltalogPs(:,:)=sp_deltalogP(bwii,:,1:no_time_points_sp);
        CIsp = bootci(1000, @mean, these_sp_deltalogPs);
        meansp=mean(these_sp_deltalogPs,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
        [hlsp, hpsp] = boundedline(lfp_t_event(1:no_time_points_sp)',meansp', CIsp', 'r');
        
        pct99=max(meanfv);
        pct1=min(meanfv);
        
        %Odor on markers
        plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
        odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        
        
        ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
        if bwii==length(handles_choice.lowF)
            xlabel('Time (sec)')
            legend([hlsp hlsm],'S+','S-')
        end
        ylabel('dB')
        xlim([-10 19.8])
        title(handles_choice.bw_names(bwii))
    end
    
    sgtitle('delta power for S+ and S- aligned to odor on')
    
    
    
    if do_warp==1
        savefig([fnamerhd(1:end-4) '_dropc_warp_spmdeltaLFP.fig'])
    else
        savefig([fnamerhd(1:end-4) '_dropc_batch_spmdeltaLFP.fig'])
    end
    
    if handles_choice.No_images~=1
        %Plot Hit, CR, et al for Ca traces
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
        
        
        %Save the calculated data
        if do_warp==1
            save_name=[fnameca(1:end-4) '_warp_pre_per.mat'];
        else
            save_name=[fnameca(1:end-4) '_batch_pre_per.mat']
        end
        
        if read_CaImAn_mat==1
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
                ,'dCR_lick_traces','dMiss_lick_traces','dFA_lick_traces','um_per_pixel',...
                'CR_deltalogP','CRelii','FA_deltalogP','FAelii','Hit_deltalogP','Hitelii','Miss_deltalogP','Misselii',...
                'sp_deltalogP','spelii','sm_deltalogP','smelii','sp_deltalogP','spelii')
            
        end
        
        %                 ,'dCR_lick_traces','dMiss_lick_traces','dFA_lick_traces','coms','num_coms','um_per_pixel',...
        
        %Plot Ca traces for Hit, CR, Miss, FA
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
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
        
        title("Ca changes aligned to odor onset, all events")
        
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
    end
    %Plot the licks
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
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
        savefig([fnamerhd(1:end-4) '_dropc_warp_Fig5.fig'])
    else
        savefig([fnamerhd(1:end-4) '_dropc_batch_Fig5.fig'])
    end
    
    %Plot lick frequency
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
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
        savefig([fnamerhd(1:end-4) '_dropc_warp_Fig6.fig'])
    else
        savefig([fnamerhd(1:end-4) '_dropc_batch_Fig6.fig'])
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
        hFig=figure(figNo);
        plot(time_p_lick,log10(p_val_Hit_CR))
        hold on
        plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)])
        title('log(p value) for the difference in licks Hit vs CR')
        xlabel('Time (sec)')
        ylabel('log10(p value)')
        if do_warp==1
            savefig([fnamerhd(1:end-4) '_dropc_warp_Fig7.fig'])
        else
            savefig([fnamerhd(1:end-4) '_dropc_batch_Fig7.fig'])
        end
    catch
    end
     
    if handles_choice.No_images~=1
        per_file_save_name=[fnameca(1:end-4) '_batch_per_file.mat'];
    else
        per_file_save_name=[fnamerhd(1:end-4) '_batch_per_file.mat'];
    end
    save(per_file_save_name, 'handles_per_file','-v7.3')
    
    handles_per_trial.time_licks=time_licks;
    handles_per_trial.acq_rate=acq_rate;
    handles_per_trial.dt_before=dt_before;
    handles_per_trial.delta_odor=delta_odor;
    handles_per_trial.delta_odor_on_reinf_on=delta_odor_on_reinf_on;
    handles_per_trial.delta_reinf=delta_reinf;
    if handles_choice.No_images~=1
        per_trial_save_name=[fnameca(1:end-4) '_batch_per_trial.mat'];
    else
        per_trial_save_name=[fnamerhd(1:end-4) '_batch_per_trial.mat'];
    end
    save(per_trial_save_name, 'handles_per_trial','-v7.3')
    
    pffft=1;
    
end

handles_per_trial.time_licks=time_licks;
handles_per_trial.acq_rate=acq_rate;
handles_per_trial.dt_before=dt_before;
handles_per_trial.delta_odor=delta_odor;
handles_per_trial.delta_odor_on_reinf_on=delta_odor_on_reinf_on;
handles_per_trial.delta_reinf=delta_reinf;
if handles_choice.No_images~=1
    per_trial_save_name=[fnameca(1:end-4) '_batch_per_trial.mat'];
else
    per_trial_save_name=[fnamerhd(1:end-4) '_batch_per_trial.mat'];
end
save(per_trial_save_name, 'handles_per_trial','-v7.3')

pffft=1;
















