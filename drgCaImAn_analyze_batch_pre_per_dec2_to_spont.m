function drgCaImAn_analyze_batch_pre_per_dec2_to_spont(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_LDA_fsdz run for ' choiceFileName '\n\n']);


fileID = fopen([choiceBatchPathName 'drgCaImAn_analyze_batch_pre_per_dec2_to_spont.txt'],'w');

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.no_files;

if isfield(handles,'processing_algo')
    processing_algo=handles.processing_algo;
else
    processing_algo=1;
end

% if isfield(handles,'suffix_out')
%     suffix_out=handles.suffix_out;
% else
%     suffix_out='_dec.mat';
% end

suffix_out='_spont.mat';

if isfield(handles,'first_file')
    first_file=handles.first_file;
else
    first_file=1;
end

mouse_names=[];
mouseNo_per_file=[];
for fileNo=1:handles.no_files
    if fileNo==1
        ii_mouse=1;
        mouse_names{1}=handles.mouse{1};
        mouseNo_per_file(1)=1;
    else
        %Find whether this mouse is already in the list
        mouse_found=0;
        for this_ii=1:length(mouse_names)
            if strcmp(handles.mouse{fileNo},mouse_names{this_ii})
                mouse_found=1;
                mouse_found_ii=this_ii;
            end
        end
        if mouse_found==0
            %This mouse's name is not in the list
            ii_mouse=ii_mouse+1;
            mouse_names{ii_mouse}=handles.mouse{fileNo};
            mouseNo_per_file(fileNo)=ii_mouse;
        else
            %This mouse's name is in the list
            mouseNo_per_file(fileNo)=mouse_found_ii;
        end
    end
end

%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:handles.no_files


    %Make sure that all the files exist
    pre_per_FileName=handles.FileName_pre_per{filNum};
    if iscell(handles.PathName_pre_per)
        pre_per_PathName=handles.PathName_pre_per{filNum};
    else
        pre_per_PathName=handles.PathName_pre_per;
    end

    if exist([pre_per_PathName pre_per_FileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' pre_per_FileName ' does not exist\n'],filNum);
        all_files_present=0;
    end

end


% if exist([handles.PathName_out handles.FileName_out])==0
%     handles_out=[];
%     ii_out=0;
% else
%     load([handles.PathName_out handles.FileName_out])
%     ii_out=handles_out.last_ii_out;
%     first_file=handles_out.last_file_processed+1;
% end

figNo=0;
show_figures=1;
min_acc=0.7;
time_span=[-7:0.05:15];
%Note that some calcium imaging was taken as 20 Hz and some at ~4 Hz. Here
%we move all data to -7:0.05:15

if all_files_present==1

    %Retrieve the data saved by drgCaImAn_analyze_SVZ_entire_session_between(handles_choices)
    %Process each file separately
    moving_mean_per_session_sm_timecourse=[];
    moving_mean_per_session_sp_timecourse=[];
    these_mice=[];
    ii_mouse=0;

    between_spont_prediction=[];
    between_no_change_prediction=[];
    these_mice_sp=[];
    ii_mouse_sp=0;
    ii_mouse_nc=0;
    these_mice_nc=[];

    corr_spont_prediction_sp=[];
    corr_spont_prediction_sm=[];
    these_mice_csp=[];
    ii_mouse_csp=0;
    ii_mouse_csm=0;
    these_mice_csm=[];

    corr_spont_prediction_sp_nc=[];
    corr_spont_prediction_sm_nc=[];
    these_mice_csp_nc=[];
    ii_mouse_csp_nc=0;
    ii_mouse_csm_nc=0;
    these_mice_csm_nc=[];

    frac_sp_before=[];
    frac_sp_after=[];
    these_mice_sp_before=[];
    ii_mouse_sp_before=0;
    ii_mouse_sp_after=0;
    these_mice_sp_after=[];
    min_ba_trials=10;


    frac_pred_between=[];
    ii_mouse_fpb=0;
    these_mice_fpb=[];

    frac_pred_between_sh=[];
    ii_mouse_fpb_sh=0;
    these_mice_fpb_sh=[];


    frac_pred_sm_odor=[];
    frac_pred_sp_odor=[];
    frac_pred_sm_baseline=[];
    frac_pred_sp_baseline=[];
    ii_mouse_fspm=0;
    these_mice_fspm=[];

    frac_pred_spon_baseline=[];
    frac_pred_spon_onset=[];
    ii_mouse_fspon=0;
    these_mice_fspon=[];

    frac_pred_nc_baseline=[];
    frac_pred_nc_onset=[];
    ii_mouse_fnc=0;
    these_mice_fnc=[];

    resample_time_span=[-7:0.05:15];
    for fileNo=first_file:length(handles.FileName_pre_per)
        if sum(fileNo==handles.skip_files)==0
            tic
            first_toc=toc;

            handles_out=[];
            ii_out=0;

            pre_per_PathName=handles.PathName_pre_per{fileNo};
            pre_per_FileName=handles.FileName_pre_per{fileNo};

            load([pre_per_PathName pre_per_FileName(1:end-4) suffix_out])

            if (handles_out.handles_out3.percent_correct>=80)&(handles_out.handles_out3.odor_acc>=min_acc)

                %Prediction within trials
                this_time_span=handles_out.handles_out3.moving_time_span;

                this_in_timecourse=mean(handles_out.handles_out3.moving_mean_per_trial_sm_timecourse);
                this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                moving_mean_per_session_sm_timecourse=[moving_mean_per_session_sm_timecourse this_out_timecourse'];

                this_in_timecourse=mean(handles_out.handles_out3.moving_mean_per_trial_sp_timecourse);
                this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                moving_mean_per_session_sp_timecourse=[moving_mean_per_session_sp_timecourse this_out_timecourse'];

                ii_mouse=ii_mouse+1;
                these_mice(ii_mouse)=mouseNo_per_file(fileNo);

                %Spontaneous change in prediction between trials
                this_time_span=handles_out.handles_out3.this_time_span_spon;

                if ~isempty(this_time_span)

                    if size(handles_out.handles_out3.spontaneous_delta_prediction,1)>0
                        if size(handles_out.handles_out3.spontaneous_delta_prediction,1)==1
                            this_in_timecourse=handles_out.handles_out3.spontaneous_delta_prediction;
                        else
                            this_in_timecourse=mean(handles_out.handles_out3.spontaneous_delta_prediction);
                        end
                        this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                        between_spont_prediction=[between_spont_prediction this_out_timecourse'];
                        ii_mouse_sp=ii_mouse_sp+1;
                        these_mice_sp(ii_mouse_sp)=mouseNo_per_file(fileNo);
                    end

                    if size(handles_out.handles_out3.no_increase_delta_prediction,1)>0
                        if size(handles_out.handles_out3.no_increase_delta_prediction,1)==1
                            this_in_timecourse=handles_out.handles_out3.no_increase_delta_prediction;
                        else
                            this_in_timecourse=mean(handles_out.handles_out3.no_increase_delta_prediction);
                        end

                        this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                        between_no_change_prediction=[between_no_change_prediction this_out_timecourse'];
                        ii_mouse_nc=ii_mouse_nc+1;
                        these_mice_nc(ii_mouse_nc)=mouseNo_per_file(fileNo);
                    end


                end




                %Correlation dFF at tie time of between trial spontaneous change in prediction and dFF within
                %trials
                this_time_span=handles_out.handles_out3.this_time_span_corr;

                if ~isempty(this_time_span)

                    if size(handles_out.handles_out3.corr_sp,1)>0
                        if size(handles_out.handles_out3.corr_sp,1)==1
                            this_in_timecourse=handles_out.handles_out3.corr_sp;
                        else
                            this_in_timecourse=mean(handles_out.handles_out3.corr_sp);
                        end
                        this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                        corr_spont_prediction_sp=[corr_spont_prediction_sp this_out_timecourse'];
                        ii_mouse_csp=ii_mouse_csp+1;
                        these_mice_csp(ii_mouse_csp)=mouseNo_per_file(fileNo);
                    end

                    if size(handles_out.handles_out3.corr_sm,1)>0
                        if size(handles_out.handles_out3.corr_sm,1)==1
                            this_in_timecourse=handles_out.handles_out3.corr_sm;
                        else
                            this_in_timecourse=mean(handles_out.handles_out3.corr_sm);
                        end

                        this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                        corr_spont_prediction_sm=[corr_spont_prediction_sm this_out_timecourse'];
                        ii_mouse_csm=ii_mouse_csm+1;
                        these_mice_csm(ii_mouse_csm)=mouseNo_per_file(fileNo);
                    end


                end

                %Correlation dFF at tie time of between trial no change in prediction and dFF within
                %trial
                this_time_span=handles_out.handles_out3.this_time_span_nc_corr;

                if ~isempty(this_time_span)

                    if size(handles_out.handles_out3.nc_corr_sp,1)>0
                        if size(handles_out.handles_out3.nc_corr_sp,1)==1
                            this_in_timecourse=handles_out.handles_out3.nc_corr_sp;
                        else
                            this_in_timecourse=mean(handles_out.handles_out3.nc_corr_sp);
                        end
                        this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                        corr_spont_prediction_sp_nc=[corr_spont_prediction_sp_nc this_out_timecourse'];
                        ii_mouse_csp_nc=ii_mouse_csp_nc+1;
                        these_mice_csp_nc(ii_mouse_csp_nc)=mouseNo_per_file(fileNo);
                    end

                    if size(handles_out.handles_out3.nc_corr_sm,1)>0
                        if size(handles_out.handles_out3.nc_corr_sm,1)==1
                            this_in_timecourse=handles_out.handles_out3.nc_corr_sm;
                        else
                            this_in_timecourse=mean(handles_out.handles_out3.nc_corr_sm);
                        end


                        this_out_timecourse=drgCaImAnTimeResample(this_time_span,this_in_timecourse,resample_time_span);
                        corr_spont_prediction_sm_nc=[corr_spont_prediction_sm_nc this_out_timecourse'];
                        ii_mouse_csm_nc=ii_mouse_csm_nc+1;
                        these_mice_csm_nc(ii_mouse_csm_nc)=mouseNo_per_file(fileNo);

                    end


                end

                %Do accounting for before and after trials
                this_sp_before=handles_out.handles_out3.before_spontaneous_spm;
                if length(this_sp_before)>=min_ba_trials
                    ii_mouse_sp_before=ii_mouse_sp_before+1;
                    frac_sp_before(ii_mouse_sp_before)=sum(this_sp_before)/length(this_sp_before);
                    these_mice_sp_before(ii_mouse_sp_before)=mouseNo_per_file(fileNo);
                end

                this_sp_after=handles_out.handles_out3.after_spontaneous_spm;
                if length(this_sp_after)>=min_ba_trials
                    ii_mouse_sp_after=ii_mouse_sp_after+1;
                    frac_sp_after(ii_mouse_sp_after)=sum(this_sp_after)/length(this_sp_after);
                    these_mice_sp_after(ii_mouse_sp_after)=mouseNo_per_file(fileNo);
                end



                pffft=1;

                %Now get the fraction of prediction per time interval

                %Between for entire between interval
                if ~isempty(handles_out.handles_out3.frac_pred_between)
                    frac_pred_between=[frac_pred_between mean(handles_out.handles_out3.frac_pred_between)];
                    frac_pred_between_sh=[frac_pred_between_sh mean(handles_out.handles_out3.frac_pred_between_sh)];
                    ii_mouse_fpb=ii_mouse_fpb+1;
                    these_mice_fpb(ii_mouse_fpb)=mouseNo_per_file(fileNo);
                end

                %Within for S+ and S-
                if ~isempty(handles_out.handles_out3.frac_pred_sm_odor)
                    frac_pred_sm_odor=[frac_pred_sm_odor mean(handles_out.handles_out3.frac_pred_sm_odor)];
                    frac_pred_sp_odor=[frac_pred_sp_odor mean(handles_out.handles_out3.frac_pred_sp_odor)];
                    frac_pred_sm_baseline=[frac_pred_sm_baseline mean(handles_out.handles_out3.frac_pred_sm_baseline)];
                    frac_pred_sp_baseline=[frac_pred_sp_baseline mean(handles_out.handles_out3.frac_pred_sp_baseline)];

                    ii_mouse_fspm=ii_mouse_fspm+1;
                    these_mice_fspm(ii_mouse_fspm)=mouseNo_per_file(fileNo);
                end

                %Between for spontaneous betweeen intervals
                if ~isempty(handles_out.handles_out3.frac_pred_spon_baseline)
                    frac_pred_spon_baseline=[frac_pred_spon_baseline mean(handles_out.handles_out3.frac_pred_spon_baseline)];
                    frac_pred_spon_onset=[frac_pred_spon_onset mean(handles_out.handles_out3.frac_pred_spon_onset)];
                    ii_mouse_fspon=ii_mouse_fspon+1;
                    these_mice_fspon(ii_mouse_fspon)=mouseNo_per_file(fileNo);
                end

                %Between for no change betweeen intervals
                if ~isempty(handles_out.handles_out3.frac_pred_nc_onset)
                    frac_pred_nc_baseline=[frac_pred_nc_baseline mean(handles_out.handles_out3.frac_pred_nc_baseline)];
                    frac_pred_nc_onset=[frac_pred_nc_onset mean(handles_out.handles_out3.frac_pred_nc_onset)];
                    ii_mouse_fnc=ii_mouse_fnc+1;
                    these_mice_fnc(ii_mouse_fnc)=mouseNo_per_file(fileNo);
                end

            end
            pffft=1;
            %             handles_choices.pre_per_PathName=pre_per_PathName;
            %             handles_choices.pre_per_FileName=pre_per_FileName;
            %             handles_choices.processing_algorithm=handles.processing_algorithm;
            %             handles_choices.MLalgo_to_use=6; %Please note I am doing glm
            %             handles_choices.dt_p_threshold=handles.dt_p_threshold;
            %             handles_choices.show_figures=handles.show_figures;
            %             handles_choices.post_time=handles.post_time;
            %             handles_choices.k_fold=handles.k_fold;
            %             handles_choices.post_shift=handles.post_shift;
            %             handles_choices.pre_time=handles.pre_time;
            %             handles_choices.ii_cost=handles.ii_cost;
            %             handles_choices.PathName_Out=handles.PathName_out{fileNo};
            %
            %             ii_thr=4;
            %
            %             handles_choices.p_threshold=handles.p_threshold(ii_thr);
            %
            %             ii_out=ii_out+1;
            %             handles_out.handles_choices=handles_choices;
            %             handles_out.grNo=handles.group(fileNo);
            %             handles_out.fileNo=fileNo;
            %
            %             start_toc=toc;
            %
            %             handles_out.handles_out.handles_out3=drgCaImAn_analyze_SVZ_entire_session_between(handles_choices);
            %
            %             fprintf(1, ['Data processed for file number %d, condition number %d\n'],fileNo,ii_thr);
            %
            %             fprintf(1,'Processing time for drgCaImAn_pre_per_to_LDA_fsdz_new %d hours\n',(toc-start_toc)/(60*60));
            %
            %
            %             fprintf(1,'\n\nProcessing time for file No %d is %d hours\n',fileNo,(toc-first_toc)/(60*60));
            %
            %             %Save output file
            %             handles_out.last_file_processed=fileNo;
            %             handles_out.last_ii_out=ii_out;
            %             handles_out.handles=handles;
            %             save([pre_per_FileName suffix_out],'handles_out','handles_choices','-v7.3')
        end
    end

    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));




    %Plot the prediction for S+ and S- and the accuracy
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    CIsm = bootci(1000, @mean, moving_mean_per_session_sm_timecourse');
    meansm=mean(moving_mean_per_session_sm_timecourse',1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(resample_time_span',mean(moving_mean_per_session_sm_timecourse,2)', CIsm', 'cmap',[158/255 31/255 99/255]);

    CIsp = bootci(1000, @mean, moving_mean_per_session_sp_timecourse');
    meansp=mean(moving_mean_per_session_sp_timecourse',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;


    [hlsp, hpsp] = boundedline(resample_time_span',mean(moving_mean_per_session_sp_timecourse,2)', CIsp', 'cmap',[0 114/255 178/255]);

    plot(resample_time_span',mean(moving_mean_per_session_sm_timecourse,2)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
    plot(resample_time_span',mean(moving_mean_per_session_sp_timecourse,2)', '-','Color',[0 114/255 178/255],'DisplayName','S+');

    text(30,0.75,'S-','Color',[158/255 31/255 99/255])
    text(30,0.85,'S+','Color',[0 114/255 178/255])

    ylim([0 1.1])
    this_ylim=ylim;

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
    title(['Label prediction'])
    xlabel('Time(sec)')
    ylabel('Label prediction, S+=1, S-=0')








    %Plot the prediction between trials
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

    hold on


    if ~isempty(between_no_change_prediction')
        CIsm = bootci(1000, @mean, between_no_change_prediction');
        meansm=mean(between_no_change_prediction',1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;

        [hlsm, hpsm] = boundedline(resample_time_span',mean(between_no_change_prediction,2)', CIsm', 'cmap',[158/255 31/255 99/255]);
    end

    if ~isempty(between_spont_prediction')
        CIsm = bootci(1000, @mean, between_spont_prediction');
        meansm=mean(between_spont_prediction',1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;

        [hlsm, hpsm] = boundedline(resample_time_span',mean(between_spont_prediction,2)', CIsm', 'cmap',[0 114/255 178/255]);
    end

    %     for ii=1:size(spontaneous_delta_prediction,1)
    %         plot(this_time_span',spontaneous_delta_prediction(ii,:),  'Color',[150/255 150/255 150/255])
    %     end
    if ~isempty(between_no_change_prediction)
        plot(resample_time_span',mean(between_no_change_prediction,2)', 'Color',[158/255 31/255 99/255]);
    end

    if ~isempty(between_spont_prediction)
        plot(resample_time_span',mean(between_spont_prediction,2)', 'Color',[0 114/255 178/255]);
    end

    text(10,0.75,'No change','Color',[158/255 31/255 99/255])
    text(10,0.85,'Spontaneous','Color',[0 114/255 178/255])

    ylim([0 1.1])
    this_ylim=ylim;
    plot([0 0],[this_ylim],'-k')
    xlim([-7 15])
    title(['Label prediction between trials'])
    xlabel('Time(sec)')
    ylabel('Label prediction')


    %Plot the correlation for spontaneous S+ prediction vs. S+ and S-
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    CIsm = bootci(1000, @mean, corr_spont_prediction_sm');
    meansm=mean(corr_spont_prediction_sm',1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(resample_time_span',mean(corr_spont_prediction_sm,2), CIsm', 'cmap',[158/255 31/255 99/255]);

    CIsp = bootci(1000, @mean, corr_spont_prediction_sp');
    meansp=mean(corr_spont_prediction_sp',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;


    [hlsp, hpsp] = boundedline(resample_time_span',mean(corr_spont_prediction_sp,2), CIsp', 'cmap',[0 114/255 178/255]);

    plot(resample_time_span',mean(corr_spont_prediction_sm,2)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
    plot(resample_time_span',mean(corr_spont_prediction_sp,2)', '-','Color',[0 114/255 178/255],'DisplayName','S+');


    ylim([-0.04 0.1])
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
    title(['Correlation for spontaneous between vs. actual within dFF'])
    xlabel('Time(sec)')
    ylabel('Rho')



    %Plot the correlation for no change S+ prediction vs. S+ and S-
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    CIsm = bootci(1000, @mean, corr_spont_prediction_sm_nc');
    meansm=mean(corr_spont_prediction_sm_nc',1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(resample_time_span',mean(corr_spont_prediction_sm_nc,2), CIsm', 'cmap',[158/255 31/255 99/255]);

    CIsp = bootci(1000, @mean, corr_spont_prediction_sp_nc');
    meansp=mean(corr_spont_prediction_sp_nc',1);
    CIsp(1,:)=meansp-CIsp(1,:);
    CIsp(2,:)=CIsp(2,:)-meansp;


    [hlsp, hpsp] = boundedline(resample_time_span',mean(corr_spont_prediction_sp_nc,2), CIsp', 'cmap',[0 114/255 178/255]);

    plot(resample_time_span',mean(corr_spont_prediction_sm_nc,2)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
    plot(resample_time_span',mean(corr_spont_prediction_sp_nc,2)', '-','Color',[0 114/255 178/255],'DisplayName','S+');


    ylim([-0.04 0.1])
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
    title(['Correlation for no change between vs. actual within dFF'])
    xlabel('Time(sec)')
    ylabel('Rho')



    %Bar graph plot for sp before and after


    figNo = figNo + 1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .3 .25])

    hold on

    edges=[0:0.05:1];
    rand_offset=0.8;

    bar_offset=0;



    id_ii=0;
    input_data=[];

    glm_frac=[];
    glm_ii=0;

    %Fraction of S- before
    bar_offset=bar_offset+1;

    bar(bar_offset,mean(1-frac_sp_before),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])

    if length(frac_sp_before)>2
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(1-frac_sp_before...
            ,edges,bar_offset,rand_offset,'k','k',3);
    end

    id_ii=id_ii+1;
    input_data(id_ii).data=1-frac_sp_before;
    input_data(id_ii).description=['S-'];

    %Fraction of S+ before
    bar_offset=bar_offset+1;

    bar(bar_offset,mean(frac_sp_before),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])

    if length(frac_sp_before)>2
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(frac_sp_before...
            ,edges,bar_offset,rand_offset,'k','k',3);
    end

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_sp_before;
    input_data(id_ii).description=['S+'];

    for mouseNo=1:max(these_mice_sp_before)
        this_frac_sp_before=mean(frac_sp_before(these_mice_sp_before==mouseNo));
        this_frac_sm_before=mean(1-frac_sp_before(these_mice_sp_before==mouseNo));
        plot([bar_offset-1 bar_offset],[this_frac_sm_before this_frac_sp_before],'LineWidth',3,'Color',[70/255 70/255 70/255])

        glm_frac.data(glm_ii+1:glm_ii+sum(these_mice_sp_before==mouseNo))=1-frac_sp_before(these_mice_sp_before==mouseNo);
        glm_frac.spm(glm_ii+1:glm_ii+sum(these_mice_sp_before==mouseNo))=zeros(1,sum(these_mice_sp_before==mouseNo));
        glm_frac.mice(glm_ii+1:glm_ii+sum(these_mice_sp_before==mouseNo))=mouseNo*ones(1,sum(these_mice_sp_before==mouseNo));
        glm_ii=glm_ii+sum(these_mice_sp_before==mouseNo);

        glm_frac.data(glm_ii+1:glm_ii+sum(these_mice_sp_before==mouseNo))=frac_sp_before(these_mice_sp_before==mouseNo);
        glm_frac.spm(glm_ii+1:glm_ii+sum(these_mice_sp_before==mouseNo))=ones(1,sum(these_mice_sp_before==mouseNo));
        glm_frac.mice(glm_ii+1:glm_ii+sum(these_mice_sp_before==mouseNo))=mouseNo*ones(1,sum(these_mice_sp_before==mouseNo));
        glm_ii=glm_ii+sum(these_mice_sp_before==mouseNo);
    end

    xticks([1 2])
    xticklabels({'S-','S+'})

    title(['Fraction of S+ or S- trials before spontaneous S+ prediction'])
    ylabel('Fraction')
    ylim([0 0.7])
    xlim([0 3])
    pffft=1;


    %Perform the glm
    fprintf(1, ['\nglm for for trials before\n'])
    fprintf(fileID, ['\nglm for for trials before\n']);

    tbl = table(glm_frac.data',glm_frac.spm',glm_frac.mice',...
        'VariableNames',{'fraction_before','sp_vs_sm','mice'});
    mdl = fitglm(tbl,'fraction_before~sp_vs_sm+mice+sp_vs_sm*mice'...
        ,'CategoricalVars',[2,3])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);

    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for trials before\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for trials before\n']);


    [output_data] = drgMutiRanksumorTtest(input_data, fileID,0);


    pffft=1;

    %Now plot a bar graph  comparing prediction between S= and S- trials
    id_ii=0;
    input_data=[];

    glm_spm=[];
    glm_ii=0;

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

    %S- baseline
    bar(bar_offset,mean(frac_pred_sm_baseline),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sm_baseline...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_sp_baseline),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sp_baseline...
        ,edges,bar_offset,rand_offset,'k','k',3);

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_sm_baseline;
    input_data(id_ii).description=['S- base'];

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_sp_baseline;
    input_data(id_ii).description=['S+ base'];

    for mouseNo=1:max(these_mice_fspm)

        plot([bar_offset-1 bar_offset],[mean(frac_pred_sm_baseline(these_mice_fspm==mouseNo)) mean(frac_pred_sp_baseline(these_mice_fspm==mouseNo))],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=frac_pred_sm_baseline(these_mice_fspm==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=zeros(1,sum(these_mice_fspm==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=zeros(1,sum(these_mice_fspm==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=mouseNo*ones(1,sum(these_mice_fspm==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fspm==mouseNo);

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=frac_pred_sp_baseline(these_mice_fspm==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=ones(1,sum(these_mice_fspm==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=zeros(1,sum(these_mice_fspm==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=mouseNo*ones(1,sum(these_mice_fspm==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fspm==mouseNo);

    end

    bar_offset=bar_offset+2;

    bar(bar_offset,mean(frac_pred_sm_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sm_odor...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_sp_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_sp_odor...
        ,edges,bar_offset,rand_offset,'k','k',3);

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_sm_odor;
    input_data(id_ii).description=['S- odor'];

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_sp_odor;
    input_data(id_ii).description=['S+ odor'];

    for mouseNo=1:max(these_mice_fspm)

        plot([bar_offset-1 bar_offset],[mean(frac_pred_sm_odor(these_mice_fspm==mouseNo)) mean(frac_pred_sp_odor(these_mice_fspm==mouseNo))],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=frac_pred_sm_odor(these_mice_fspm==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=zeros(1,sum(these_mice_fspm==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=ones(1,sum(these_mice_fspm==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=mouseNo*ones(1,sum(these_mice_fspm==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fspm==mouseNo);

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=frac_pred_sp_odor(these_mice_fspm==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=ones(1,sum(these_mice_fspm==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=ones(1,sum(these_mice_fspm==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fspm==mouseNo))=mouseNo*ones(1,sum(these_mice_fspm==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fspm==mouseNo);
    end


    xticks([0 1 3 4])
    xticklabels({'base S-', 'base S+', 'odor S-', 'odor S+'})

    ylabel('Prediction')



    title(['Prediction for S+ and S- trials'])


    %Perform the glm
    fprintf(1, ['\nglm for prediction for S+ and S- trials\n'])
    fprintf(fileID, ['\nglm prediction for S+ and S- trials\n']);

    tbl = table(glm_spm.data',glm_spm.spm',glm_spm.epoch',glm_spm.mice',...
        'VariableNames',{'prediction','sp_vs_sm','time_interval','mice'});
    mdl = fitglm(tbl,'prediction~sp_vs_sm+mice+time_interval+sp_vs_sm*mice*time_interval'...
        ,'CategoricalVars',[2,3,4])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);

    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for prediction for S+ and S- trialse\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for prediction for S+ and S- trials\n']);

    %Now plot a bar graph  comparing prediction between trials
    id_ii=0;
    input_data=[];

    glm_spm=[];
    glm_ii=0;

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

    %S- baseline
    bar(bar_offset,mean(frac_pred_nc_baseline),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_nc_baseline...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_spon_baseline),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_spon_baseline...
        ,edges,bar_offset,rand_offset,'k','k',3);

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_nc_baseline;
    input_data(id_ii).description=['No change base'];

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_spon_baseline;
    input_data(id_ii).description=['Spontaneous base'];

    for mouseNo=1:max(these_mice_fspon)

        plot([bar_offset-1 bar_offset],[mean(frac_pred_nc_baseline(these_mice_fnc==mouseNo)) mean(frac_pred_spon_baseline(these_mice_fspon==mouseNo))],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=frac_pred_nc_baseline(these_mice_fnc==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=zeros(1,sum(these_mice_fnc==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=zeros(1,sum(these_mice_fnc==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=mouseNo*ones(1,sum(these_mice_fnc==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fnc==mouseNo);

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=frac_pred_spon_baseline(these_mice_fspon==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=ones(1,sum(these_mice_fspon==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=zeros(1,sum(these_mice_fspon==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=mouseNo*ones(1,sum(these_mice_fspon==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fspon==mouseNo);

    end

    bar_offset=bar_offset+2;

    bar(bar_offset,mean(frac_pred_nc_onset),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_nc_onset...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_spon_onset),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_spon_onset...
        ,edges,bar_offset,rand_offset,'k','k',3);

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_nc_onset;
    input_data(id_ii).description=['No change onset'];

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_spon_onset;
    input_data(id_ii).description=['Spontaneous onset'];
 
    for mouseNo=1:max(these_mice_fspm)

        plot([bar_offset-1 bar_offset],[mean(frac_pred_nc_onset(these_mice_fnc==mouseNo)) mean(frac_pred_spon_onset(these_mice_fspon==mouseNo))],'-','Color',[0/255 0/255 0/255],'LineWidth',2 )

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=frac_pred_nc_onset(these_mice_fnc==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=zeros(1,sum(these_mice_fnc==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=ones(1,sum(these_mice_fnc==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fnc==mouseNo))=mouseNo*ones(1,sum(these_mice_fnc==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fnc==mouseNo);

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=frac_pred_spon_onset(these_mice_fspon==mouseNo);
        glm_spm.spm(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=ones(1,sum(these_mice_fspon==mouseNo));
        glm_spm.epoch(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=ones(1,sum(these_mice_fspon==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fspon==mouseNo))=mouseNo*ones(1,sum(these_mice_fspon==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fspon==mouseNo);
    end


    xticks([0 1 3 4])
    xticklabels({'base NC', 'base SP', 'onset NC', 'onset SP'})

    ylabel('Prediction')



    title(['Prediction for between trials'])


    %Perform the glm
    fprintf(1, ['\nglm for prediction for between trials\n'])
    fprintf(fileID, ['\nglm prediction for between trials\n']);

    tbl = table(glm_spm.data',glm_spm.spm',glm_spm.epoch',glm_spm.mice',...
        'VariableNames',{'prediction','sp_vs_sm','time_interval','mice'});
    mdl = fitglm(tbl,'prediction~sp_vs_sm+mice+time_interval+sp_vs_sm*mice*time_interval'...
        ,'CategoricalVars',[2,3,4])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);

    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for prediction for between trials\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for prediction for between trials\n']);

    %Now plot a bar graph comparing entire between interval with shuffled
    id_ii=0;
    input_data=[];

    glm_spm=[];
    glm_ii=0;

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

    %S- baseline
    bar(bar_offset,mean(frac_pred_between),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_between...
        ,edges,bar_offset,rand_offset,'k','k',3);

    bar_offset=bar_offset+1;
    bar(bar_offset,mean(frac_pred_between_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
    %Violin plot
    [mean_out, CIout]=drgViolinPoint(frac_pred_between_sh...
        ,edges,bar_offset,rand_offset,'k','k',3);

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_between;
    input_data(id_ii).description=['Between'];

    id_ii=id_ii+1;
    input_data(id_ii).data=frac_pred_between_sh;
    input_data(id_ii).description=['Between shuffled'];

    for mouseNo=1:max(these_mice_fpb)

        plot([bar_offset-1 bar_offset],[mean(frac_pred_between(these_mice_fpb==mouseNo)) mean(frac_pred_between_sh(these_mice_fspon==mouseNo))],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fpb==mouseNo))=frac_pred_between(these_mice_fpb==mouseNo);
        glm_spm.shuffled(glm_ii+1:glm_ii+sum(these_mice_fpb==mouseNo))=zeros(1,sum(these_mice_fpb==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fpb==mouseNo))=mouseNo*ones(1,sum(these_mice_fpb==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fpb==mouseNo);

        glm_spm.data(glm_ii+1:glm_ii+sum(these_mice_fpb==mouseNo))=frac_pred_between_sh(these_mice_fpb==mouseNo);
        glm_spm.shuffled(glm_ii+1:glm_ii+sum(these_mice_fpb==mouseNo))=ones(1,sum(these_mice_fpb==mouseNo));
        glm_spm.mice(glm_ii+1:glm_ii+sum(these_mice_fpb==mouseNo))=mouseNo*ones(1,sum(these_mice_fpb==mouseNo));
        glm_ii=glm_ii+sum(these_mice_fpb==mouseNo);

    end


    xticks([0 1])
    xticklabels({'Between', 'Shuffled'})

    ylabel('Prediction')



    title(['Prediction for all of between'])


    %Perform the glm
    fprintf(1, ['\nglm for prediction for all of between\n'])
    fprintf(fileID, ['\nglm prediction for all of between\n']);

    tbl = table(glm_spm.data',glm_spm.shuffled',glm_spm.mice',...
        'VariableNames',{'prediction','shuffled','mice'});
    mdl = fitglm(tbl,'prediction~shuffled+mice+shuffled*mice'...
        ,'CategoricalVars',[2,3])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);

    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for prediction for all of between\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for prediction all of between\n']);

    pffft=1;
end

fclose(fileID)





