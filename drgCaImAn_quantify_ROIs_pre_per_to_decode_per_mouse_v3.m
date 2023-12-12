%drgCaImAn_quantify_ROIs_pre_per_to_decode_per_mouse_v3
close all
clear all

is_Fabio=0;
%0 Choices for Ming's go-no go processing
%1 Choices for Fabio
%2 Choices for Ming's passive

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m';'drgCaImAn_LDAfsdz_choices*.mat'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdzv2 run for ' choiceFileName '\n\n']);


%     is_mat=0;
tempDirName=['temp' choiceFileName(12:end-2)];
addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;


% fileID = fopen([choiceBatchPathName 'drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse.txt'],'w');

no_files=handles.no_files;
moving_mean_n=10;

time_periods_eu=[
    -1 0;
    3.1 4.1;
    4.5 5.5];

group_sets(1,:)=[2 6];
group_sets(2,:)=[3 7];
group_sets(3,:)=[4 8];

%Which threshold value should we use?
ii_thr=find(handles.p_threshold==1.1); %Note that I am using 1.1

%Show the per mouse graphs?
show_per_mouse=0;



switch is_Fabio
    case 0

        %Choices for Ming's go-no go processing
        no_pcorr=4;

        %groups to be shown in the zoomed figures for Ming's data
        grNo1=4; %Forward proficient
        grNo1_label='forward proficient';
        grNo2=8; %Forward proficient
        grNo2_label='reversed proficient';

        %         window_label{1}='Base';
        window_label{1}='Pre';
        window_label{2}='Odor';
        window_label{3}='Reinf';


        per_corr_set_label{1}='learning';
        per_corr_set_label{2}='intermediate';
        per_corr_set_label{3}='proficient';


    case 1

        %Choices for Fabio
        %groups to be shown in the zoomed figures for Fabio's data
        % grNo1=1; %AAAP
        % grNo1_label='AAAP';
        grNo1=2; %female bedding
        grNo1_label='female bedding';
        grNo2=4; %male bedding
        grNo2_label='male bedding';

        %Choices for Fabio's passive odorant exposure processing
        no_pcorr=1;

    case 2

        %Choices for Ming's passive
        %groups to be shown in the zoomed figures for Fabio's data
        % grNo1=1; %AAAP
        % grNo1_label='AAAP';
        grNo1=1; %passive
        grNo1_label='passive';
        grNo2=1; %passive
        grNo2_label='passive';

        %Choices for Fabio's passive odorant exposure processing
        no_pcorr=1;

        %          window_label{1}='Base';
        window_label{2}='Pre';
        window_label{3}='Odor';
        window_label{4}='Reinf';


        per_corr_set_label{1}='';
    case 3
        no_pcorr=1;

end

fprintf(1, ['\nData were processed with p value threshold = ' num2str(handles.p_threshold(ii_thr)) '\n\n']);

t_from=-10;
t_to=20;

%Now plot for each algorithm the prediction accuracy
handles_out2=[];

handles_out2.classifier_names{1}='Linear Discriminant';
handles_out2.classifier_names{2}='Support Vector Machine';
handles_out2.classifier_names{3}='Naive Bayes Classifier';
handles_out2.classifier_names{4}='Neural Network';
handles_out2.classifier_names{5}='Decision tree';
handles_out2.classifier_names{6}='Binomial glm';

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


%
% %First and last sessions per mouse
% first_last=[20 2;
%     5 15;
%     25 1;
%     13 11;
%     26 12];
%
% handles_out2.first_last=first_last;
%
% %Mouse names
% handles_out2.mouse_names{1}='GRIN1';
% handles_out2.mouse_names{2}='GRIN3';
% handles_out2.mouse_names{3}='GRIN4';
% handles_out2.mouse_names{4}='GRIN6';
% handles_out2.mouse_names{5}='GRIN7';

per_names{1}='<40%';
per_names{2}='40-65%';
per_names{3}='65-80%';
per_names{4}='>=80%';

fr_per_names{1}='Fwd <40%%';
fr_per_names{2}='Fwd 40-65%%';
fr_per_names{3}='Fwd 65-80%%';
fr_per_names{4}='Fwd >=80%%';
fr_per_names{5}='Rev <40%%';
fr_per_names{6}='Rev 40-65%%';
fr_per_names{7}='Rev 65-80%%';
fr_per_names{8}='Rev >=80%%';

if exist('no_pcorr')~=0
    if no_pcorr==1
        per_names{1}='';
    end
end

resample_time_span=[-7:0.05:15];

these_groups=unique(handles.group);

% if is_mat==0
%Process each file separately
switch is_Fabio
    case 0
        for grNo=1:no_pcorr*length(these_groups)
            handles_out2.group_no(grNo).ii_euclid=0;
            for iiMLalgo=handles.MLalgo_to_use
                for ii_out=1:length(handles.p_threshold)
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=0;
                end
            end
        end
    case 3
        for grNo=1:length(these_groups)
            handles_out2.group_no(grNo).ii_euclid=0;
            for iiMLalgo=handles.MLalgo_to_use
                for ii_out=1:length(handles.p_threshold)
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=0;
                end
            end
        end
end

%Find the mouse numbers
handles_out2.mouseNo_per_file=[];
handles_out2.mouse_names=[];

for fileNo=1:no_files
    if fileNo==1
        ii_mouse=1;
        handles_out2.mouse_names{1}=handles.mouse{1};
        handles_out2.mouseNo_per_file(1)=1;
    else
        %Find whether this mouse is already in the list
        mouse_found=0;
        for this_ii=1:length(handles_out2.mouse_names)
            if strcmp(handles.mouse{fileNo},handles_out2.mouse_names{this_ii})
                mouse_found=1;
                mouse_found_ii=this_ii;
            end
        end
        if mouse_found==0
            %This mouse's name is not in the list
            ii_mouse=ii_mouse+1;
            handles_out2.mouse_names{ii_mouse}=handles.mouse{fileNo};
            handles_out2.mouseNo_per_file(fileNo)=ii_mouse;
        else
            %This mouse's name is in the list
            handles_out2.mouseNo_per_file(fileNo)=mouse_found_ii;
        end
    end
end

mouse_numbers=[3 4 1 2]; %Note: This maps the mouse numbers to the Figure 2 numbers

ROI_accounting=[];
for mouseNo=handles_out2.mouseNo_per_file
    for grNo=1:4
        ROI_accounting.group(grNo).mouse(mouseNo).no_sessions=0;
    end
end
group_per_file=handles.group;
handles_out2.pcorr_per_file=zeros(1,length(group_per_file));
for fileNo=1:no_files
    tic
    if iscell(handles.PathName_out)
        pre_per_outPathName=handles.PathName_out{fileNo};
    else
        pre_per_outPathName=handles.PathName_out;
    end

    pre_per_FileName=handles.FileName_pre_per{fileNo};
    this_group=handles.group(fileNo);
    this_grNo=find(these_groups==this_group);

    if exist([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])~=0
        load([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])

        if isfield(handles_out.ii_out(1).handles_out,'percent_correct')
            pCorr=handles_out.ii_out(1).handles_out.percent_correct;
            handles_out2.pcorr_per_file(fileNo)=pCorr;
            ii_pCorr=1;
            if (pCorr>=40)&(pCorr<=65)
                ii_pCorr=2;
            else
                if (pCorr>65)&(pCorr<80)
                    ii_pCorr=3;
                else
                    if pCorr>=80
                        ii_pCorr=4;
                    end
                end
            end
        else
            ii_pCorr=1;
        end

        if (is_Fabio==1)||(is_Fabio==2)
            ii_pCorr=1;
        end

        switch is_Fabio
            case 0
                grNo=(this_grNo-1)*no_pcorr+ii_pCorr;
            case 3
                grNo=this_grNo;
        end
        mouseNo=handles_out2.mouseNo_per_file(fileNo);
        if  handles_out.ii_out(1).handles_out.no_spm_odor_trials>=20
            ROI_accounting.group(ii_pCorr).mouse(mouseNo).no_sessions=ROI_accounting.group(ii_pCorr).mouse(mouseNo).no_sessions+1;
            ROI_accounting.group(ii_pCorr).mouse(mouseNo).session(ROI_accounting.group(ii_pCorr).mouse(mouseNo).no_sessions).no_ROIS=handles_out.ii_out(1).handles_out.no_components;
        end

        %             for iiMLalgo=handles.MLalgo_to_use
        %                 if iiMLalgo==handles.MLalgo_to_use(1)
        %                     handles_out2.group_no(grNo).ii_euclid=handles_out2.group_no(grNo).ii_euclid+1;
        %                     ii_euclid=handles_out2.group_no(grNo).ii_euclid;
        % %                     handles_out2.group_no(grNo).dist_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.dist_euclid))=handles_out.ii_out(1).handles_out.dist_euclid-handles_out.ii_out(1).handles_out.dist_euclid_zero;
        % %                     handles_out2.group_no(grNo).KLdivergence(ii_euclid,1:length(handles_out.ii_out(1).handles_out.KLdivergence))=handles_out.ii_out(1).handles_out.KLdivergence;
        %                     handles_out2.group_no(grNo).time_span_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.time_span;
        %                     handles_out2.group_no(grNo).ii_time_span(ii_euclid,1)=length(handles_out.ii_out(1).handles_out.time_span);
        %                     handles_out2.group_no(grNo).meandFFsp(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsp;
        %                     handles_out2.group_no(grNo).meandFFsm(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsm;
        %                     handles_out2.group_no(grNo).mouseNo_euclid(ii_euclid)= handles_out2.mouseNo_per_file(fileNo);
        %                 end
        %                 for ii_out=1:length(handles_out.ii_out)
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii+1;
        %                     ii=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii;
        %                     accuracy_tr=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr;
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr(ii)=accuracy_tr;
        %
        %                     accuracy_tr_sh=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh;
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_sh(ii)=accuracy_tr_sh;
        %
        %                     accuracy_tr_sh2=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh2;
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_sh2(ii)=accuracy_tr_sh2;
        %
        %                     this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(ii,1:length(this_mean_correct_predict))=this_mean_correct_predict;
        %
        %                     this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(ii,1:length(this_mean_correct_predict))=this_mean_correct_predict;
        %
        %                     this_mean_correct_predict_sh=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh);
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh(ii,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;
        %
        %                     this_per_trial_sp_timecourse=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_sp_timecourse;
        %                     this_per_trial_sm_timecourse=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_sm_timecourse;
        %
        %                     this_prediction_mean_per_trial_sp_timecourse = movmean(this_per_trial_sp_timecourse',moving_mean_n)';
        %                     this_mean_prediction_mean_per_trial_sp_timecourse=mean(this_prediction_mean_per_trial_sp_timecourse);
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sp_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_sp_timecourse))=this_mean_prediction_mean_per_trial_sp_timecourse;
        %
        %                     this_prediction_mean_per_trial_sm_timecourse = movmean(this_per_trial_sm_timecourse',moving_mean_n)';
        %                     this_mean_prediction_mean_per_trial_sm_timecourse=mean(this_prediction_mean_per_trial_sm_timecourse);
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sm_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_sm_timecourse))=this_mean_prediction_mean_per_trial_sm_timecourse;
        %
        %                     per_trial_scores_sp=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_scores_sp;
        %                     mean_per_trial_scores_sp=zeros(size(per_trial_scores_sp,2),size(per_trial_scores_sp,3));
        %                     mean_per_trial_scores_sp(:,:)=mean(per_trial_scores_sp,1);
        %                     if ~isempty(per_trial_scores_sp)
        %                         handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii,1:2,1:size(per_trial_scores_sp,3))=mean_per_trial_scores_sp;
        %                     end
        %
        %                     per_trial_scores_sm=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_scores_sm;
        %                     mean_per_trial_scores_sm=zeros(size(per_trial_scores_sm,2),size(per_trial_scores_sm,3));
        %                     mean_per_trial_scores_sm(:,:)=mean(per_trial_scores_sm,1);
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii,1:2,1:size(per_trial_scores_sm,3))=mean_per_trial_scores_sm;
        %
        %                     handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mouseNo(ii)=handles_out2.mouseNo_per_file(fileNo);
        %                 end
        %             end
        fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
    else
        fprintf(1, ['Import for file No ' num2str(fileNo) ' failed\n'])
        if fileNo==69
            pffft=1;
        end
    end
end


%Now print the output (ROI+/-SD, no sessions)
naive_pro{1}='After reversal';
naive_pro{2}='Learning';
naive_pro{4}='Proficient';
for ii_pc=[2 4]
    all_mice_no_ROIS=[];
    for mouseNo=unique(handles_out2.mouseNo_per_file)
        noROIs=zeros(1,ROI_accounting.group(ii_pc).mouse(mouseNo).no_sessions);
        if ROI_accounting.group(ii_pc).mouse(mouseNo).no_sessions>0
            for sessionNo=1:ROI_accounting.group(ii_pc).mouse(mouseNo).no_sessions
                noROIs(sessionNo)=ROI_accounting.group(ii_pc).mouse(mouseNo).session(sessionNo).no_ROIS;
            end
            all_mice_no_ROIS=[all_mice_no_ROIS noROIs];
            fprintf(1,['For ' naive_pro{ii_pc} ' mouse number ' num2str(mouse_numbers(mouseNo)) ' mean number of ROIs ' num2str(mean(noROIs)) ' STD ' num2str(std(noROIs)) ' n= ' num2str(length(noROIs)) '\n'])
        else
            fprintf(1,['For ' naive_pro{ii_pc} ' mouse number ' num2str(mouse_numbers(mouseNo)) ' no sessions\n'])

        end

    end
    fprintf(1,['For ' naive_pro{ii_pc} ' all mice mean number of ROIs ' num2str(mean(all_mice_no_ROIS)) ' STD ' num2str(std(all_mice_no_ROIS)) ' n= ' num2str(length(all_mice_no_ROIS)) '\n'])
end
% else
%     load([choiceBatchPathName choiceFileName])
% end

% %  ii_thr=length(handles.p_threshold);
%  ii_thr=1;

pfft=1;