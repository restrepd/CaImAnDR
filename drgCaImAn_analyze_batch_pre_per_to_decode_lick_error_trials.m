%drgCaImAn_analyze_batch_pre_per_to_decode_lick_error_trials
close all
clear all

is_Fabio=3;
%0 Choices for Ming's go-no go processing
%1 Choices for Fabio
%2 Choices for Ming's passive
 
[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_lickdec_choices*.m';'drgCaImAn_LDAfsdz_choices*.mat'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_analyze_batch_pre_per_to_decode_lick_error_trials run for ' choiceFileName '\n\n']);


if strcmp(choiceFileName(end-3:end),'.mat')
    is_mat=1;
    tempDirName=['temp' choiceFileName(12:end-4)];
    addpath(choiceBatchPathName)
    eval(['handles=' choiceFileName(1:end-4) ';'])
    handles.choiceFileName=[choiceFileName(1:end-4) '.m'];
    handles.choiceBatchPathName=choiceBatchPathName;
else
    is_mat=0;
    tempDirName=['temp' choiceFileName(12:end-2)];
    addpath(choiceBatchPathName)
    eval(['handles=' choiceFileName(1:end-2) ';'])
    handles.choiceFileName=choiceFileName;
    handles.choiceBatchPathName=choiceBatchPathName;
end

fileID = fopen([choiceBatchPathName 'drgCaImAn_analyze_batch_pre_per_to_decode_lick_error_trials.txt'],'w');

no_files=handles.no_files;
moving_mean_n=10;

lick_frac_thr=0.1; %This is a lick fraction threshold used to classify a response area as lick vs no lick


time_periods_eu=[
            -1 0;
            3.1 4.1;
            4.5 5.5];

% time_periods_eu=[
%             0 2;
%             2 4]; %Note: Here we are interested in the two response areas where the animal must lick


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

        window_label{1}='Base';
        window_label{2}='PreFV';
        window_label{3}='PreOdor';
        window_label{4}='Odor';


        per_corr_set_label{1}='naive';
        per_corr_set_label{2}='intermediate';
        per_corr_set_label{3}='proficient';


        group_sets(1,:)=[2 6];
group_sets(2,:)=[3 7];
group_sets(3,:)=[4 8];

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

        window_label{1}='Base';
        window_label{2}='PreFV';
        window_label{3}='PreOdor';
        window_label{4}='Odor';


        per_corr_set_label{1}='';

    case 3
        %lick error processing
        %Choices for Ming's go-no go processing
        no_pcorr=1;

        %groups to be shown in the zoomed figures for Ming's data
        grNo1=1; %Forward proficient
        grNo1_label='forward proficient';
        grNo2=2; %Forward proficient
        grNo2_label='reversed proficient';

        window_label{1}='Base';
        window_label{2}='PreFV';
        window_label{3}='PreOdor';
        window_label{4}='Odor';


        per_corr_set_label{1}='naive';
        per_corr_set_label{2}='intermediate';
        per_corr_set_label{3}='proficient';

        group_sets(1,:)=[1 2];
        group_sets(2,:)=[1 2];
        group_sets(3,:)=[1 2];


end

fprintf(1, ['\nData were processed with p value threshold = ' num2str(handles.p_threshold(ii_thr)) '\n\n']);



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

if no_pcorr==1
    per_names{1}='';
end

resample_time_span=[-7:0.05:15];

these_groups=unique(handles.group);

if is_mat==0
    %Process each file separately
    for grNo=1:no_pcorr*length(these_groups)
        handles_out2.ii_euclid=0;
        for iiMLalgo=handles.MLalgo_to_use
            for ii_out=1:length(handles.p_threshold)
                handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).ii=0;
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

    group_per_file=handles.group;
    handles_out2.pcorr_per_file=zeros(1,length(group_per_file));
    first_file=1;
    for fileNo=1:no_files
        tic
        if iscell(handles.PathName_out)
            pre_per_outPathName=handles.PathName_out{fileNo};
        else
            pre_per_outPathName=handles.PathName_out;
        end

        pre_per_FileName=handles.FileName_pre_per{fileNo};
        %         grNo=handles.group(fileNo);
        %         this_grNo=find(these_groups==this_group);

        if exist([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])~=0
            load([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])
            this_handles_choices.pre_perFileName=[pre_per_FileName(1:end-4) handles.suffix_out];
            this_handles_choices.pre_perPathName=pre_per_outPathName;
            this_handles_choices.show_images=0;
            handles_outdl=drgCaImAn_analyze_one_decode_licks_entire_sessionv2(this_handles_choices);

            if first_file==1
                standard_time_span_ii=length(handles_outdl.trimmed_time_span);
                time_span=handles_outdl.trimmed_time_span;
                first_file=0;
            end

            %Note, for one file the time span was not the same
            if standard_time_span_ii==length(handles_outdl.trimmed_time_span)


                %             if isfield(handles_out.ii_out(1).handles_out,'percent_correct')
                %                 pCorr=handles_out.ii_out(1).handles_out.percent_correct;
                %                 handles_out2.pcorr_per_file(fileNo)=pCorr;
                %                 ii_pCorr=1;
                %                 if (pCorr>=40)&(pCorr<=65)
                %                     ii_pCorr=2;
                %                 else
                %                     if (pCorr>65)&(pCorr<80)
                %                         ii_pCorr=3;
                %                     else
                %                         if pCorr>=80
                %                             ii_pCorr=4;
                %                         end
                %                     end
                %                 end
                %             else
                %                 ii_pCorr=1;
                %             end
                %
                %             if (is_Fabio==1)||(is_Fabio==2)
                %                 ii_pCorr=1;
                %             end
                %
                %             grNo=(this_grNo-1)*no_pcorr+ii_pCorr;



                for iiMLalgo=handles.MLalgo_to_use
                    %                 if iiMLalgo==handles.MLalgo_to_use(1)
                    %                     handles_out2.ii_euclid=handles_out2.ii_euclid+1;
                    %                     ii_euclid=handles_out2.ii_euclid;
                    % %                     handles_out2.dist_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.dist_euclid))=handles_out.ii_out(1).handles_out.dist_euclid-handles_out.ii_out(1).handles_out.dist_euclid_zero;
                    % %                     handles_out2.KLdivergence(ii_euclid,1:length(handles_out.ii_out(1).handles_out.KLdivergence))=handles_out.ii_out(1).handles_out.KLdivergence;
                    %                     handles_out2.time_span_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.time_span;
                    %                     handles_out2.ii_time_span(ii_euclid,1)=length(handles_out.ii_out(1).handles_out.time_span);
                    %                     handles_out2.meandFFsp(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsp;
                    %                     handles_out2.meandFFsm(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsm;
                    %                     handles_out2.mouseNo_euclid(ii_euclid)= handles_out2.mouseNo_per_file(fileNo);
                    %                 end


                    ii_out=1;
                    handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).ii=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).ii+1;
                    ii=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).ii;
                    %                     accuracy_tr=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr;
                    %                     handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr(ii)=accuracy_tr;
                    %
                    %                     accuracy_tr_sh=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh;
                    %                     handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_sh(ii)=accuracy_tr_sh;
                    %
                    %                     accuracy_tr_sh2=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh2;
                    %                     handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_sh2(ii)=accuracy_tr_sh2;

                    handles_out2.ii_time_span(ii)=length(handles_outdl.trimmed_time_span);
                    handles_out2.time_span_euclid(ii,1:length(handles_outdl.trimmed_time_span))=handles_outdl.trimmed_time_span;

                    this_mean_correct_predict=mean(handles_outdl.this_correct_predict);
                    handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(ii,1:length(this_mean_correct_predict))=this_mean_correct_predict;

                    this_mean_correct_predict_sh=mean(handles_outdl.this_correct_predict_sh);
                    handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh(ii,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;

                    %S+ prediction
                    this_prediction_mean_per_trial_sp_timecourse = handles_outdl.this_prediction_mean_per_trial_sp_timecourse;
                    this_mean_prediction_mean_per_trial_sp_timecourse=mean(this_prediction_mean_per_trial_sp_timecourse);
                    handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sp_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_sp_timecourse))=this_mean_prediction_mean_per_trial_sp_timecourse;

                    %Hits
                    these_sp_hits=handles_outdl.st_hit_per_trial((handles_outdl.st_hit_per_trial==1)|(handles_outdl.st_miss_per_trial==1));
                    if  sum(these_sp_hits)>0
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).hit_included(ii)=1;

                        %Prediction
                        this_prediction_mean_per_trial_hit_timecourse = handles_outdl.this_prediction_mean_per_trial_sp_timecourse(logical(these_sp_hits),:);
                        if size(this_prediction_mean_per_trial_hit_timecourse,1)~=1
                            this_mean_prediction_mean_per_trial_hit_timecourse=mean(this_prediction_mean_per_trial_hit_timecourse);
                        else
                            this_mean_prediction_mean_per_trial_hit_timecourse=this_prediction_mean_per_trial_hit_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_hit_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_hit_timecourse))=this_mean_prediction_mean_per_trial_hit_timecourse;

                        %licks
                        this_mean_per_trial_hit_lick_timecourse = handles_outdl.this_per_trial_mean_sp_licks_post(logical(these_sp_hits),:);
                        if size(this_mean_per_trial_hit_lick_timecourse,1)~=1
                            this_mean_mean_per_trial_hit_lick_timecourse=mean(this_mean_per_trial_hit_lick_timecourse);
                        else
                            this_mean_mean_per_trial_hit_lick_timecourse=this_mean_per_trial_hit_lick_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_hit_lick_timecourse(ii,1:length(this_mean_mean_per_trial_hit_lick_timecourse))=this_mean_mean_per_trial_hit_lick_timecourse;
                    else
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).hit_included(ii)=0;
                    end

                    %Miss
                    these_sp_miss=handles_outdl.st_miss_per_trial((handles_outdl.st_hit_per_trial==1)|(handles_outdl.st_miss_per_trial==1));
                    if  sum(these_sp_miss)>0
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).miss_included(ii)=1;

                        %Prediction
                        this_prediction_mean_per_trial_miss_timecourse = handles_outdl.this_prediction_mean_per_trial_sp_timecourse(logical(these_sp_miss),:);
                        if size(this_prediction_mean_per_trial_miss_timecourse,1)~=1
                            this_mean_prediction_mean_per_trial_miss_timecourse=mean(this_prediction_mean_per_trial_miss_timecourse);
                        else
                            this_mean_prediction_mean_per_trial_miss_timecourse=this_prediction_mean_per_trial_miss_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_miss_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_miss_timecourse))=this_mean_prediction_mean_per_trial_miss_timecourse;

                        %licks
                        this_mean_per_trial_miss_lick_timecourse = handles_outdl.this_per_trial_mean_sp_licks_post(logical(these_sp_miss),:);
                        if size(this_mean_per_trial_miss_lick_timecourse,1)~=1
                            this_mean_mean_per_trial_miss_lick_timecourse=mean(this_mean_per_trial_miss_lick_timecourse);
                        else
                            this_mean_mean_per_trial_miss_lick_timecourse=this_mean_per_trial_miss_lick_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_miss_lick_timecourse(ii,1:length(this_mean_mean_per_trial_miss_lick_timecourse))=this_mean_mean_per_trial_miss_lick_timecourse;
                    else
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).hit_included(ii)=0;
                    end

                    %S- prediction
                    this_prediction_mean_per_trial_sm_timecourse = handles_outdl.this_prediction_mean_per_trial_sm_timecourse;
                    this_mean_prediction_mean_per_trial_sm_timecourse=mean(this_prediction_mean_per_trial_sm_timecourse);
                    handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sm_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_sm_timecourse))=this_mean_prediction_mean_per_trial_sm_timecourse;

                    %CR
                    these_sm_crs=handles_outdl.st_cr_per_trial((handles_outdl.st_cr_per_trial==1)|(handles_outdl.st_fa_per_trial==1));
                    if  sum(these_sm_crs)>0
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).cr_included(ii)=1;

                        %Prediction
                        this_prediction_mean_per_trial_cr_timecourse = handles_outdl.this_prediction_mean_per_trial_sm_timecourse(logical(these_sm_crs),:);
                        if size(this_prediction_mean_per_trial_cr_timecourse,1)~=1
                            this_mean_prediction_mean_per_trial_cr_timecourse=mean(this_prediction_mean_per_trial_cr_timecourse);
                        else
                            this_mean_prediction_mean_per_trial_cr_timecourse=this_prediction_mean_per_trial_cr_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_cr_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_cr_timecourse))=this_mean_prediction_mean_per_trial_cr_timecourse;

                        %licks
                        this_mean_per_trial_cr_lick_timecourse = handles_outdl.this_per_trial_mean_sm_licks_post(logical(these_sm_crs),:);
                        if size(this_mean_per_trial_cr_lick_timecourse,1)~=1
                            this_mean_mean_per_trial_cr_lick_timecourse=mean(this_mean_per_trial_cr_lick_timecourse);
                        else
                            this_mean_mean_per_trial_cr_lick_timecourse=this_mean_per_trial_cr_lick_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_cr_lick_timecourse(ii,1:length(this_mean_mean_per_trial_cr_lick_timecourse))=this_mean_mean_per_trial_cr_lick_timecourse;
                    else
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).cr_included=0;
                    end

                    %FA
                    these_sm_fas=handles_outdl.st_fa_per_trial((handles_outdl.st_cr_per_trial==1)|(handles_outdl.st_fa_per_trial==1));
                    if  sum(these_sm_fas)>0
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).fa_included(ii)=1;

                        %Prediction
                        this_prediction_mean_per_trial_fa_timecourse = handles_outdl.this_prediction_mean_per_trial_sm_timecourse(logical(these_sm_fas),:);
                        if size(this_prediction_mean_per_trial_fa_timecourse,1)~=1
                            this_mean_prediction_mean_per_trial_fa_timecourse=mean(this_prediction_mean_per_trial_fa_timecourse);
                        else
                            this_mean_prediction_mean_per_trial_fa_timecourse=this_prediction_mean_per_trial_fa_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_fa_timecourse(ii,1:length(this_mean_prediction_mean_per_trial_fa_timecourse))=this_mean_prediction_mean_per_trial_fa_timecourse;

                        %licks
                        this_mean_per_trial_fa_lick_timecourse = handles_outdl.this_per_trial_mean_sm_licks_post(logical(these_sm_fas),:);
                        if size(this_mean_per_trial_fa_lick_timecourse,1)~=1
                            this_mean_mean_per_trial_fa_lick_timecourse=mean(this_mean_per_trial_fa_lick_timecourse);
                        else
                            this_mean_mean_per_trial_fa_lick_timecourse=this_mean_per_trial_fa_lick_timecourse;
                        end
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_fa_lick_timecourse(ii,1:length(this_mean_mean_per_trial_fa_lick_timecourse))=this_mean_mean_per_trial_fa_lick_timecourse;
                    else
                        handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).fa_included(ii)=0;
                    end

                    %                 per_trial_scores_sp=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_scores_sp;
                    %                 mean_per_trial_scores_sp=zeros(size(per_trial_scores_sp,2),size(per_trial_scores_sp,3));
                    %                 mean_per_trial_scores_sp(:,:)=mean(per_trial_scores_sp,1);
                    %                 if ~isempty(per_trial_scores_sp)
                    %                     handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii,1:2,1:size(per_trial_scores_sp,3))=mean_per_trial_scores_sp;
                    %                 end
                    %
                    %                 per_trial_scores_sm=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_scores_sm;
                    %                 mean_per_trial_scores_sm=zeros(size(per_trial_scores_sm,2),size(per_trial_scores_sm,3));
                    %                 mean_per_trial_scores_sm(:,:)=mean(per_trial_scores_sm,1);
                    %                 handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii,1:2,1:size(per_trial_scores_sm,3))=mean_per_trial_scores_sm;

                    handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mouseNo(ii)=handles_out2.mouseNo_per_file(fileNo);

                end
                fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
            end
        else
            fprintf(1, ['Import for file No ' num2str(fileNo) ' failed\n'])
            if fileNo==69
                pffft=1;
            end
        end
        pffft=1;
    end
else
    load([choiceBatchPathName choiceFileName])
end

% Bar graph plot for accuracy
figureNo=0;


%Plot decoding accuracy for each group per mouse
per_mouse_acc=[];
these_mice=handles_out2.ii_thr(ii_thr).MLalgo(iiMLalgo).mouseNo;
for mouseNo=1:length(handles_out2.mouse_names)

    per_mouse_acc.mouse(mouseNo).mean_correct_predict=[];
    per_mouse_acc.mouse(mouseNo).mean_correct_predict_sh=[];


    if show_per_mouse==1
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    end



    %         handles_out2.ii_thr(ii_thr).MLalgo(iiMLalgo).mouseNo


    these_per_corr=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict,2));
    these_per_corr(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_correct_predict=mean(these_per_corr);
    per_mouse_acc.mouse(mouseNo).correct_predict=these_per_corr;


    if show_per_mouse==1
        if size(these_per_corr,1)>0
            if size(these_per_corr,1)>2

                CIpv = bootci(1000, @mean, these_per_corr);
                meanpv=mean(these_per_corr,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;


                [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'k');
            else
                if size(these_per_corr,1)>1
                    plot(time_span',mean(these_per_corr,1)', 'k');
                else
                    plot(time_span',these_per_corr', 'k');
                end

            end
        end
    end


    these_per_corr=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh,2));
    these_per_corr(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_correct_predict_sh=mean(these_per_corr);
    per_mouse_acc.mouse(mouseNo).correct_predict_sh=these_per_corr;


    if show_per_mouse==1
        if size(these_per_corr,1)>0
            if size(these_per_corr,1)>2

                CIpv = bootci(1000, @mean, these_per_corr);
                meanpv=mean(these_per_corr,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;


                [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'r');
            else
                if size(these_per_corr,1)>1
                    plot(time_span',mean(these_per_corr,1)', 'r');
                else
                    plot(time_span',these_per_corr', 'r');
                end

            end

            this_ylim=ylim;
            plot([0 0],this_ylim,'-k')

        end
    end





    if show_per_mouse==1


        ylim([0.3 1.2])
        this_ylim=ylim;
        plot([0 0],this_ylim,'-k')


        xlim([-5 15])

        xlabel('Time(sec)')




        title(['Decoding accuracy for mouse no ' num2str(mouseNo)])

    end

    %Now calculate per mouse prediction
    these_sp_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sp_timecourse,2));
    these_sp_prediction(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sp_timecourse(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_sp_prediction=mean(these_sp_prediction);
    per_mouse_acc.mouse(mouseNo).sp_prediction=these_sp_prediction;

    these_sm_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sm_timecourse,2));
    these_sm_prediction(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sm_timecourse(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_sm_prediction=mean(these_sm_prediction);
    per_mouse_acc.mouse(mouseNo).sm_prediction=these_sm_prediction;

    %Error prediciton

    %Hit
    these_hit_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_hit_timecourse,2));
    these_hit_prediction(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_hit_timecourse(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_hit_prediction=mean(these_hit_prediction);
    per_mouse_acc.mouse(mouseNo).hit_prediction=these_hit_prediction;

    %Miss
    these_miss_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_miss_timecourse,2));
    ii_included=0;
    these_ii=find(these_mice==mouseNo);
    for ii=these_ii
        try
            this_miss_prediction=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_miss_timecourse(ii,:);
            found_it=1;
        catch
            found_it=0;
        end
        if found_it==1
            ii_included=ii_included+1;
            these_miss_prediction(ii_included,:)=this_miss_prediction;
        end
    end
    these_miss_prediction=these_miss_prediction(1:ii_included,:);
    per_mouse_acc.mouse(mouseNo).mean_miss_prediction=mean(these_miss_prediction);
    per_mouse_acc.mouse(mouseNo).miss_prediction=these_miss_prediction;

    these_cr_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_cr_timecourse,2));
    these_cr_prediction(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_cr_timecourse(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_cr_prediction=mean(these_cr_prediction);
    per_mouse_acc.mouse(mouseNo).cr_prediction=these_cr_prediction;

    these_fa_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_fa_timecourse,2));
    ii_included=0;
    these_ii=find(these_mice==mouseNo);
    for ii=these_ii
        try
            this_fa_prediction=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_fa_timecourse(ii,:);
            found_it=1;
        catch
            found_it=0;
        end
        if found_it==1
            ii_included=ii_included+1;
            these_fa_prediction(ii_included,:)=this_fa_prediction;
        end
    end
    these_fa_prediction=these_fa_prediction(1:ii_included,:);
    per_mouse_acc.mouse(mouseNo).mean_fa_prediction=mean(these_fa_prediction);
    per_mouse_acc.mouse(mouseNo).fa_prediction=these_fa_prediction;

    %Licks

    %Hit
    these_hit_lick=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_hit_lick_timecourse,2));
    these_hit_lick(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_hit_lick_timecourse(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_hit_lick=mean(these_hit_lick);
    per_mouse_acc.mouse(mouseNo).hit_lick=these_hit_lick;

    %Miss
    these_miss_lick=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_miss_lick_timecourse,2));
    ii_included=0;
    these_ii=find(these_mice==mouseNo);
    for ii=these_ii
        try
            this_miss_lick=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_miss_lick_timecourse(ii,:);
            found_it=1;
        catch
            found_it=0;
        end
        if found_it==1
            ii_included=ii_included+1;
            these_miss_lick(ii_included,:)=this_miss_lick;
        end
    end
    these_miss_lick=these_miss_lick(1:ii_included,:);
    per_mouse_acc.mouse(mouseNo).mean_miss_lick=mean(these_miss_lick);
    per_mouse_acc.mouse(mouseNo).miss_lick=these_miss_lick;

    these_cr_lick=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_cr_lick_timecourse,2));
    these_cr_lick(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_cr_lick_timecourse(these_mice==mouseNo,:);
    per_mouse_acc.mouse(mouseNo).mean_cr_lick=mean(these_cr_lick);
    per_mouse_acc.mouse(mouseNo).cr_lick=these_cr_lick;

     %fa
    these_fa_lick=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_fa_lick_timecourse,2));
    ii_included=0;
    these_ii=find(these_mice==mouseNo);
    for ii=these_ii
        try
            this_fa_lick=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_fa_lick_timecourse(ii,:);
            found_it=1;
        catch
            found_it=0;
        end
        if found_it==1
            ii_included=ii_included+1;
            these_fa_lick(ii_included,:)=this_fa_lick;
        end
    end
    these_fa_lick=these_fa_lick(1:ii_included,:);
    per_mouse_acc.mouse(mouseNo).mean_fa_lick=mean(these_fa_lick);
    per_mouse_acc.mouse(mouseNo).fa_lick=these_fa_lick;

end

%Plot the overall accuracy calculated from per mouse decoding

these_mean_per_mouse_accuracy=[];
these_mean_per_mouse_accuracy_sh=[];

for mouseNo=1:length(handles_out2.mouse_names)
    these_mean_per_mouse_accuracy(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_correct_predict;
    these_mean_per_mouse_accuracy_sh(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_correct_predict_sh;
end


if size(these_mean_per_mouse_accuracy,1)>2


    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    CIpv = bootci(1000, @mean, these_mean_per_mouse_accuracy_sh);
    meanpv=mean(these_mean_per_mouse_accuracy_sh,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;

    [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_accuracy_sh,1)', CIpv', 'r');

    CIpv = bootci(1000, @mean, these_mean_per_mouse_accuracy);
    meanpv=mean(these_mean_per_mouse_accuracy,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;


    [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_accuracy,1)', CIpv', 'k');

    for ii_mouse=1:size(these_mean_per_mouse_accuracy,1)
        plot(time_span',smoothdata(these_mean_per_mouse_accuracy(ii_mouse,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
    end

    plot(time_span',mean(these_mean_per_mouse_accuracy,1)', 'k','LineWidth',1.5);

    ylim([0.3 1.2])
    this_ylim=ylim;

    %FV
    rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])


    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 delta_odor],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
    plot([delta_odor delta_odor],this_ylim,'-k')

    %Reinforcement markers
    plot([delta_odor_on_reinf_on delta_odor_on_reinf_on],this_ylim,'-r')
    reinfhl=plot([delta_odor_on_reinf_on delta_odor_on_reinf_on+delta_reinf],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
    plot([delta_odor_on_reinf_on+delta_reinf delta_odor_on_reinf_on+delta_reinf],this_ylim,'-r')


    xlim([-5 15])

    title(['Decoding accuracy calculated per mouse for ' handles_out2.classifier_names{iiMLalgo}])
    pffft=1;
end

%Now do a bar graph
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

glm_acc=[];
glm_acc_ii=0;

glm_acc_no_sh=[];
glm_acc_no_sh_ii=0;


id_acc_ii=0;
input_acc_data=[];


bar_offset=0;



edges=[0:0.05:1];
rand_offset=0.8;



handles_out2.all_accs=[];

for window_no=1:size(time_periods_eu,1)

    t_from=time_periods_eu(window_no,1);
    t_to=time_periods_eu(window_no,2);

    per_mouse_mean_accuracy=zeros(1,length(handles_out2.mouse_names));
    per_mouse_mean_accuracy_sh=zeros(1,length(handles_out2.mouse_names));
    all_session_accuracy=[];
    all_session_accuracy_sh=[];

    %Get per session and per mouse
    ii_m_included=0;
    for mouseNo=1:length(handles_out2.mouse_names)
        these_accs=[];
        these_accs_sh=[];

        these_accs=mean(per_mouse_acc.mouse(mouseNo).correct_predict(:,(time_span>=t_from)&(time_span<=t_to)),2);
        these_accs_sh=mean(per_mouse_acc.mouse(mouseNo).correct_predict_sh(:,(time_span>=t_from)&(time_span<=t_to)),2);


        per_mouse_mean_accuracy(1,mouseNo)=mean(these_accs);
        per_mouse_mean_accuracy_sh(1,mouseNo)=mean(these_accs_sh);


        all_session_accuracy=[all_session_accuracy these_accs'];
        all_session_accuracy_sh=[all_session_accuracy_sh these_accs_sh'];

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=these_accs_sh;
        %             glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=perCorr_no*ones(1,length(these_accs_sh));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=window_no*ones(1,length(these_accs_sh));
        glm_acc.shuffled(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=ones(1,length(these_accs_sh));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=mouseNo*ones(1,length(these_accs_sh));
        glm_acc_ii=glm_acc_ii+length(these_accs_sh);

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        %             glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=perCorr_no*ones(1,length(these_accs));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=window_no*ones(1,length(these_accs));
        glm_acc.shuffled(glm_acc_ii+1:glm_acc_ii+length(these_accs))=zeros(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouseNo*ones(1,length(these_accs));
        glm_acc_ii=glm_acc_ii+length(these_accs);

        glm_acc_no_sh.data(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=these_accs;
        %             glm_acc_no_sh.pcorr(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=perCorr_no*ones(1,length(these_accs));
        glm_acc_no_sh.window(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=window_no*ones(1,length(these_accs));
        glm_acc_no_sh.shuffled(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=ones(1,length(these_accs));
        glm_acc_no_sh.mouse_nos(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=mouseNo*ones(1,length(these_accs));
        glm_acc_no_sh_ii=glm_acc_no_sh_ii+length(these_accs);


    end

    id_acc_ii=id_acc_ii+1;
    input_acc_data(id_acc_ii).data=all_session_accuracy_sh;
    input_acc_data(id_acc_ii).description=['Shuffled '  window_label{window_no}];


    id_acc_ii=id_acc_ii+1;
    input_acc_data(id_acc_ii).data=all_session_accuracy;
    input_acc_data(id_acc_ii).description=[ window_label{window_no}];


    %Accuracy shuffled
    bar(bar_offset,mean(all_session_accuracy_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[200/255 200/255 200/255])

    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_session_accuracy_sh...
        ,edges,bar_offset,rand_offset,'k','k',2);



    bar_offset=bar_offset+1;

    %Accuracy
    bar(bar_offset,mean(all_session_accuracy),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])

    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_session_accuracy...
        ,edges,bar_offset,rand_offset,'k','k',2);



    for ii_mouse=1:length(handles_out2.mouse_names)
        plot([bar_offset-1 bar_offset],[per_mouse_mean_accuracy_sh(ii_mouse) per_mouse_mean_accuracy(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+1;


end

%     xticks([1 3 5 7 9 11 13 15])
%     xticklabels({'Base','PreFV','PreOdor','Odor','Base','PreFV','PreOdor','Odor'})
xticks([0.5 2.5 4.5])
xticklabels({'Pre','Odor','Reinf'})


title(['Lick decoding accuracy for ' handles_out2.classifier_names{iiMLalgo}])
ylabel('Accuracy')
ylim([0.4 1.1])
xlim([-1 6])


%Perform the glm not including shuffled
fprintf(1, ['\nglm for lick decoding accuracy\n'])
fprintf(fileID, ['\nglm for lick decoding accuracy\n']);

%         fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
%         fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_acc.data',glm_acc.shuffled',glm_acc.window',glm_acc.mouse_nos',...
    'VariableNames',{'accuracy','shuffled','window','mice'});
mdl = fitglm(tbl,'accuracy~window+shuffled+mice+shuffled*window'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for lick decoding accuracy\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for lick decoding accuracy\n']);


[output_data] = drgMutiRanksumorTtest(input_acc_data, fileID,0);


fprintf(fileID, '\n\n');


pffft=1;

%Now show prediction graphs
%Plot the overall accuracy calculated from per mouse decoding

% per_mouse_acc.mouse(mouseNo).mean_sp_prediction=mean(these_sp_prediction);
% per_mouse_acc.mouse(mouseNo).sp_prediction=these_sp_prediction;
% 
% these_sm_prediction=zeros(sum(these_mice==mouseNo),size(handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sm_timecourse,2));
% these_sm_prediction(:,:)=handles_out2.ii_thr(ii_out).MLalgo(iiMLalgo).mean_prediction_mean_per_trial_sm_timecourse(these_mice==mouseNo,:);
% per_mouse_acc.mouse(mouseNo).mean_sm_prediction=mean(these_sm_prediction);
% per_mouse_acc.mouse(mouseNo).sm_prediction=these_sm_prediction;

these_mean_per_mouse_sp_prediciton=[];
these_mean_per_mouse_sm_prediction=[];

for mouseNo=1:length(handles_out2.mouse_names)
    these_mean_per_mouse_sp_prediciton(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_sp_prediction;
    these_mean_per_mouse_sm_prediciton(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_sm_prediction;
end


if size(these_mean_per_mouse_sp_prediciton,1)>2


    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    CIpv = bootci(1000, @mean, these_mean_per_mouse_sm_prediciton);
    meanpv=mean(these_mean_per_mouse_sm_prediciton,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;

    [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_sm_prediciton,1)', CIpv', 'cmap',[158/255 31/255 99/255]);

    CIpv = bootci(1000, @mean, these_mean_per_mouse_sp_prediciton);
    meanpv=mean(these_mean_per_mouse_sp_prediciton,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;


    [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_sp_prediciton,1)', CIpv', 'cmap',[0 114/255 178/255]);

    %     for ii_mouse=1:size(these_mean_per_mouse_accuracy,1)
    %         plot(time_span',smoothdata(these_mean_per_mouse_accuracy(ii_mouse,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
    %     end

    %     plot(time_span',mean(these_mean_per_mouse_accuracy,1)', 'k','LineWidth',1.5);
    ylim([-0.1 1.1])
    this_ylim=ylim;
    xlim([-5 15])

    %FV
    rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 delta_odor],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
    plot([delta_odor delta_odor],this_ylim,'-k')

    %Reinforcement markers
    plot([delta_odor_on_reinf_on delta_odor_on_reinf_on],this_ylim,'-r')
    reinfhl=plot([delta_odor_on_reinf_on delta_odor_on_reinf_on+delta_reinf],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
    plot([delta_odor_on_reinf_on+delta_reinf delta_odor_on_reinf_on+delta_reinf],this_ylim,'-r')
    

    title(['Lick prediction calculated per mouse for ' handles_out2.classifier_names{iiMLalgo}])
    pffft=1;
end

%Now do a bar graph
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

glm_acc=[];
glm_acc_ii=0;

glm_acc_no_sh=[];
glm_acc_no_sh_ii=0;


id_acc_ii=0;
input_acc_data=[];


bar_offset=0;



edges=[0:0.05:1];
rand_offset=0.8;



handles_out2.all_accs=[];

for window_no=1:size(time_periods_eu,1)

    t_from=time_periods_eu(window_no,1);
    t_to=time_periods_eu(window_no,2);

    per_mouse_mean_accuracy=zeros(1,length(handles_out2.mouse_names));
    per_mouse_mean_accuracy_sh=zeros(1,length(handles_out2.mouse_names));
    all_session_accuracy=[];
    all_session_accuracy_sh=[];

    %Get per session and per mouse
    ii_m_included=0;
    for mouseNo=1:length(handles_out2.mouse_names)
        these_accs=[];
        these_accs_sh=[];

        these_accs=mean(per_mouse_acc.mouse(mouseNo).sp_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        these_accs_sh=mean(per_mouse_acc.mouse(mouseNo).sm_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);


        per_mouse_mean_accuracy(1,mouseNo)=mean(these_accs);
        per_mouse_mean_accuracy_sh(1,mouseNo)=mean(these_accs_sh);


        all_session_accuracy=[all_session_accuracy these_accs'];
        all_session_accuracy_sh=[all_session_accuracy_sh these_accs_sh'];

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=these_accs_sh;
        %             glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=perCorr_no*ones(1,length(these_accs_sh));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=window_no*ones(1,length(these_accs_sh));
        glm_acc.shuffled(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=ones(1,length(these_accs_sh));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=mouseNo*ones(1,length(these_accs_sh));
        glm_acc_ii=glm_acc_ii+length(these_accs_sh);

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        %             glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=perCorr_no*ones(1,length(these_accs));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=window_no*ones(1,length(these_accs));
        glm_acc.shuffled(glm_acc_ii+1:glm_acc_ii+length(these_accs))=zeros(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouseNo*ones(1,length(these_accs));
        glm_acc_ii=glm_acc_ii+length(these_accs);

        glm_acc_no_sh.data(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=these_accs;
        %             glm_acc_no_sh.pcorr(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=perCorr_no*ones(1,length(these_accs));
        glm_acc_no_sh.window(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=window_no*ones(1,length(these_accs));
        glm_acc_no_sh.shuffled(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=ones(1,length(these_accs));
        glm_acc_no_sh.mouse_nos(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=mouseNo*ones(1,length(these_accs));
        glm_acc_no_sh_ii=glm_acc_no_sh_ii+length(these_accs);


    end

    id_acc_ii=id_acc_ii+1;
    input_acc_data(id_acc_ii).data=all_session_accuracy_sh;
    input_acc_data(id_acc_ii).description=['S+ '  window_label{window_no}];


    id_acc_ii=id_acc_ii+1;
    input_acc_data(id_acc_ii).data=all_session_accuracy;
    input_acc_data(id_acc_ii).description=['S- ' window_label{window_no}];


    %Accuracy shuffled
    bar(bar_offset,mean(all_session_accuracy_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[200/255 200/255 200/255])

    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_session_accuracy_sh...
        ,edges,bar_offset,rand_offset,'k','k',2);



    bar_offset=bar_offset+1;

    %Accuracy
    bar(bar_offset,mean(all_session_accuracy),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])

    %Violin plot
    [mean_out, CIout]=drgViolinPoint(all_session_accuracy...
        ,edges,bar_offset,rand_offset,'k','k',2);



    for ii_mouse=1:length(handles_out2.mouse_names)
        plot([bar_offset-1 bar_offset],[per_mouse_mean_accuracy_sh(ii_mouse) per_mouse_mean_accuracy(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+1;


end

%     xticks([1 3 5 7 9 11 13 15])
%     xticklabels({'Base','PreFV','PreOdor','Odor','Base','PreFV','PreOdor','Odor'})
xticks([0.5 2.5 4.5])
xticklabels({'Pre','Odor','Reinf'})


title(['Lick prediction for ' handles_out2.classifier_names{iiMLalgo}])
ylabel('Prediction')
ylim([0 1])
xlim([-1 6])


%Perform the glm not including shuffled
fprintf(1, ['\nglm for lick prediction\n'])
fprintf(fileID, ['\nglm for lick prediction\n']);

%         fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
%         fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_acc.data',glm_acc.shuffled',glm_acc.window',glm_acc.mouse_nos',...
    'VariableNames',{'accuracy','spm','window','mice'});
mdl = fitglm(tbl,'accuracy~window+spm+mice+spm*window'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for lick prediction\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for lick prediction\n']);


[output_data] = drgMutiRanksumorTtest(input_acc_data, fileID,0);


fprintf(fileID, '\n\n');

pfft=1;

%Now show prediction graph for errors

%Now show prediction graphs for errors

these_mean_per_mouse_hit_prediciton=[];
these_mean_per_mouse_miss_prediction=[];
these_mean_per_mouse_cr_prediciton=[];
these_mean_per_mouse_fa_prediction=[];

per_mouse_licks_sp_hit=[];
per_mouse_licks_sp_miss=[];
per_mouse_licks_sm_cr=[];
per_mouse_licks_sm_fa=[];

for mouseNo=1:length(handles_out2.mouse_names)
    these_mean_per_mouse_hit_prediciton(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_hit_prediction;
    these_mean_per_mouse_miss_prediciton(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_miss_prediction;
    these_mean_per_mouse_cr_prediciton(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_cr_prediction;
    these_mean_per_mouse_fa_prediciton(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_fa_prediction;

    per_mouse_licks_sp_hit(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_hit_lick;
    per_mouse_licks_sp_miss(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_miss_lick;
    per_mouse_licks_sm_cr(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_cr_lick;
    per_mouse_licks_sm_fa(mouseNo,:)=per_mouse_acc.mouse(mouseNo).mean_fa_lick;
end



figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

CIpv = bootci(1000, @mean, these_mean_per_mouse_miss_prediciton);
meanpv=mean(these_mean_per_mouse_miss_prediciton,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;

[hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_miss_prediciton,1)', CIpv', 'c');

CIpv = bootci(1000, @mean, these_mean_per_mouse_fa_prediciton);
meanpv=mean(these_mean_per_mouse_fa_prediciton,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;


[hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_fa_prediciton,1)', CIpv', 'm');

CIpv = bootci(1000, @mean, these_mean_per_mouse_cr_prediciton);
meanpv=mean(these_mean_per_mouse_cr_prediciton,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;

[hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_cr_prediciton,1)', CIpv', 'b');

CIpv = bootci(1000, @mean, these_mean_per_mouse_hit_prediciton);
meanpv=mean(these_mean_per_mouse_hit_prediciton,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;


[hlpvl, hppvl] = boundedline(time_span',mean(these_mean_per_mouse_hit_prediciton,1)', CIpv', 'r');

plot(time_span',mean(these_mean_per_mouse_miss_prediciton,1)', 'c');
plot(time_span',mean(these_mean_per_mouse_fa_prediciton,1)', 'm');

%     for ii_mouse=1:size(these_mean_per_mouse_accuracy,1)
%         plot(time_span',smoothdata(these_mean_per_mouse_accuracy(ii_mouse,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
%     end

%     plot(time_span',mean(these_mean_per_mouse_accuracy,1)', 'k','LineWidth',1.5);
ylim([-0.1 1.1])
this_ylim=ylim;
xlim([-5 15])

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor on markers
plot([0 0],this_ylim,'-k')
odorhl=plot([0 delta_odor],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
plot([delta_odor delta_odor],this_ylim,'-k')

%Reinforcement markers
plot([delta_odor_on_reinf_on delta_odor_on_reinf_on],this_ylim,'-r')
reinfhl=plot([delta_odor_on_reinf_on delta_odor_on_reinf_on+delta_reinf],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
plot([delta_odor_on_reinf_on+delta_reinf delta_odor_on_reinf_on+delta_reinf],this_ylim,'-r')


title(['Lick prediction calculated per mouse for ' handles_out2.classifier_names{iiMLalgo}])
pffft=1;
 

%Bar graph for glm prediction
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

glm_pred=[];
glm_pred_ii=0;


id_pred_ii=0;
input_pred_data=[];

% for grNo=1:no_pcorr*length(these_groups)
bar_offset=0;


per_session=[];

mean_pred=[];
for window_no=1:size(time_periods_eu,1)

    per_mouse_mean_pred_sp_hit=[];
    per_mouse_mean_pred_sp_miss=[];
    per_mouse_mean_pred_sm_cr=[];
    per_mouse_mean_pred_sm_fa=[];
    all_session_pred_sp_hit=[];
    all_session_pred_sp_miss=[];
    all_session_pred_sm_cr=[];
    all_session_pred_sm_fa=[];

     t_from=time_periods_eu(window_no,1);
                    t_to=time_periods_eu(window_no,2);

    %Get per session and per mouse

    %fa
    for mouseNo=1:length(handles_out2.mouse_names)
        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).fa_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).fa_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_pred_sm_fa(1,mouseNo)=mean(these_preds);
        all_session_pred_sm_fa=[all_session_pred_sm_fa these_preds'];

        glm_pred.data(glm_pred_ii+1:glm_pred_ii+length(these_preds))=these_preds;
        glm_pred.window(glm_pred_ii+1:glm_pred_ii+length(these_preds))=window_no*ones(1,length(these_preds));
        glm_pred.trial_type(glm_pred_ii+1:glm_pred_ii+length(these_preds))=0*ones(1,length(these_preds));
        glm_pred.mouse_nos(glm_pred_ii+1:glm_pred_ii+length(these_preds))=mouseNo*ones(1,length(these_preds));
        glm_pred_ii=glm_pred_ii+length(these_preds);
    end

    id_pred_ii=id_pred_ii+1;
    input_pred_data(id_pred_ii).data=all_session_pred_sm_fa;
    input_pred_data(id_pred_ii).description=['FA ' window_label{window_no}];


    %FA
    bar(bar_offset,mean(all_session_pred_sm_fa),'LineWidth', 3,'EdgeColor','none','FaceColor','m')

    %Violin plot
    if ~isempty(all_session_pred_sm_fa)
        [mean_out, CIout]=drgViolinPoint(all_session_pred_sm_fa...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end
    
    mean_pred.window(window_no).fa_mean_pred=mean_out;
    mean_pred.window(window_no).fa_CI_pred=CIout;

    per_session.window(window_no).all_session_pred_sm_fa=all_session_pred_sm_fa;

    bar_offset=bar_offset+1;

    %cr
    for mouseNo=1:length(handles_out2.mouse_names)
        these_preds=[];
        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).cr_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).cr_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_pred_sm_cr(1,mouseNo)=mean(these_preds);
        all_session_pred_sm_cr=[all_session_pred_sm_cr these_preds'];

        glm_pred.data(glm_pred_ii+1:glm_pred_ii+length(these_preds))=these_preds;
        glm_pred.window(glm_pred_ii+1:glm_pred_ii+length(these_preds))=window_no*ones(1,length(these_preds));
        glm_pred.trial_type(glm_pred_ii+1:glm_pred_ii+length(these_preds))=1*ones(1,length(these_preds));
        glm_pred.mouse_nos(glm_pred_ii+1:glm_pred_ii+length(these_preds))=mouseNo*ones(1,length(these_preds));
        glm_pred_ii=glm_pred_ii+length(these_preds);
    end

    id_pred_ii=id_pred_ii+1;
    input_pred_data(id_pred_ii).data=all_session_pred_sm_cr;
    input_pred_data(id_pred_ii).description=['CR ' window_label{window_no}];


    %CR
    bar(bar_offset,mean(all_session_pred_sm_cr),'LineWidth', 3,'EdgeColor','none','FaceColor','b')

    %Violin plot
    if ~isempty(all_session_pred_sm_cr)
        [mean_out, CIout]=drgViolinPoint(all_session_pred_sm_cr...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).cr_mean_pred=mean_out;
    mean_pred.window(window_no).cr_CI_pred=CIout;

    per_session.window(window_no).all_session_pred_sm_cr=all_session_pred_sm_cr;

    bar_offset=bar_offset+1;

    %miss
    for mouseNo=1:length(handles_out2.mouse_names)
        these_preds=[];
        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).miss_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).miss_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_pred_sp_miss(1,mouseNo)=mean(these_preds);
        all_session_pred_sp_miss=[all_session_pred_sp_miss these_preds'];

        glm_pred.data(glm_pred_ii+1:glm_pred_ii+length(these_preds))=these_preds;
        glm_pred.window(glm_pred_ii+1:glm_pred_ii+length(these_preds))=window_no*ones(1,length(these_preds));
        glm_pred.trial_type(glm_pred_ii+1:glm_pred_ii+length(these_preds))=2*ones(1,length(these_preds));
        glm_pred.mouse_nos(glm_pred_ii+1:glm_pred_ii+length(these_preds))=mouseNo*ones(1,length(these_preds));
        glm_pred_ii=glm_pred_ii+length(these_preds);
    end

    id_pred_ii=id_pred_ii+1;
    input_pred_data(id_pred_ii).data=all_session_pred_sp_miss;
    input_pred_data(id_pred_ii).description=['Miss ' window_label{window_no}];


    %Miss
    bar(bar_offset,mean(all_session_pred_sp_miss),'LineWidth', 3,'EdgeColor','none','FaceColor','c')

    %Violin plot
    if ~isempty(all_session_pred_sp_miss)
        [mean_out, CIout]=drgViolinPoint(all_session_pred_sp_miss...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).miss_mean_pred=mean_out;
    mean_pred.window(window_no).miss_CI_pred=CIout;

    per_session.window(window_no).all_session_pred_sp_miss=all_session_pred_sp_miss;

    bar_offset=bar_offset+1;

    %hit

    for mouseNo=1:length(handles_out2.mouse_names)
        these_preds=[];
        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).hit_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).hit_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_pred_sp_hit(1,mouseNo)=mean(these_preds');
        all_session_pred_sp_hit=[all_session_pred_sp_hit these_preds'];

        glm_pred.data(glm_pred_ii+1:glm_pred_ii+length(these_preds))=these_preds;
        glm_pred.window(glm_pred_ii+1:glm_pred_ii+length(these_preds))=window_no*ones(1,length(these_preds));
        glm_pred.trial_type(glm_pred_ii+1:glm_pred_ii+length(these_preds))=3*ones(1,length(these_preds));
        glm_pred.mouse_nos(glm_pred_ii+1:glm_pred_ii+length(these_preds))=mouseNo*ones(1,length(these_preds));
        glm_pred_ii=glm_pred_ii+length(these_preds);
    end

    id_pred_ii=id_pred_ii+1;
    input_pred_data(id_pred_ii).data=all_session_pred_sp_hit;
    input_pred_data(id_pred_ii).description=['Hit ' window_label{window_no}];


    %Hit
    bar(bar_offset,mean(all_session_pred_sp_hit),'LineWidth', 3,'EdgeColor','none','FaceColor','r')

    %Violin plot
    if ~isempty(all_session_pred_sp_hit)
        [mean_out, CIout]=drgViolinPoint(all_session_pred_sp_hit...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).hit_mean_pred=mean_out;
    mean_pred.window(window_no).hit_CI_pred=CIout;

    per_session.window(window_no).all_session_pred_sp_hit=all_session_pred_sp_hit;

    bar_offset=bar_offset+1;

    for ii_mouse=1:4
        plot([bar_offset-4 bar_offset-3],[per_mouse_mean_pred_sm_fa(ii_mouse) per_mouse_mean_pred_sm_cr(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
        plot([bar_offset-3 bar_offset-2],[per_mouse_mean_pred_sm_cr(ii_mouse) per_mouse_mean_pred_sp_miss(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
        plot([bar_offset-2 bar_offset-1],[per_mouse_mean_pred_sp_miss(ii_mouse) per_mouse_mean_pred_sp_hit(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+1;


end


xticks([0 1 2 3 5 6 7 8])
xticklabels({'FA RA1','CR RA1','Miss RA1','Hit RA1','FA RA2','CR RA2','Miss RA2','Hit RA2'})


title(['Prediction for proficient, hit=red, miss=cyan, fa=magenta, cr=blue'])
ylabel('Prediction, S+=1, S-=0')
ylim([0 1.1])
xlim([-1 9])


%Perform the glm  for errors
fprintf(1, ['\nglm for lick prediction\n'])
fprintf(fileID, ['\nglm for lick prediction\n']);
% 
% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_pred.data',glm_pred.window',glm_pred.trial_type',glm_pred.mouse_nos',...
    'VariableNames',{'prediction','window','trial_type','mice'});
mdl = fitglm(tbl,'prediction~window+trial_type+mice+trial_type*window'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for lick prediction\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for lick prediction\n']);
 

[output_data] = drgMutiRanksumorTtest(input_pred_data, fileID,0);

%Now do the lick fraction plots
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

%Plot hits
if size(per_mouse_licks_sp_hit,1)>=3
    CIpv = bootci(1000, @mean, per_mouse_licks_sp_hit);
    meanpv=mean(per_mouse_licks_sp_hit,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;

    [hlpvl, hppvl] = boundedline(time_span',mean(per_mouse_licks_sp_hit,1)', CIpv', '-r');
else
    plot(time_span',mean(per_mouse_licks_sp_hit,1)', '-r');
end

%Plot miss
if size(per_mouse_licks_sp_miss,1)>=3
    CIpv = bootci(1000, @mean, per_mouse_licks_sp_miss);
    meanpv=mean(per_mouse_licks_sp_miss,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;

    [hlpvl, hppvl] = boundedline(time_span',mean(per_mouse_licks_sp_miss,1)', CIpv', '-c');
else
    plot(time_span',mean(per_mouse_licks_sp_miss,1)', '-c');
end

%Plot cr
if size(per_mouse_licks_sm_cr,1)>=3
    CIpv = bootci(1000, @mean, per_mouse_licks_sm_cr);
    meanpv=mean(per_mouse_licks_sm_cr,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;

    [hlpvl, hppvl] = boundedline(time_span',mean(per_mouse_licks_sm_cr,1)', CIpv', '-b');
else
    plot(time_span',mean(per_mouse_licks_sm_cr,1)', '-b');
end

%Plot fa
if size(per_mouse_licks_sm_fa,1)>=3
    CIpv = bootci(1000, @mean, per_mouse_licks_sm_fa);
    meanpv=mean(per_mouse_licks_sm_fa,1);
    CIpv(1,:)=meanpv-CIpv(1,:);
    CIpv(2,:)=CIpv(2,:)-meanpv;

    [hlpvl, hppvl] = boundedline(time_span',mean(per_mouse_licks_sm_fa,1)', CIpv', '-m');
else
    plot(time_span',mean(per_mouse_licks_sm_fa,1)', '-m');
end

plot(time_span',mean(per_mouse_licks_sp_miss,1)', '-c');
plot(time_span',mean(per_mouse_licks_sp_hit,1)', '-r');
plot(time_span',mean(per_mouse_licks_sm_fa,1)', '-m');
plot(time_span',mean(per_mouse_licks_sm_cr,1)', '-b');

ylim([-0.1 1.1])
this_ylim=ylim;


xlim([-7 15])

this_grNo=floor((grNo-1)/4)+1;
ii_pcorr=grNo-4*(this_grNo-1);

%FV
plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'LineWidth',5, Color=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor on markers
plot([0 0],this_ylim,'-k')
odorhl=plot([0 delta_odor],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
plot([delta_odor delta_odor],this_ylim,'-k')

%Reinforcement markers
plot([delta_odor_on_reinf_on delta_odor_on_reinf_on],this_ylim,'-r')
reinfhl=plot([delta_odor_on_reinf_on delta_odor_on_reinf_on+delta_reinf],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
plot([delta_odor_on_reinf_on+delta_reinf delta_odor_on_reinf_on+delta_reinf],this_ylim,'-r')


title(['Lick fraction for proficient, hit=red, miss=cyan, fa=magenta, cr=blue'])
xlabel('Time(sec)')
ylabel('Lick fraction')



% end 

%Bar graph for glm prediction
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

glm_lick=[];
glm_lick_ii=0;


id_lick_ii=0;
input_lick_data=[];

% for grNo=1:no_pcorr*length(these_groups)
bar_offset=0;

RA_lick_stats=[];


for window_no=1:size(time_periods_eu,1)
    
    RA_lick_stats.window(window_no).ii=0;

    per_mouse_mean_lick_sp_hit=[];
    per_mouse_mean_lick_sp_miss=[];
    per_mouse_mean_lick_sm_cr=[];
    per_mouse_mean_lick_sm_fa=[];
    all_session_lick_sp_hit=[];
    all_session_lick_sp_miss=[];
    all_session_lick_sm_cr=[];
    all_session_lick_sm_fa=[];

    t_from=time_periods_eu(window_no,1);
    t_to=time_periods_eu(window_no,2);

    %Get per session and per mouse

    %fa

    for mouseNo=1:length(handles_out2.mouse_names)
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).fa_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).fa_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_lick_sm_fa(1,mouseNo)=mean(these_licks');
        all_session_lick_sm_fa=[all_session_lick_sm_fa these_licks'];

        RA_lick_stats.window(window_no).lick_fracs(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=these_licks;
        RA_lick_stats.window(window_no).lick_types(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=1*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).mouse_sm(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).ii=RA_lick_stats.window(window_no).ii+length(these_licks);

        glm_lick.data(glm_lick_ii+1:glm_lick_ii+length(these_licks))=these_licks;
        glm_lick.window(glm_lick_ii+1:glm_lick_ii+length(these_licks))=window_no*ones(1,length(these_licks));
        glm_lick.trial_type(glm_lick_ii+1:glm_lick_ii+length(these_licks))=0*ones(1,length(these_licks));
        glm_lick.mouse_nos(glm_lick_ii+1:glm_lick_ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        glm_lick_ii=glm_lick_ii+length(these_licks);

    end

    id_lick_ii=id_lick_ii+1;
    input_lick_data(id_lick_ii).data=all_session_lick_sm_fa;
    input_lick_data(id_lick_ii).description=['FA ' window_label{window_no}];


    %FA
    bar(bar_offset,mean(all_session_lick_sm_fa),'LineWidth', 3,'EdgeColor','none','FaceColor','m')

    %Violin plot
    if ~isempty(all_session_lick_sm_fa)
        [mean_out, CIout]=drgViolinPoint(all_session_lick_sm_fa...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).fa_mean_lick_frac=mean_out;
    mean_pred.window(window_no).fa_CI_lick_frac=CIout;

    per_session.window(window_no).all_session_lick_sm_fa=all_session_lick_sm_fa;
    bar_offset=bar_offset+1;

    %cr
    for mouseNo=1:length(handles_out2.mouse_names)
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).cr_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).cr_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_lick_sm_cr(1,mouseNo)=mean(these_licks');
        all_session_lick_sm_cr=[all_session_lick_sm_cr these_licks'];

        RA_lick_stats.window(window_no).lick_fracs(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=these_licks;
        RA_lick_stats.window(window_no).lick_types(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=2*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).mouse_sm(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).ii=RA_lick_stats.window(window_no).ii+length(these_licks);

        glm_lick.data(glm_lick_ii+1:glm_lick_ii+length(these_licks))=these_licks;
        glm_lick.window(glm_lick_ii+1:glm_lick_ii+length(these_licks))=window_no*ones(1,length(these_licks));
        glm_lick.trial_type(glm_lick_ii+1:glm_lick_ii+length(these_licks))=1*ones(1,length(these_licks));
        glm_lick.mouse_nos(glm_lick_ii+1:glm_lick_ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        glm_lick_ii=glm_lick_ii+length(these_licks);


    end

    id_lick_ii=id_lick_ii+1;
    input_lick_data(id_lick_ii).data=all_session_lick_sm_cr;
    input_lick_data(id_lick_ii).description=['CR ' window_label{window_no}];


    %CR
    bar(bar_offset,mean(all_session_lick_sm_cr),'LineWidth', 3,'EdgeColor','none','FaceColor','b')

    %Violin plot
    if ~isempty(all_session_lick_sm_cr)
        [mean_out, CIout]=drgViolinPoint(all_session_lick_sm_cr...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).cr_mean_lick_frac=mean_out;
    mean_pred.window(window_no).cr_CI_lick_frac=CIout;

    per_session.window(window_no).all_session_lick_sm_cr=all_session_lick_sm_cr;
    bar_offset=bar_offset+1;

    %miss
    for mouseNo=1:length(handles_out2.mouse_names)
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).miss_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).miss_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_lick_sp_miss(1,mouseNo)=mean(these_licks');
        all_session_lick_sp_miss=[all_session_lick_sp_miss these_licks'];

        RA_lick_stats.window(window_no).lick_fracs(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=these_licks;
        RA_lick_stats.window(window_no).lick_types(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=3*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).mouse_sm(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).ii=RA_lick_stats.window(window_no).ii+length(these_licks);

        glm_lick.data(glm_lick_ii+1:glm_lick_ii+length(these_licks))=these_licks;
        glm_lick.window(glm_lick_ii+1:glm_lick_ii+length(these_licks))=window_no*ones(1,length(these_licks));
        glm_lick.trial_type(glm_lick_ii+1:glm_lick_ii+length(these_licks))=2*ones(1,length(these_licks));
        glm_lick.mouse_nos(glm_lick_ii+1:glm_lick_ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        glm_lick_ii=glm_lick_ii+length(these_licks);

    end

    id_lick_ii=id_lick_ii+1;
    input_lick_data(id_lick_ii).data=all_session_lick_sp_miss;
    input_lick_data(id_lick_ii).description=['Miss ' window_label{window_no}];


    %Miss
    bar(bar_offset,mean(all_session_lick_sp_miss),'LineWidth', 3,'EdgeColor','none','FaceColor','c')

    %Violin plot
    if ~isempty(all_session_lick_sp_miss)
        [mean_out, CIout]=drgViolinPoint(all_session_lick_sp_miss...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).miss_mean_lick_frac=mean_out;
    mean_pred.window(window_no).miss_CI_lick_frac=CIout;

    per_session.window(window_no).all_session_lick_sp_miss=all_session_lick_sp_miss;
    bar_offset=bar_offset+1;

    %hit
    for mouseNo=1:length(handles_out2.mouse_names)

        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).hit_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).hit_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        per_mouse_mean_lick_sp_hit(1,mouseNo)=mean(these_licks');
        all_session_lick_sp_hit=[all_session_lick_sp_hit these_licks'];


        RA_lick_stats.window(window_no).lick_fracs(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=these_licks;
        RA_lick_stats.window(window_no).lick_types(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=4*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).mouse_sm(RA_lick_stats.window(window_no).ii+1:RA_lick_stats.window(window_no).ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        RA_lick_stats.window(window_no).ii=RA_lick_stats.window(window_no).ii+length(these_licks);

        glm_lick.data(glm_lick_ii+1:glm_lick_ii+length(these_licks))=these_licks;
        glm_lick.window(glm_lick_ii+1:glm_lick_ii+length(these_licks))=window_no*ones(1,length(these_licks));
        glm_lick.trial_type(glm_lick_ii+1:glm_lick_ii+length(these_licks))=3*ones(1,length(these_licks));
        glm_lick.mouse_nos(glm_lick_ii+1:glm_lick_ii+length(these_licks))=mouseNo*ones(1,length(these_licks));
        glm_lick_ii=glm_lick_ii+length(these_licks);


    end

    id_lick_ii=id_lick_ii+1;
    input_lick_data(id_lick_ii).data=all_session_lick_sp_hit;
    input_lick_data(id_lick_ii).description=['Hit ' window_label{window_no}];


    %Hit
    bar(bar_offset,mean(all_session_lick_sp_hit),'LineWidth', 3,'EdgeColor','none','FaceColor','r')

    %Violin plot
    if ~isempty(all_session_lick_sp_hit)
        [mean_out, CIout]=drgViolinPoint(all_session_lick_sp_hit...
            ,edges,bar_offset,rand_offset,'k','k',2);
    end

    mean_pred.window(window_no).hit_mean_lick_frac=mean_out;
    mean_pred.window(window_no).hit_CI_lick_frac=CIout;
    per_session.window(window_no).all_session_lick_sp_hit=all_session_lick_sp_hit;
    bar_offset=bar_offset+1;

    for ii_mouse=1:4
        plot([bar_offset-4 bar_offset-3],[per_mouse_mean_lick_sm_fa(ii_mouse) per_mouse_mean_lick_sm_cr(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
        plot([bar_offset-3 bar_offset-2],[per_mouse_mean_lick_sm_cr(ii_mouse) per_mouse_mean_lick_sp_miss(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
        plot([bar_offset-2 bar_offset-1],[per_mouse_mean_lick_sp_miss(ii_mouse) per_mouse_mean_lick_sp_hit(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
    end

    bar_offset=bar_offset+1;


end
bar_offset=bar_offset+1;


% xticks([1 3 5 7 9 11 13 15])
% xticklabels({'Base','PreFV','PreOdor','Odor','Base','PreFV','PreOdor','Odor'})

xticks([0 1 2 3 5 6 7 8])
xticklabels({'FA RA1','CR RA1','Miss RA1','Hit RA1','FA RA2','CR RA2','Miss RA2','Hit RA2'})


title(['Lick fraction proficient, hit=red, miss=cyan, fa=magenta, cr=blue'])
ylabel('Lick fraction')
ylim([0 0.8])
xlim([-1 9])

%     text(1,1,'45-65%')
%     text(6,1,'65-80%')
%     text(11,1,'>80%')

 

%Perform the glm for licks
fprintf(1, ['\nglm for decoding lick fraction\n'])
fprintf(fileID, ['\nglm for decoding lick fraction\n']);

% fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
% fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

tbl = table(glm_lick.data',glm_lick.window',glm_lick.trial_type',glm_lick.mouse_nos',...
    'VariableNames',{'lick_f','window','trial_type','mice'});
mdl = fitglm(tbl,'lick_f~window+trial_type+mice+trial_type*window'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for decoding lick fraction\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for decoding lick fraction\n']);
 

[output_data] = drgMutiRanksumorTtest(input_lick_data, fileID,0);


%Now plot the relationship between prediction and lick fraction
% window_labels{1}='Pre';
window_labels{1}='OdorRA1';
window_labels{2}='OdorRA2';
for window_no=1:2
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    hold on
    plot(per_session.window(window_no).all_session_lick_sm_cr,per_session.window(window_no).all_session_pred_sm_cr,'ob')
    plot(per_session.window(window_no).all_session_lick_sm_fa,per_session.window(window_no).all_session_pred_sm_fa,'om')
    plot(per_session.window(window_no).all_session_lick_sp_hit,per_session.window(window_no).all_session_pred_sp_hit,'or')
    plot(per_session.window(window_no).all_session_lick_sp_miss,per_session.window(window_no).all_session_pred_sp_miss,'oc')

    xlim([-0.05 0.7])
    ylim([0 1])
    xlabel('lick fraction')
    ylabel('prediction')
    title(['Prediction vs. lick fraction ' window_labels{window_no}])
end


%Plot the overall lick vs prediction
%Now plot mean predictions for both windows
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

hold on
spmx=zeros(1,length(mean_pred.window)*2);
spmy=zeros(1,length(mean_pred.window)*2);
spm_ii=0;



for window_no=1:length(mean_pred.window)
    plot([mean_pred.window(window_no).cr_CI_pred(1) mean_pred.window(window_no).cr_CI_pred(2)],[mean_pred.window(window_no).cr_mean_lick_frac mean_pred.window(window_no).cr_mean_lick_frac],'-k','LineWidth',2)
    plot([mean_pred.window(window_no).cr_mean_pred mean_pred.window(window_no).cr_mean_pred],[mean_pred.window(window_no).cr_CI_lick_frac(1) mean_pred.window(window_no).cr_CI_lick_frac(2)],'-k','LineWidth',2)
    plot(mean_pred.window(window_no).cr_mean_pred,mean_pred.window(window_no).cr_mean_lick_frac,'ob','MarkerFaceColor','b','MarkerSize',8)
    spm_ii=spm_ii+1;
    spm_y(spm_ii)=mean_pred.window(window_no).cr_mean_lick_frac;
    spm_x(spm_ii)=mean_pred.window(window_no).cr_mean_pred;

    plot([mean_pred.window(window_no).fa_CI_pred(1) mean_pred.window(window_no).fa_CI_pred(2)],[mean_pred.window(window_no).fa_mean_lick_frac mean_pred.window(window_no).fa_mean_lick_frac],'-k','LineWidth',2)
    plot([mean_pred.window(window_no).fa_mean_pred mean_pred.window(window_no).fa_mean_pred],[mean_pred.window(window_no).fa_CI_lick_frac(1) mean_pred.window(window_no).fa_CI_lick_frac(2)],'-k','LineWidth',2)
    plot(mean_pred.window(window_no).fa_mean_pred,mean_pred.window(window_no).fa_mean_lick_frac,'om','MarkerFaceColor','m','MarkerSize',8)
    spm_ii=spm_ii+1;
    spm_y(spm_ii)=mean_pred.window(window_no).fa_mean_lick_frac;
    spm_x(spm_ii)=mean_pred.window(window_no).fa_mean_pred;

    plot([mean_pred.window(window_no).hit_CI_pred(1) mean_pred.window(window_no).hit_CI_pred(2)],[mean_pred.window(window_no).hit_mean_lick_frac mean_pred.window(window_no).hit_mean_lick_frac],'-k','LineWidth',2)
    plot([mean_pred.window(window_no).hit_mean_pred mean_pred.window(window_no).hit_mean_pred],[mean_pred.window(window_no).hit_CI_lick_frac(1) mean_pred.window(window_no).hit_CI_lick_frac(2)],'-k','LineWidth',2)
    plot(mean_pred.window(window_no).hit_mean_pred,mean_pred.window(window_no).hit_mean_lick_frac,'or','MarkerFaceColor','r','MarkerSize',8)
    spm_ii=spm_ii+1;
    spm_y(spm_ii)=mean_pred.window(window_no).hit_mean_lick_frac;
    spm_x(spm_ii)=mean_pred.window(window_no).hit_mean_pred;

    plot([mean_pred.window(window_no).miss_CI_pred(1) mean_pred.window(window_no).miss_CI_pred(2)],[mean_pred.window(window_no).miss_mean_lick_frac mean_pred.window(window_no).miss_mean_lick_frac],'-k','LineWidth',2)
    plot([mean_pred.window(window_no).miss_mean_pred mean_pred.window(window_no).miss_mean_pred],[mean_pred.window(window_no).miss_CI_lick_frac(1) mean_pred.window(window_no).miss_CI_lick_frac(2)],'-k','LineWidth',2)
    plot(mean_pred.window(window_no).miss_mean_pred,mean_pred.window(window_no).miss_mean_lick_frac,'oc','MarkerFaceColor','c','MarkerSize',8)
    spm_ii=spm_ii+1;
    spm_y(spm_ii)=mean_pred.window(window_no).miss_mean_lick_frac;
    spm_x(spm_ii)=mean_pred.window(window_no).miss_mean_pred;
end

p=polyfit(spm_x,spm_y,1);
x1=[0:0.05:0.4];
f1 = polyval(p,x1);
plot(x1,f1,'-k','LineWidth',2)


ylim([-0.05 0.5])
xlim([0 0.4])
ylabel('lick fraction')
xlabel('prediction')
title(['Lick fraction vs. predcition '])

%How about the licks/predictions for every 0.5 sec period
dt_sample=0.5;

this_t=-5;

all_pred_sm_fa=[];
all_pred_sm_cr=[];
all_pred_sp_hit=[];
all_pred_sp_miss=[];

all_lick_sm_fa=[];
all_lick_sm_cr=[];
all_lick_sp_hit=[];
all_lick_sp_miss=[];



while this_t<=7-dt_sample

    t_from=this_t;
    t_to=this_t+dt_sample;

    %hit
    for mouseNo=1:length(handles_out2.mouse_names)

        %fa
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).fa_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).fa_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_lick_sm_fa=[all_lick_sm_fa mean(these_licks)'];

        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).fa_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).fa_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_pred_sm_fa=[all_pred_sm_fa mean(these_preds)'];

        %cr
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).cr_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).cr_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_lick_sm_cr=[all_lick_sm_cr mean(these_licks)'];

        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).cr_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).cr_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_pred_sm_cr=[all_pred_sm_cr mean(these_preds)'];

        %miss
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).miss_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).miss_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_lick_sp_miss=[all_lick_sp_miss mean(these_licks)'];

        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).miss_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).miss_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_pred_sp_miss=[all_pred_sp_miss mean(these_preds)'];

        %hit
        these_licks=[];
        these_licks=zeros(1,size(per_mouse_acc.mouse(mouseNo).hit_lick,1));
        these_licks=mean(per_mouse_acc.mouse(mouseNo).hit_lick(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_lick_sp_hit=[all_lick_sp_hit mean(these_licks)'];

        these_preds=zeros(1,size(per_mouse_acc.mouse(mouseNo).hit_prediction,1));
        these_preds=mean(per_mouse_acc.mouse(mouseNo).hit_prediction(:,(time_span>=t_from)&(time_span<=t_to)),2);
        all_pred_sp_hit=[all_pred_sp_hit mean(these_preds)'];

        

    end

    this_t=this_t+dt_sample;

end

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

hold on
% plot(all_pred_sp_hit,all_lick_sp_hit,'or')
% plot(all_pred_sm_cr,all_lick_sm_cr,'ob')
% plot(all_pred_sm_fa,all_lick_sm_fa,'om')
% plot(all_pred_sp_miss,all_lick_sp_miss,'oc')

plot(all_lick_sp_hit,all_pred_sp_hit,'or')
plot(all_lick_sm_cr,all_pred_sm_cr,'ob')
plot(all_lick_sm_fa,all_pred_sm_fa,'om')
plot(all_lick_sp_miss,all_pred_sp_miss,'oc')


spm_y=[all_pred_sp_hit];
spm_x=[all_lick_sp_hit];


p=polyfit(spm_x,spm_y,1);
x1=[0:0.05:0.6];
f1 = polyval(p,x1);
plot(x1,f1,'-k','LineWidth',2)

ylim([0 0.6])
xlim([0 0.7])
xlabel('lick fraction')
ylabel('prediction')
title(['Lick prediction vs. lick fraction ' ])

fclose(fileID);

pfft=1;