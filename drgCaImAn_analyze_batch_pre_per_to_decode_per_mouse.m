%drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse
close all
clear all

 
[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdzv2 run for ' choiceFileName '\n\n']);


tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

fileID = fopen([choiceBatchPathName 'drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse.txt'],'w');

no_files=handles.no_files;
moving_mean_n=10;



%Which threshold value should we use?
ii_thr=4;

%Show the per mouse graphs?
show_per_mouse=0;

is_Fabio=0;

switch is_Fabio
    case 0
        
        %Choices for Ming's go-no go processing
        no_pcorr=4;
        
        %groups to be shown in the zoomed figures for Ming's data
        grNo1=4; %Forward proficient
        grNo1_label='forward proficient';
        grNo2=8; %Forward proficient
        grNo2_label='reversed proficient';
        
    case 1
        
        %Chaoices for Fabio
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

 
these_groups=unique(handles.group);

if exist([choiceBatchPathName choiceFileName(1:end-2) '.mat'])==0
    %Process each file separately
    for grNo=1:no_pcorr*length(these_groups)
        handles_out2.group_no(grNo).ii_euclid=0;
        for iiMLalgo=handles.MLalgo_to_use
            for ii_out=1:length(handles.p_threshold)
                handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=0;
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

            grNo=(this_grNo-1)*no_pcorr+ii_pCorr;

%             %Is this the first session?
%             if sum(first_last(:,1)==fileNo)==1
%                 mouseNo=find(first_last(:,1)==fileNo);
%                 for iiMLalgo=handles.MLalgo_to_use
%                     for ii_out=1:length(handles_out.ii_out)
%                         this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
%                         handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict(1,1:length(this_mean_correct_predict))=this_mean_correct_predict;
% 
%                         this_mean_correct_predict_sh=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh);
%                         handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict_sh(1,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;
%                     end
%                 end
%             end
% 
%             %Is this the last session?
%             if sum(first_last(:,2)==fileNo)
%                 mouseNo=find(first_last(:,2)==fileNo);
%                 for iiMLalgo=handles.MLalgo_to_use
%                     for ii_out=1:length(handles_out.ii_out)
%                         this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
%                         handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict(1,1:length(this_mean_correct_predict))=this_mean_correct_predict;
% 
%                         this_mean_correct_predict_sh=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh);
%                         handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict_sh(1,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;
%                     end
%                 end
%             end

            for iiMLalgo=handles.MLalgo_to_use
                if iiMLalgo==handles.MLalgo_to_use(1)
                    handles_out2.group_no(grNo).ii_euclid=handles_out2.group_no(grNo).ii_euclid+1;
                    ii_euclid=handles_out2.group_no(grNo).ii_euclid;
                    handles_out2.group_no(grNo).dist_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.dist_euclid))=handles_out.ii_out(1).handles_out.dist_euclid-handles_out.ii_out(1).handles_out.dist_euclid_zero;
                    handles_out2.group_no(grNo).KLdivergence(ii_euclid,1:length(handles_out.ii_out(1).handles_out.KLdivergence))=handles_out.ii_out(1).handles_out.KLdivergence;
                    handles_out2.group_no(grNo).time_span_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.time_span;
                    handles_out2.group_no(grNo).ii_time_span(ii_euclid,1)=length(handles_out.ii_out(1).handles_out.time_span);
                    handles_out2.group_no(grNo).meandFFsp(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsp;
                    handles_out2.group_no(grNo).meandFFsm(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsm;
                    handles_out2.group_no(grNo).mouseNo_euclid(ii_euclid)= handles_out2.mouseNo_per_file(fileNo);
                end
                for ii_out=1:length(handles_out.ii_out)
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii+1;
                    ii=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii;
                    accuracy_tr=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr;
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr(ii)=accuracy_tr;

                    accuracy_tr_sh=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh;
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_sh(ii)=accuracy_tr_sh;

                    accuracy_tr_sh2=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).accuracy_tr_sh2;
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_sh2(ii)=accuracy_tr_sh2;

                    this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(ii,1:length(this_mean_correct_predict))=this_mean_correct_predict;

                    this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(ii,1:length(this_mean_correct_predict))=this_mean_correct_predict;

                    this_mean_correct_predict_sh=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh(ii,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;

                    this_per_trial_sp_timecourse=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_sp_timecourse;
                    this_per_trial_sm_timecourse=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_sm_timecourse;

                    this_moving_mean_per_trial_sp_timecourse = movmean(this_per_trial_sp_timecourse',moving_mean_n)';
                    this_mean_moving_mean_per_trial_sp_timecourse=mean(this_moving_mean_per_trial_sp_timecourse);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii,1:length(this_mean_moving_mean_per_trial_sp_timecourse))=this_mean_moving_mean_per_trial_sp_timecourse;

                    this_moving_mean_per_trial_sm_timecourse = movmean(this_per_trial_sm_timecourse',moving_mean_n)';
                    this_mean_moving_mean_per_trial_sm_timecourse=mean(this_moving_mean_per_trial_sm_timecourse);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii,1:length(this_mean_moving_mean_per_trial_sm_timecourse))=this_mean_moving_mean_per_trial_sm_timecourse;

                    per_trial_scores_sp=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_scores_sp;
                    mean_per_trial_scores_sp=zeros(size(per_trial_scores_sp,2),size(per_trial_scores_sp,3));
                    mean_per_trial_scores_sp(:,:)=mean(per_trial_scores_sp);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii,1:2,1:size(per_trial_scores_sp,3))=mean_per_trial_scores_sp;

                    per_trial_scores_sm=handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).per_trial_scores_sm;
                    mean_per_trial_scores_sm=zeros(size(per_trial_scores_sm,2),size(per_trial_scores_sm,3));
                    mean_per_trial_scores_sm(:,:)=mean(per_trial_scores_sm);
                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii,1:2,1:size(per_trial_scores_sm,3))=mean_per_trial_scores_sm;

                    handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mouseNo(ii)=handles_out2.mouseNo_per_file(fileNo);
                end
            end
            fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
        else
            fprintf(1, ['Import for file No ' num2str(fileNo) ' failed\n'])
            if fileNo==69
                pffft=1;
            end
        end
    end
else
    load([choiceBatchPathName choiceFileName(1:end-2) '.mat'])
end

% %  ii_thr=length(handles.p_threshold);
%  ii_thr=1;

%Bar graph plot for accuracy
figureNo=0;
for iiMLalgo=handles.MLalgo_to_use
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .5 .25])
    
    hold on
    
    edges=[0:0.05:1];
    rand_offset=0.8;
    
    bar_offset=0;
    
    for grNo=1:no_pcorr*length(these_groups)
        
       
        
        if handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).ii>0
            %Shuffled
            bar_offset=bar_offset+1;
            these_accuracy_tr_sh2=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).accuracy_tr_sh2;
            bar(bar_offset,mean(these_accuracy_tr_sh2),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            
            if length(these_accuracy_tr_sh2)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(these_accuracy_tr_sh2...
                    ,edges,bar_offset,rand_offset,'k','k',3);
            end
            
            
            %Accuracy
            bar_offset=bar_offset+1;
            these_accuracy_tr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).accuracy_tr;
            bar(bar_offset,mean(these_accuracy_tr),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            
            if length(these_accuracy_tr)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(these_accuracy_tr...
                    ,edges,bar_offset,rand_offset,'k','k',3);
            end
        else
            bar_offset=bar_offset+2;
        end
        
        bar_offset=bar_offset+1;
    end
    xticks([1.5:3:22.5])
    labels='xticklabels({';
    for ii_label=1:length(these_groups)
        for ii_pcorr=1:no_pcorr
            labels=[labels '''' handles.group_names{these_groups(ii_label)} per_names{ii_pcorr} ''', '];
        end
    end
    labels=[labels(1:end-2) '})'];
    eval(labels)
    title(['Prediction accuracy for ' handles_out2.classifier_names{iiMLalgo}])
    ylabel('Accuracy')
    ylim([0 1])
    xlim([0 24])
    pffft=1;
    
end

%Plot the timecourse of mean dFF

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

for grNo=1:no_pcorr*length(these_groups)
     
    subplot(8,1,grNo)
    hold on
         
    [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
    if ~isempty(ii_tspan)

        %Extrapolate all points onto the longest ii_tspan
        time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);

        these_meandFFsp=zeros(size(handles_out2.group_no(grNo).ii_time_span,1),ii_tspan);
        these_meandFFsm=handles_out2.group_no(grNo).meandFFsm(1,1:ii_tspan);
        for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
            if ii_f==ii_file
                these_meandFFsp=handles_out2.group_no(grNo).meandFFsp(ii_file,1:ii_tspan);
                these_meandFFsm=handles_out2.group_no(grNo).meandFFsm(ii_file,1:ii_tspan);
            else
                this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
                this_meandFFsp=zeros(1,ii_tspan);
                for ii_tsp=1:ii_tspan
                    if time_span(ii_tsp)<this_time_span(1)
                        these_meandFFsp(ii_f,ii_tsp)=handles_out2.group_no(grNo).meandFFsp(ii_f,1);
                        these_meandFFsm(ii_f,ii_tsp)=handles_out2.group_no(grNo).meandFFsm(ii_f,1);
                    else
                        if time_span(ii_tsp)>this_time_span(end)
                            these_meandFFsp(ii_f,ii_tsp)=handles_out2.group_no(grNo).meandFFsp(ii_f,end);
                            these_meandFFsm(ii_f,ii_tsp)=handles_out2.group_no(grNo).meandFFsm(ii_f,end);
                        else
                            ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                            ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                            these_meandFFsp(ii_f,ii_tsp)=handles_out2.group_no(grNo).meandFFsp(ii_f,ii_0)+...
                                (handles_out2.group_no(grNo).meandFFsp(ii_f,ii_1)-handles_out2.group_no(grNo).meandFFsp(ii_f,ii_0))...
                                *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                             these_meandFFsm(ii_f,ii_tsp)=handles_out2.group_no(grNo).meandFFsm(ii_f,ii_0)+...
                                (handles_out2.group_no(grNo).meandFFsm(ii_f,ii_1)-handles_out2.group_no(grNo).meandFFsm(ii_f,ii_0))...
                                *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                        end
                    end
                end

            end
                
        end
        
        
         

        if size(these_meandFFsm,1)>2
            
            CIpv = bootci(1000, @mean, these_meandFFsm);
            meanpv=mean(these_meandFFsm,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
             
            
            [hlpvl, hppvl] = boundedline(time_span',mean(these_meandFFsm,1)', CIpv', 'b');
        else
            if size(these_meandFFsm,1)>0
                plot(time_span',mean(these_meandFFsm,1)', 'b');
            end
            
        end
        
        
        if size(these_meandFFsp,1)>2
            
            CIpv = bootci(1000, @mean, these_meandFFsp);
            meanpv=mean(these_meandFFsp,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
            
            
            [hlpvl, hppvl] = boundedline(time_span',mean(these_meandFFsp,1)', CIpv', 'r');
        else
            if size(these_meandFFsp,1)>0
                plot(time_span',mean(these_meandFFsp,1)', 'r');
            end
            
        end
        
    end
    xlim([-10 20])
    xlabel('Time(sec)')

    if no_pcorr==1
        title([handles.group_names{grNo}])
    else

        this_grNo=floor((grNo-1)/4)+1;
        ii_pcorr=grNo-4*(this_grNo-1);
        title([handles.group_names{these_groups(this_grNo)} ' ' per_names{ii_pcorr}])
    end
    pffft=1;
end
sgtitle('Mean dFF')
  
%Now plot the timecourse for euclidean distance
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

for grNo=1:no_pcorr*length(these_groups)
    
    [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
    subplot(8,1,grNo)
    hold on
    
    if ~isempty(ii_tspan)
        time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
        
       
        
        %Extrapolate all points onto the longest ii_tspan
        time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);

        these_dist_euclid=zeros(size(handles_out2.group_no(grNo).ii_time_span,1),ii_tspan);
        
        for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
            if ii_f==ii_file
                these_dist_euclid=handles_out2.group_no(grNo).dist_euclid(ii_file,1:ii_tspan);
            else
                this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
                this_dist_euclid=zeros(1,ii_tspan);
                for ii_tsp=1:ii_tspan
                    if time_span(ii_tsp)<this_time_span(1)
                        these_dist_euclid(ii_f,ii_tsp)=handles_out2.group_no(grNo).dist_euclid(ii_f,1);
                        
                    else
                        if time_span(ii_tsp)>this_time_span(end)
                            these_dist_euclid(ii_f,ii_tsp)=handles_out2.group_no(grNo).dist_euclid(ii_f,end);
                        else
                            ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                            ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                            these_dist_euclid(ii_f,ii_tsp)=handles_out2.group_no(grNo).dist_euclid(ii_f,ii_0)+...
                                (handles_out2.group_no(grNo).dist_euclid(ii_f,ii_1)-handles_out2.group_no(grNo).dist_euclid(ii_f,ii_0))...
                                *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                        end
                    end
                end

            end
                
        end
        
        if size(these_dist_euclid,1)>2
            
            CIpv = bootci(1000, @mean, these_dist_euclid);
            meanpv=mean(these_dist_euclid,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
            
            
            [hlpvl, hppvl] = boundedline(time_span',mean(these_dist_euclid,1)', CIpv', 'm');
        else
            if size(these_dist_euclid,1)>0
                plot(time_span',mean(these_dist_euclid,1)', 'm');
            end
            
        end
        this_ylim=ylim;
        plot([0 0],this_ylim,'-k')
    end
    xlim([-10 20])
    xlabel('Time(sec)')
    if no_pcorr==1
        title([handles.group_names{grNo}])
    else
        this_grNo=floor((grNo-1)/4)+1;
        ii_pcorr=grNo-4*(this_grNo-1);
        title([handles.group_names{these_groups(this_grNo)} ' ' per_names{ii_pcorr}])
    end
end
sgtitle('Euclidean distance')

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

delta_KLdiv_post=[];
for grNo=1:no_pcorr*length(these_groups)
    
    [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
    subplot(8,1,grNo)
    hold on
        
    if ~isempty(ii_tspan)
        time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);

        %Extrapolate all points onto the longest ii_tspan
        time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
        
        for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
            if ii_f==ii_file
                these_KLdivergence=handles_out2.group_no(grNo).KLdivergence(ii_file,1:ii_tspan);
            else
                this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
                this_KLdivergence=zeros(1,ii_tspan);
                for ii_tsp=1:ii_tspan
                    if time_span(ii_tsp)<this_time_span(1)
                        these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,1);
                        
                    else
                        if time_span(ii_tsp)>this_time_span(end)
                            these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,end);
                        else
                            ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                            ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                            these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0)+...
                                (handles_out2.group_no(grNo).KLdivergence(ii_f,ii_1)-handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0))...
                                *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                        end
                    end
                end

            end
                
        end
        
        %Calculate the delta_KLdiv
        delta_KLdiv=[];
        these_delta_KLdiv_post=[];
        for ii_file=1:ii_f
            this_delta_KLdiv=[];
            this_delta_KLdiv=these_KLdivergence(ii_file,:)-mean(these_KLdivergence(ii_file,(time_span'>-20)&(time_span'<-2)));
            delta_KLdiv(ii_file,:)=this_delta_KLdiv;
            these_delta_KLdiv_post=[these_delta_KLdiv_post mean(this_delta_KLdiv((time_span'>=0)&(time_span'<=handles.post_time)))];
        end
        delta_KLdiv_post.group(grNo).delta_KLdiv_post=these_delta_KLdiv_post;

        if size(these_KLdivergence,1)>2
            
            CIpv = bootci(1000, @mean, delta_KLdiv);
            meanpv=mean(delta_KLdiv,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
            
            delta_mean_KLdiv=mean(delta_KLdiv,1)';
            
            [hlpvl, hppvl] = boundedline(time_span',delta_mean_KLdiv, CIpv', 'm');
        else
            if size(these_KLdivergence,1)>0
                delta_mean_KLdiv=mean(delta_KLdiv,1)';
                plot(time_span',delta_mean_KLdiv, 'm');
            end
            
        end
        this_ylim=ylim;
        plot([0 0],this_ylim,'-k')
    end
    
        xlim([-10 20])
    
    xlabel('Time(sec)')

    if no_pcorr==1
        title([handles.group_names{grNo}])
    else
        this_grNo=floor((grNo-1)/4)+1;
        ii_pcorr=grNo-4*(this_grNo-1);
        title([handles.group_names{these_groups(this_grNo)} ' ' per_names{ii_pcorr}])
    end

end
sgtitle('KL divergence')

%Plot a bar graph for post_time KL divergence
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.3 .3 .5 .25])

hold on

edges=[0:0.05:1];
rand_offset=0.8;

bar_offset=0;

for grNo=1:no_pcorr*length(these_groups)



    if length(delta_KLdiv_post.group(grNo).delta_KLdiv_post)>0
     
        bar_offset=bar_offset+1;
        these_delta_KLdiv=delta_KLdiv_post.group(grNo).delta_KLdiv_post;
        bar(bar_offset,mean(these_delta_KLdiv),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])

        if length(these_delta_KLdiv)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_delta_KLdiv...
                ,edges,bar_offset,rand_offset,'k','k',3);
        end
    else
        bar_offset=bar_offset+1;
    end

    bar_offset=bar_offset+1;
end
xticks([1:2:15])
labels='xticklabels({';
for ii_label=1:length(these_groups)
    for ii_pcorr=1:no_pcorr
        labels=[labels '''' handles.group_names{these_groups(ii_label)} per_names{ii_pcorr} ''', '];
    end
end
labels=[labels(1:end-2) '})'];
eval(labels)
title('delta KL divergence ')
ylabel('delta KL div')




%Plot the delta KL divergence for grNo1 and grNo2
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])




%Do forward
grNo=grNo1;

[ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);


if ~isempty(ii_tspan)
    time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
    
    %Extrapolate all points onto the longest ii_tspan
    time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
    
    for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
        if ii_f==ii_file
            these_KLdivergence=handles_out2.group_no(grNo).KLdivergence(ii_file,1:ii_tspan);
        else
            this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
            this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
            this_KLdivergence=zeros(1,ii_tspan);
            for ii_tsp=1:ii_tspan
                if time_span(ii_tsp)<this_time_span(1)
                    these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,1);
                    
                else
                    if time_span(ii_tsp)>this_time_span(end)
                        these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,end);
                    else
                        ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                        ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                        these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0)+...
                            (handles_out2.group_no(grNo).KLdivergence(ii_f,ii_1)-handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0))...
                            *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                    end
                end
            end
            
        end
        
    end
     
    %Calculate the delta_KLdiv
    delta_KLdiv=[];
    for ii_file=1:ii_f
        delta_KLdiv(ii_file,:)=these_KLdivergence(ii_file,:)-mean(these_KLdivergence(ii_file,(time_span'>-20)&(time_span'<-2)));
    end
    
    if size(these_KLdivergence,1)>2
        
        CIpv = bootci(1000, @mean, delta_KLdiv);
        meanpv=mean(delta_KLdiv,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
        
        delta_mean_KLdiv=mean(delta_KLdiv,1)';
        [hlpvl, hppvl] = boundedline(time_span',delta_mean_KLdiv, CIpv', 'b');
    else
        if size(these_KLdivergence,1)>0
            delta_mean_KLdiv=mean(delta_KLdiv,1)';
            plot(time_span',delta_mean_KLdiv, 'b');
        end
        
    end
    this_ylim=ylim;
        plot([0 0],this_ylim,'-k')
end

%Do reversed
grNo=grNo2;

[ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);


if ~isempty(ii_tspan)
    time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
    
    %Extrapolate all points onto the longest ii_tspan
    time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
    
    for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
        if ii_f==ii_file
            these_KLdivergence=handles_out2.group_no(grNo).KLdivergence(ii_file,1:ii_tspan);
        else
            this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
            this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
            this_KLdivergence=zeros(1,ii_tspan);
            for ii_tsp=1:ii_tspan
                if time_span(ii_tsp)<this_time_span(1)
                    these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,1);
                    
                else
                    if time_span(ii_tsp)>this_time_span(end)
                        these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,end);
                    else
                        ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                        ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                        these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0)+...
                            (handles_out2.group_no(grNo).KLdivergence(ii_f,ii_1)-handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0))...
                            *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                    end
                end
            end
            
        end
        
    end
    
      %Calculate the delta_KLdiv
    delta_KLdiv=[];
    for ii_file=1:ii_f
        delta_KLdiv(ii_file,:)=these_KLdivergence(ii_file,:)-mean(these_KLdivergence(ii_file,(time_span'>-20)&(time_span'<-2)));
    end
    
    if size(these_KLdivergence,1)>2
        
        CIpv = bootci(1000, @mean, delta_KLdiv);
        meanpv=mean(delta_KLdiv,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;
        
        
        delta_mean_KLdiv=mean(delta_KLdiv,1)';
        [hlpvl, hppvl] = boundedline(time_span',delta_mean_KLdiv, CIpv', 'r');
    else
        if size(these_KLdivergence,1)>0
            delta_mean_KLdiv=mean(delta_KLdiv,1)';
            plot(time_span',delta_mean_KLdiv, 'r');
        end
        
    end
    this_ylim=ylim;
        plot([0 0],this_ylim,'-k')
end

% ylim([-10 10])

    xlim([-10 20])

xlabel('Time(sec)')


title(['delta KL distance for blue=' grNo1_label ' red= '  grNo2_label])

% %Plot the delta KL divergence for forward proficient and naive
% figureNo = figureNo + 1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% %Do forward naive
% grNo=2;
% 
% [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
% 
% 
% if ~isempty(ii_tspan)
%     time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
%     
%     %Extrapolate all points onto the longest ii_tspan
%     time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
%     
%     for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%         if ii_f==ii_file
%             these_KLdivergence=handles_out2.group_no(grNo).KLdivergence(ii_file,1:ii_tspan);
%         else
%             this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%             this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
%             this_KLdivergence=zeros(1,ii_tspan);
%             for ii_tsp=1:ii_tspan
%                 if time_span(ii_tsp)<this_time_span(1)
%                     these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,1);
%                     
%                 else
%                     if time_span(ii_tsp)>this_time_span(end)
%                         these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,end);
%                     else
%                         ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                         ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                         these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0)+...
%                             (handles_out2.group_no(grNo).KLdivergence(ii_f,ii_1)-handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0))...
%                             *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                     end
%                 end
%             end
%             
%         end
%         
%     end
%      
%     %Calculate the delta_KLdiv
%     delta_KLdiv=[];
%     for ii_file=1:ii_f
%         delta_KLdiv(ii_file,:)=these_KLdivergence(ii_file,:)-mean(these_KLdivergence(ii_file,(time_span'>-20)&(time_span'<-2)));
%     end
%     
%     if size(these_KLdivergence,1)>2
%         
%         CIpv = bootci(1000, @mean, delta_KLdiv);
%         meanpv=mean(delta_KLdiv,1);
%         CIpv(1,:)=meanpv-CIpv(1,:);
%         CIpv(2,:)=CIpv(2,:)-meanpv;
%         
%         delta_mean_KLdiv=mean(delta_KLdiv,1)';
%         [hlpvl, hppvl] = boundedline(time_span',delta_mean_KLdiv, CIpv', 'cmap',[80/255 194/255 255/255]);
%     else
%         if size(these_KLdivergence,1)>0
%             delta_mean_KLdiv=mean(delta_KLdiv,1)';
%             plot(time_span',delta_mean_KLdiv, 'Color' ,[80/255 194/255 255/255]);
%         end
%         
%     end
%     
% end
% 
% %Do forward proficient
% grNo=4;
% 
% [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
% 
% 
% if ~isempty(ii_tspan)
%     time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
%     
%     %Extrapolate all points onto the longest ii_tspan
%     time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
%     
%     for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%         if ii_f==ii_file
%             these_KLdivergence=handles_out2.group_no(grNo).KLdivergence(ii_file,1:ii_tspan);
%         else
%             this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%             this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
%             this_KLdivergence=zeros(1,ii_tspan);
%             for ii_tsp=1:ii_tspan
%                 if time_span(ii_tsp)<this_time_span(1)
%                     these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,1);
%                     
%                 else
%                     if time_span(ii_tsp)>this_time_span(end)
%                         these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,end);
%                     else
%                         ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                         ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                         these_KLdivergence(ii_f,ii_tsp)=handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0)+...
%                             (handles_out2.group_no(grNo).KLdivergence(ii_f,ii_1)-handles_out2.group_no(grNo).KLdivergence(ii_f,ii_0))...
%                             *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                     end
%                 end
%             end
%             
%         end
%         
%     end
%     
%       %Calculate the delta_KLdiv
%     delta_KLdiv=[];
%     for ii_file=1:ii_f
%         delta_KLdiv(ii_file,:)=these_KLdivergence(ii_file,:)-mean(these_KLdivergence(ii_file,(time_span'>-20)&(time_span'<-2)));
%     end
%     
%     if size(these_KLdivergence,1)>2
%         
%         CIpv = bootci(1000, @mean, delta_KLdiv);
%         meanpv=mean(delta_KLdiv,1);
%         CIpv(1,:)=meanpv-CIpv(1,:);
%         CIpv(2,:)=CIpv(2,:)-meanpv;
%         
%         
%         delta_mean_KLdiv=mean(delta_KLdiv,1)';
%         [hlpvl, hppvl] = boundedline(time_span',delta_mean_KLdiv, CIpv', 'b');
%         
%     else
%         if size(these_KLdivergence,1)>0
%             delta_mean_KLdiv=mean(delta_KLdiv,1)';
%             plot(time_span',delta_mean_KLdiv, 'b');
%         end
%         
%     end
%     
% end
% 
% ylim([-10 10])
% xlim([-10 20])
% xlabel('Time(sec)')
% 
% 
% title('delta KL distance for naive vs. proficient forward')


% 
% 
% % %Choose threshold
% % ii_out=4;
% 
% %Decoding accuracy
% for iiMLalgo=handles.MLalgo_to_use
%     figureNo = figureNo + 1;
%     try
%         close(figureNo)
%     catch
%     end
%     hFig=figure(figureNo);
%     
%     ax=gca;ax.LineWidth=3;
%     set(hFig, 'units','normalized','position',[.2 .2 .3 .6])
%     
%     for grNo=1:no_pcorr*length(these_groups)
%         
%         [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
%         
%         subplot(8,1,grNo)
%         hold on
%         
%         if ~isempty(ii_tspan)
%             time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
%             
%            
% 
%             %Extrapolate all points onto the longest ii_tspan
%             time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
% 
%             for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%                 if ii_f==ii_file
%                     these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_file,1:ii_tspan);
%                 else
%                     this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%                     this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
%                     
%                     for ii_tsp=1:ii_tspan
%                         if time_span(ii_tsp)<this_time_span(1)
%                             these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,1);
% 
%                         else
%                             if time_span(ii_tsp)>this_time_span(end)
%                                 these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,end);
%                             else
%                                 ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                                 ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                                 these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_0)+...
%                                     (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_0))...
%                                     *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                             end
%                         end
%                     end
% 
%                 end
% 
%             end
% 
%             if size(these_per_corr,1)>2
%                 
%                 CIpv = bootci(1000, @mean, these_per_corr);
%                 meanpv=mean(these_per_corr,1);
%                 CIpv(1,:)=meanpv-CIpv(1,:);
%                 CIpv(2,:)=CIpv(2,:)-meanpv;
%                 
%                 
%                 [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'k');
%             else
%                 if size(these_per_corr,1)>1
%                     plot(time_span',mean(these_per_corr,1)', 'k');
%                 else 
%                     plot(time_span',these_per_corr', 'k');
%                 end
%                 
%             end
%             
% 
%              for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%                 if ii_f==ii_file
%                     these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_file,1:ii_tspan);
%                 else
%                     this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%                     this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
%                     
%                     for ii_tsp=1:ii_tspan
%                         if time_span(ii_tsp)<this_time_span(1)
%                             these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,1);
% 
%                         else
%                             if time_span(ii_tsp)>this_time_span(end)
%                                 these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,end);
%                             else
%                                 ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                                 ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                                 these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_0)+...
%                                     (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_0))...
%                                     *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                             end
%                         end
%                     end
% 
%                 end
% 
%              end
%       
%             if size(these_per_corr,1)>2
%                 
%                 CIpv = bootci(1000, @mean, these_per_corr);
%                 meanpv=mean(these_per_corr,1);
%                 CIpv(1,:)=meanpv-CIpv(1,:);
%                 CIpv(2,:)=CIpv(2,:)-meanpv;
%                 
%                 
%                 [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'm');
%             else
%                 if size(these_per_corr,1)>1
%                     plot(time_span',mean(these_per_corr,1)', 'm');
%                 else
%                     plot(time_span',these_per_corr', 'm');
%                 end
%                 
%             end
%             this_ylim=ylim;
%         plot([0 0],this_ylim,'-k')
%             xlim([-10 20])
%             xlabel('Time(sec)')
% 
%             if no_pcorr==1
%                 title([handles.group_names{grNo}])
%             else
%                 this_grNo=floor((grNo-1)/4)+1;
%                 ii_pcorr=grNo-4*(this_grNo-1);
%                 title([handles.group_names{these_groups(this_grNo)} ' ' per_names{ii_pcorr}])
%             end
%         end
%         
%     end
%     sgtitle(['Percent correct ' handles_out2.classifier_names{iiMLalgo}])
% end
% 
% %Decoding accuracy
% for iiMLalgo=handles.MLalgo_to_use
%     figureNo = figureNo + 1;
%     try
%         close(figureNo)
%     catch
%     end
%     hFig=figure(figureNo);
%     
%     ax=gca;ax.LineWidth=3;
%     set(hFig, 'units','normalized','position',[.2 .2 .3 .6])
%     
%     for grNo=1:length(these_groups)
%         ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
%         if ~isempty(ii_tspan)
%             time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
%             subplot(6,1,grNo)
%             hold on
%             
%             these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(:,1:ii_tspan);
%             
%             if size(these_per_corr,1)>2
%                 
%                 CIpv = bootci(1000, @mean, these_per_corr);
%                 meanpv=mean(these_per_corr,1);
%                 CIpv(1,:)=meanpv-CIpv(1,:);
%                 CIpv(2,:)=CIpv(2,:)-meanpv;
%                 
%                 
%                 [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'k');
%             else
%                 if size(these_per_corr,1)>1
%                     plot(time_span',mean(these_per_corr,1)', 'k');
%                 else 
%                     plot(time_span',these_per_corr', 'k');
%                 end
%                 
%             end
%             
%             these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(:,1:ii_tspan);
%       
%             if size(these_per_corr,1)>2
%                 
%                 CIpv = bootci(1000, @mean, these_per_corr);
%                 meanpv=mean(these_per_corr,1);
%                 CIpv(1,:)=meanpv-CIpv(1,:);
%                 CIpv(2,:)=CIpv(2,:)-meanpv;
%                 
%                 
%                 [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'm');
%             else
%                 if size(these_per_corr,1)>1
%                     plot(time_span',mean(these_per_corr,1)', 'm');
%                 else
%                     plot(time_span',these_per_corr', 'm');
%                 end
%                 
%             end
%             xlim([-10 20])
%             xlabel('Time(sec)')
%             title(handles.group_names{these_groups(grNo)})
%         end
%     end
%     sgtitle(['Percent correct ' handles_out2.classifier_names{iiMLalgo}])
% end

%Plot decoding accuracy for each group per mouse
per_mouse_acc=[];
for mouseNo=1:length(handles_out2.mouse_names)
    for grNo=1:no_pcorr*length(these_groups)

        per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy=[];
        per_mouse_acc.group(grNo).mouse(mouseNo).t_mean_accuracy=[];
        per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy_sh=[];
        per_mouse_acc.group(grNo).mouse(mouseNo).t_mean_accuracy_sh=[];

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


        [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
        %         handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mouseNo
        if ~isempty(ii_tspan)
            time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
            these_per_corr=[];


            %Extrapolate all points onto the longest ii_tspan
            time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
            these_mice=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mouseNo;
            ii_included=0;
            for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
                if these_mice(ii_f)==mouseNo
                    ii_included=ii_included+1;
                    if ii_f==ii_file
                        these_per_corr(ii_included,:)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_file,1:ii_tspan);
                    else
                        this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                        this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);

                        for ii_tsp=1:ii_tspan
                            if time_span(ii_tsp)<this_time_span(1)
                                these_per_corr(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,1);

                            else
                                if time_span(ii_tsp)>this_time_span(end)
                                    these_per_corr(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,end);
                                else
                                    ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                                    ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                                    these_per_corr(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_0)+...
                                        (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_0))...
                                        *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                end
                            end
                        end

                    end
                end
            end


            %For some reason whe have a subset of sessions with zero accuracy
            vetted_these_per_corr=[];
            ii_s=0;
            for ii_sessions=1:size(these_per_corr,1)
                if sum(these_per_corr(ii_sessions,:)==0)<200
                    ii_s=ii_s+1;
                    vetted_these_per_corr(ii_s,:)=these_per_corr(ii_sessions,:);
                end
            end

            these_per_corr=[];
            these_per_corr=vetted_these_per_corr;

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

            per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy_sh=mean(these_per_corr,1)';
            per_mouse_acc.group(grNo).mouse(mouseNo).t_mean_accuracy_sh=time_span';

            these_per_corr=[];
            ii_included=0;
            for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
                if these_mice(ii_f)==mouseNo
                    ii_included=ii_included+1;
                    if ii_f==ii_file
                        these_per_corr(ii_included,:)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_file,1:ii_tspan);
                    else
                        this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                        this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);

                        for ii_tsp=1:ii_tspan
                            if time_span(ii_tsp)<this_time_span(1)
                                these_per_corr(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,1);

                            else
                                if time_span(ii_tsp)>this_time_span(end)
                                    these_per_corr(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,end);
                                else
                                    ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                                    ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                                    these_per_corr(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_0)+...
                                        (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_0))...
                                        *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                end
                            end
                        end

                    end
                end
            end


            %For some reason whe have a subset of sessions with zero accuracy
            vetted_these_per_corr=[];
            ii_s=0;
            for ii_sessions=1:size(these_per_corr,1)
                if sum(these_per_corr(ii_sessions,:)==0)<200
                    ii_s=ii_s+1;
                    vetted_these_per_corr(ii_s,:)=these_per_corr(ii_sessions,:);
                end
            end

            these_per_corr=[];
            these_per_corr=vetted_these_per_corr;


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

            per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy=mean(these_per_corr,1)';
            per_mouse_acc.group(grNo).mouse(mouseNo).t_mean_accuracy=time_span';

            if show_per_mouse==1
                for ii_session=1:size(these_per_corr,1)
                    plot(time_span',smoothdata(these_per_corr(ii_session,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
                end

                plot(time_span',mean(these_per_corr,1)', 'k','LineWidth',1.5);

                ylim([0.3 1.2])
                this_ylim=ylim;
                plot([0 0],this_ylim,'-k')
            end
        end

        if no_pcorr==1
            xlim([-20 30])
        else
            xlim([-10 20])
        end
        xlabel('Time(sec)')



        if no_pcorr==1
            title(['Decoding accuracy for ' handles.group_names{grNo}])
        else
            this_grNo=floor((grNo-1)/4)+1;
            ii_pcorr=grNo-4*(this_grNo-1);
            title(['Decoding accuracy for ' handles.group_names{these_groups(this_grNo)} ' mouse ' handles_out2.mouse_names{mouseNo} ' ' per_names{ii_pcorr}])
        end
    end
end


%Plot the overall accuracy calculated from per mouse decoding
for grNo=1:no_pcorr*length(these_groups)

    ii_m_included=0;
    these_mean_accuracy=[];
    these_mean_accuracy_sh=[];
    for mouseNo=1:length(handles_out2.mouse_names)
        if ~isempty(per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy)
            ii_m_included=ii_m_included+1;
            these_mean_accuracy(ii_m_included,:)=per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy;
            these_mean_accuracy_sh(ii_m_included,:)=per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy_sh;
        end
    end

    if size(these_mean_accuracy,1)>2


        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        CIpv = bootci(1000, @mean, these_mean_accuracy_sh);
        meanpv=mean(these_mean_accuracy_sh,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;

        [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_accuracy_sh,1)', CIpv', 'r');

        CIpv = bootci(1000, @mean, these_mean_accuracy);
        meanpv=mean(these_mean_accuracy,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_accuracy,1)', CIpv', 'k');

        for ii_session=1:size(these_mean_accuracy,1)
            plot(time_span',smoothdata(these_mean_accuracy(ii_session,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
        end

        plot(time_span',mean(these_mean_accuracy,1)', 'k','LineWidth',1.5);

        ylim([0.3 1.2])
        this_ylim=ylim;
        


        %Odor on markers
        plot([0 0],this_ylim,'-k')
        odorhl=plot([0 delta_odor],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
        plot([delta_odor delta_odor],this_ylim,'-k')

        %Reinforcement markers
        plot([delta_odor_on_reinf_on delta_odor_on_reinf_on],this_ylim,'-r')
        reinfhl=plot([delta_odor_on_reinf_on delta_odor_on_reinf_on+delta_reinf],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
        plot([delta_odor_on_reinf_on+delta_reinf delta_odor_on_reinf_on+delta_reinf],this_ylim,'-r')


        xlim([-10 20])

        this_grNo=floor((grNo-1)/4)+1;
        ii_pcorr=grNo-4*(this_grNo-1);

        title(['Decoding accuracy calculated per mouse for ' handles.group_names{these_groups(this_grNo)} ' ' per_names{ii_pcorr}])

    end

end

%Plot a bar graph with separate forward and reversed and calculate statistics
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

id_acc_ii=0;
input_acc_data=[];

% for grNo=1:no_pcorr*length(these_groups)
bar_offset=0;

%Time windows
pre_t=[-1 0];
odor_t=[3.1 4.1];
reinf_t=[4.5 5.5];

edges=[0:0.05:1];
rand_offset=0.7;

switch is_Fabio
    case 0
        these_groups_out=[2 3 4 5 6 7 8];
    case 1
        these_groups_out=[grNo1 grNo2];
    case 2
        these_groups_out=[1];
end


for grNo=these_groups_out

    ii_m_included=0;
    these_mean_accuracy=[];
    these_mean_accuracy_sh=[];
    mouse_nos=[];
    for mouseNo=1:length(handles_out2.mouse_names)
        if ~isempty(per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy)
            ii_m_included=ii_m_included+1;
            these_mean_accuracy(ii_m_included,:)=per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy;
            these_mean_accuracy_sh(ii_m_included,:)=per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy_sh;
            mouse_nos(ii_m_included)=mouseNo;
        end
    end

    if size(these_mean_accuracy,1)>0

        %For sh we use the odor window
        these_accs=[];
        these_accs=mean(these_mean_accuracy_sh(:,(time_span>=odor_t(1))&(time_span<=odor_t(2))),2);
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_accs...
            ,edges,bar_offset,rand_offset,'k','k',3);

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        if grNo<5
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=0*ones(1,length(these_accs));
        else
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=1*ones(1,length(these_accs));
        end
        pcorr=rem(grNo-1,4)+1;
        glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=0*ones(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
        glm_acc_ii=glm_acc_ii+length(these_accs);

        id_acc_ii=id_acc_ii+1;
        input_acc_data(id_acc_ii).data=these_accs;
        input_acc_data(id_acc_ii).description=['Shuffled ' fr_per_names{grNo}];

        bar_offset=bar_offset+1;

        %Pre window
        these_accs=[];
        these_accs=mean(these_mean_accuracy(:,(time_span>=pre_t(1))&(time_span<=pre_t(2))),2);
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_accs...
            ,edges,bar_offset,rand_offset,'k','k',3);

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        if grNo<5
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=0*ones(1,length(these_accs));
        else
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=1*ones(1,length(these_accs));
        end
        pcorr=rem(grNo-1,4)+1;
        glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=1*ones(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
        glm_acc_ii=glm_acc_ii+length(these_accs);


        id_acc_ii=id_acc_ii+1;
        input_acc_data(id_acc_ii).data=these_accs;
        input_acc_data(id_acc_ii).description=['Pre ' fr_per_names{grNo}];


        bar_offset=bar_offset+1;

        %Odor window
        these_accs=[];
        these_accs=mean(these_mean_accuracy(:,(time_span>=odor_t(1))&(time_span<=odor_t(2))),2);
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_accs...
            ,edges,bar_offset,rand_offset,'k','k',3);

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        if grNo<5
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=0*ones(1,length(these_accs));
        else
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=1*ones(1,length(these_accs));
        end
        pcorr=rem(grNo-1,4)+1;
        glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=2*ones(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
        glm_acc_ii=glm_acc_ii+length(these_accs);


        id_acc_ii=id_acc_ii+1;
        input_acc_data(id_acc_ii).data=these_accs;
        input_acc_data(id_acc_ii).description=['Odor ' fr_per_names{grNo}];


        bar_offset=bar_offset+1;

        %Reinf window
        these_accs=[];
        these_accs=mean(these_mean_accuracy(:,(time_span>=reinf_t(1))&(time_span<=reinf_t(2))),2);
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_accs...
            ,edges,bar_offset,rand_offset,'k','k',3);

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        if grNo<5
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=0*ones(1,length(these_accs));
        else
            glm_acc.fwd_rev(glm_acc_ii+1:glm_acc_ii+length(these_accs))=1*ones(1,length(these_accs));
        end
        pcorr=rem(grNo-1,4)+1;
        glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
        glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=3*ones(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
        glm_acc_ii=glm_acc_ii+length(these_accs);


        id_acc_ii=id_acc_ii+1;
        input_acc_data(id_acc_ii).data=these_accs;
        input_acc_data(id_acc_ii).description=['Reinforcement ' fr_per_names{grNo}];

        bar_offset=bar_offset+2;

    end

end

xticks([0 1 2 3 5 6 7 8 10 11 12 13 20 21 22 23 25 26 27 28])
xticklabels({'Shuffled','Pre','Odor','Reinforcement','Shuffled','Pre','Odor','Reinforcement','Shuffled','Pre','Odor','Reinforcement'...
    ,'Shuffled','Pre','Odor','Reinforcement','Shuffled','Pre','Odor','Reinforcement','Shuffled','Pre','Odor','Reinforcement'})
xtickangle(45)
title(['Prediction accuracy for ' handles_out2.classifier_names{iiMLalgo}])
ylabel('Accuracy')
ylim([0.4 1.1])
xlim([-1 30])

text(1,1,'Fwd 40-65%')
text(6,1,'Fwd 65-80%')
text(11,1,'Fwd >80%')
text(16,1,'Rev 40-65%')
text(21,1,'Rev 65-80%')
text(26,1,'Rev >80%')

%Perform the glm
fprintf(1, ['\nglm for decoding accuracy\n'])
fprintf(fileID, ['\nglm for decoding accuracy\n']);

tbl = table(glm_acc.data',glm_acc.fwd_rev',glm_acc.pcorr',glm_acc.window',...
    'VariableNames',{'accuracy','forward_vs_reversed','percent_correct','window'});
mdl = fitglm(tbl,'accuracy~forward_vs_reversed+percent_correct+window+forward_vs_reversed*percent_correct*window'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);

switch is_Fabio
    
    case 0
        
        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for decoding accuracy\n'])
        fprintf(fileID, ['\n\nRanksum or t-test p values for decoding accuracy\n']);
        
        
        [output_data] = drgMutiRanksumorTtest(input_acc_data, fileID,0);
        
        %Nested ANOVAN
        %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
        nesting=[0 0 0 0; ... % This line indicates that group factor is not nested in any other factor.
            0 0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
            0 0 0 0; ... % This line indicates that event is not nested in any other factor.
            1 1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
        % (the 1 in position 1 on the line indicates nesting under the first factor).
        figureNo=figureNo+1;
        
        [p anovanTbl stats]=anovan(glm_acc.data,{glm_acc.fwd_rev glm_acc.pcorr glm_acc.window glm_acc.mouse_nos},...
            'model','interaction',...
            'nested',nesting,...
            'varnames',{'forward_vs_reversed', 'percent_correct','window','mouse_no'});
        
        fprintf(fileID, ['\n\nNested ANOVAN for decoding accuracy\n']);
        drgWriteANOVANtbl(anovanTbl,fileID);
        fprintf(fileID, '\n\n');
    otherwise
        
end

%Plot a bar graph with merged forward and reversed and calculate statistics
if is_Fabio==0
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
    
    
    id_acc_ii=0;
    input_acc_data=[];
    
    % for grNo=1:no_pcorr*length(these_groups)
    bar_offset=0;
    
    %Time windows
    pre_t=[-1 0];
    odor_t=[3.1 4.1];
    reinf_t=[4.5 5.5];
    
    edges=[0:0.05:1];
    rand_offset=0.7;
    
    group_sets(1,:)=[2 6];
    group_sets(2,:)=[3 7];
    group_sets(3,:)=[4 8];
    
    per_corr_set_label{1}='40-65%%';
    per_corr_set_label{2}='65-80%%';
    per_corr_set_label{3}='>80%%';
    
    for group_set_no=1:3
        
        ii_m_included=0;
        these_mean_accuracy=[];
        these_mean_accuracy_sh=[];
        mouse_nos=[];
        for grNo=group_sets(group_set_no,:)
            for mouseNo=1:length(handles_out2.mouse_names)
                if ~isempty(per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy)
                    ii_m_included=ii_m_included+1;
                    these_mean_accuracy(ii_m_included,:)=per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy;
                    these_mean_accuracy_sh(ii_m_included,:)=per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy_sh;
                    mouse_nos(ii_m_included)=mouseNo;
                end
            end
        end
        
        if size(these_mean_accuracy,1)>0
            
            %For sh we use the odor window
            these_accs=[];
            these_accs=mean(these_mean_accuracy_sh(:,(time_span>=odor_t(1))&(time_span<=odor_t(2))),2);
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_accs...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
            pcorr=rem(grNo-1,4)+1
            glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
            glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=0*ones(1,length(these_accs));
            glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
            glm_acc_ii=glm_acc_ii+length(these_accs);
            
            id_acc_ii=id_acc_ii+1;
            input_acc_data(id_acc_ii).data=these_accs;
            input_acc_data(id_acc_ii).description=['Shuffled ' per_corr_set_label{group_set_no}];
            
            bar_offset=bar_offset+1;
            
            %Pre window
            these_accs=[];
            these_accs=mean(these_mean_accuracy(:,(time_span>=pre_t(1))&(time_span<=pre_t(2))),2);
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_accs...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
            pcorr=rem(grNo-1,4)+1;
            glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
            glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=1*ones(1,length(these_accs));
            glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
            glm_acc_ii=glm_acc_ii+length(these_accs);
            
            id_acc_ii=id_acc_ii+1;
            input_acc_data(id_acc_ii).data=these_accs;
            input_acc_data(id_acc_ii).description=['Pre ' per_corr_set_label{group_set_no}];
            
            bar_offset=bar_offset+1;
            
            %Odor window
            these_accs=[];
            these_accs=mean(these_mean_accuracy(:,(time_span>=odor_t(1))&(time_span<=odor_t(2))),2);
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_accs...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
            pcorr=rem(grNo-1,4)+1
            glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
            glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=2*ones(1,length(these_accs));
            glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
            glm_acc_ii=glm_acc_ii+length(these_accs);
            
            id_acc_ii=id_acc_ii+1;
            input_acc_data(id_acc_ii).data=these_accs;
            input_acc_data(id_acc_ii).description=['Odor ' per_corr_set_label{group_set_no}];
            
            bar_offset=bar_offset+1;
            
            %Reinf window
            these_accs=[];
            these_accs=mean(these_mean_accuracy(:,(time_span>=reinf_t(1))&(time_span<=reinf_t(2))),2);
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[180/255 180/255 180/255])
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_accs...
                ,edges,bar_offset,rand_offset,'k','k',3);
            
            glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
            pcorr=rem(grNo-1,4)+1
            glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=pcorr*ones(1,length(these_accs));
            glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=3*ones(1,length(these_accs));
            glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouse_nos;
            glm_acc_ii=glm_acc_ii+length(these_accs);
            
            id_acc_ii=id_acc_ii+1;
            input_acc_data(id_acc_ii).data=these_accs;
            input_acc_data(id_acc_ii).description=['Reinforcement ' per_corr_set_label{group_set_no}];
            
            bar_offset=bar_offset+2;
            
        end
        
    end
    
    xticks([0 1 2 3 5 6 7 8 10 11 12 13 ])
    xticklabels({'Shuffled','Pre','Odor','Reinforcement','Shuffled','Pre','Odor','Reinforcement','Shuffled','Pre','Odor','Reinforcement'...
        })
    xtickangle(45)
    title(['Prediction accuracy for ' handles_out2.classifier_names{iiMLalgo}])
    ylabel('Accuracy')
    ylim([0.4 1.1])
    xlim([-1 15])
    
    text(1,1,'0-65%')
    text(6,1,'65-80%')
    text(11,1,'>80%')
    
    %Perform the glm
    fprintf(1, ['\nglm for decoding accuracy\n'])
    fprintf(fileID, ['\nglm for decoding accuracy\n']);
    
    tbl = table(glm_acc.data',glm_acc.pcorr',glm_acc.window',...
        'VariableNames',{'accuracy','percent_correct','window'});
    mdl = fitglm(tbl,'accuracy~percent_correct+window+percent_correct*window'...
        ,'CategoricalVars',[2,3])
    
    
    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);
    
    fprintf(fileID,'%s\n', txt);
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for decoding accuracy\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for decoding accuracy\n']);
    
    
    [output_data] = drgMutiRanksumorTtest(input_acc_data, fileID,0);
    
    %Nested ANOVAN
    %https://www.mathworks.com/matlabcentral/answers/491365-within-between-subjects-in-anovan
    nesting=[0 0 0; ... % This line indicates that group factor is not nested in any other factor.
        0 0 0; ... % This line indicates that perCorr is not nested in any other factor.
        1 1 0];    % This line indicates that mouse_no (the last factor) is nested under group, perCorr and event
    % (the 1 in position 1 on the line indicates nesting under the first factor).
    figureNo=figureNo+1;
    
    [p anovanTbl stats]=anovan(glm_acc.data,{glm_acc.pcorr glm_acc.window glm_acc.mouse_nos},...
        'model','interaction',...
        'nested',nesting,...
        'varnames',{'percent_correct','window','mouse_no'});
    
    fprintf(fileID, ['\n\nNested ANOVAN for decoding accuracy\n']);
    drgWriteANOVANtbl(anovanTbl,fileID);
    fprintf(fileID, '\n\n');
    
end

pfft=1;
% %Plot decoding accuracy for forward proficient
% figureNo = figureNo + 1;
% try
%     close(figureNo)
% catch
% end
% hFig=figure(figureNo);
% hold on
% 
% ax=gca;ax.LineWidth=3;
% set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
% 
% 
% %Do forward
% grNo=grNo1;
% 
% [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
% 
% if ~isempty(ii_tspan)
%     time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
% 
% 
%     %Extrapolate all points onto the longest ii_tspan
%     time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
% 
%     for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%         if ii_f==ii_file
%             these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_file,1:ii_tspan);
%         else
%             this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%             this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
% 
%             for ii_tsp=1:ii_tspan
%                 if time_span(ii_tsp)<this_time_span(1)
%                     these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,1);
% 
%                 else
%                     if time_span(ii_tsp)>this_time_span(end)
%                         these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,end);
%                     else
%                         ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                         ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                         these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_0)+...
%                             (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict_sh(ii_f,ii_0))...
%                             *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                     end
%                 end
%             end
% 
%         end
% 
%     end
% 
% 
%     %For some reason whe have a subset of sessions with zero accuracy
%     vetted_these_per_corr=[];
%     ii_s=0;
%     for ii_sessions=1:size(these_per_corr,1)
%         if sum(these_per_corr(ii_sessions,:)==0)<200
%             ii_s=ii_s+1;
%             vetted_these_per_corr(ii_s,:)=these_per_corr(ii_sessions,:);
%         end
%     end
% 
%     these_per_corr=[];
%     these_per_corr=vetted_these_per_corr;
% 
%     if size(these_per_corr,1)>2
% 
%         CIpv = bootci(1000, @mean, these_per_corr);
%         meanpv=mean(these_per_corr,1);
%         CIpv(1,:)=meanpv-CIpv(1,:);
%         CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%         [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'k');
%     else
%         if size(these_per_corr,1)>1
%             plot(time_span',mean(these_per_corr,1)', 'k');
%         else
%             plot(time_span',these_per_corr', 'k');
%         end
% 
%     end
% 
% 
%     for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%         if ii_f==ii_file
%             these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_file,1:ii_tspan);
%         else
%             this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%             this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
% 
%             for ii_tsp=1:ii_tspan
%                 if time_span(ii_tsp)<this_time_span(1)
%                     these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,1);
% 
%                 else
%                     if time_span(ii_tsp)>this_time_span(end)
%                         these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,end);
%                     else
%                         ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                         ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                         these_per_corr(ii_f,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_0)+...
%                             (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_correct_predict(ii_f,ii_0))...
%                             *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                     end
%                 end
%             end
% 
%         end
% 
%     end
% 
%     %For some reason whe have a subset of sessions with zero accuracy
%     vetted_these_per_corr=[];
%     ii_s=0;
%     for ii_sessions=1:size(these_per_corr,1)
%         if sum(these_per_corr(ii_sessions,:)==0)<200
%             ii_s=ii_s+1;
%             vetted_these_per_corr(ii_s,:)=these_per_corr(ii_sessions,:);
%         end
%     end
% 
%     these_per_corr=[];
%     these_per_corr=vetted_these_per_corr;
% 
%     if size(these_per_corr,1)>2
% 
%         CIpv = bootci(1000, @mean, these_per_corr);
%         meanpv=mean(these_per_corr,1);
%         CIpv(1,:)=meanpv-CIpv(1,:);
%         CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%         [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'b');
%     else
%         if size(these_per_corr,1)>1
%             plot(time_span',mean(these_per_corr,1)', 'b');
%         else
%             plot(time_span',these_per_corr', 'b');
%         end
% 
%     end
% 
%     for ii_session=1:size(these_per_corr,1)
%         plot(time_span',smoothdata(these_per_corr(ii_session,:)','gaussian',100),'Color',[120/255 140/255 255/255],'LineWidth',1)
%     end
% 
%     plot(time_span',mean(these_per_corr,1)', 'b','LineWidth',1.5);
% 
%     ylim([0.3 1.2])
%     this_ylim=ylim;
%         plot([0 0],this_ylim,'-k')
%     if no_pcorr==1
%         xlim([-20 30])
%     else
%         xlim([-10 20])
%     end
%     xlabel('Time(sec)')
% 
% 
%     title(['Decoding accuracy for ' grNo1_label])
% end

%Label prediction for S+ and S-
per_mouse_pred=[];

for mouseNo=1:length(handles_out2.mouse_names)
    for grNo=these_groups_out

 
        per_mouse_pred.group(grNo).mouse(mouseNo).mean_accuracy=[];
        per_mouse_pred.group(grNo).mouse(mouseNo).t_mean_accuracy=[];
        per_mouse_pred.group(grNo).mouse(mouseNo).mean_accuracy_sh=[];
        per_mouse_pred.group(grNo).mouse(mouseNo).t_mean_accuracy_sh=[];
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

        [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);

        if ~isempty(ii_tspan)


            %Extrapolate all points onto the longest ii_tspan
            time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);

            these_moving_mean_sm=[];
            ii_included=0;
            time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
            these_mice=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mouseNo;
            for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
                if these_mice(ii_f)==mouseNo
                    ii_included=ii_included+1;
                    if ii_f==ii_file
                        these_moving_mean_sm(ii_included,:)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii_file,1:ii_tspan);
                    else
                        this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                        this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);

                        for ii_tsp=1:ii_tspan
                            if time_span(ii_tsp)<this_time_span(1)
                                these_moving_mean_sm(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii_f,1);

                            else
                                if time_span(ii_tsp)>this_time_span(end)
                                    these_moving_mean_sm(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii_f,end);
                                else
                                    ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                                    ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                                    these_moving_mean_sm(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii_f,ii_0)+...
                                        (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(ii_f,ii_0))...
                                        *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                end
                            end
                        end

                    end

                end

            end

            %             these_moving_mean_sm=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(:,1:ii_tspan);
            if show_per_mouse==1
                if ~isempty(these_moving_mean_sm)

                    if size(these_moving_mean_sm,1)>2

                        CIpv = bootci(1000, @mean, these_moving_mean_sm);
                        meanpv=mean(these_moving_mean_sm,1);
                        CIpv(1,:)=meanpv-CIpv(1,:);
                        CIpv(2,:)=CIpv(2,:)-meanpv;


                        [hlpvl, hppvl] = boundedline(time_span',mean(these_moving_mean_sm,1)', CIpv', 'cmap',[158/255 31/255 99/255]);
                    else
                        if size(these_moving_mean_sm,1)>1
                            plot(time_span',mean(these_moving_mean_sm,1)', '-', 'Color', [158/255 31/255 99/255]);
                        else
                            plot(time_span',these_moving_mean_sm', '-', 'Color', [158/255 31/255 99/255]);
                        end

                    end
                end
            end

            per_mouse_pred.group(grNo).mouse(mouseNo).mean_moving_mean_sm=mean(these_moving_mean_sm,1)';
            per_mouse_pred.group(grNo).mouse(mouseNo).t_mean__moving_mean_sm=time_span';

            %Extrapolate all points onto the longest ii_tspan
            time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);

            these_moving_mean_sp=[];
            ii_included=0;
            for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
                if these_mice(ii_f)==mouseNo
                    ii_included=ii_included+1;
                    if ii_f==ii_file
                        these_moving_mean_sp(ii_included,:)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii_file,1:ii_tspan);
                    else
                        this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                        this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);

                        for ii_tsp=1:ii_tspan
                            if time_span(ii_tsp)<this_time_span(1)
                                these_moving_mean_sp(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii_f,1);

                            else
                                if time_span(ii_tsp)>this_time_span(end)
                                    these_moving_mean_sp(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii_f,end);
                                else
                                    ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                                    ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                                    these_moving_mean_sp(ii_included,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii_f,ii_0)+...
                                        (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii_f,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(ii_f,ii_0))...
                                        *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                end
                            end
                        end

                    end
                end
            end

            %             these_moving_mean_sp=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(:,1:ii_tspan);

            if show_per_mouse==1
                if ~isempty(these_moving_mean_sp)
                    if size(these_moving_mean_sp,1)>2

                        CIpv = bootci(1000, @mean, these_moving_mean_sp);
                        meanpv=mean(these_moving_mean_sp,1);
                        CIpv(1,:)=meanpv-CIpv(1,:);
                        CIpv(2,:)=CIpv(2,:)-meanpv;


                        [hlpvl, hppvl] = boundedline(time_span',mean(these_moving_mean_sp,1)', CIpv', 'cmap',[0 114/255 178/255]);
                    else
                        if size(these_moving_mean_sp,1)>1
                            plot(time_span',mean(these_moving_mean_sp,1)', '-', 'Color', [0 114/255 178/255]);
                        else
                            plot(time_span',these_moving_mean_sp', '-', 'Color', [0 114/255 178/255]);
                        end

                    end
                end
            end

            per_mouse_pred.group(grNo).mouse(mouseNo).mean_moving_mean_sp=mean(these_moving_mean_sp,1)';
            per_mouse_pred.group(grNo).mouse(mouseNo).t_mean__moving_mean_sp=time_span';

            if show_per_mouse==1
                text(30,0.75,'S-','Color',[158/255 31/255 99/255])
                text(30,0.85,'S+','Color',[0 114/255 178/255])

                ylim([0 1])
                this_ylim=ylim;
                plot([0 0],this_ylim,'-k')
                if no_pcorr==1
                    xlim([-20 30])
                else
                    xlim([-10 20])
                end
                xlabel('Time(sec)')



                if no_pcorr==1
                    title(['Label prediction for ' handles.group_names{grNo}])
                else
                    this_grNo=floor((grNo-1)/4)+1;
                    ii_pcorr=grNo-4*(this_grNo-1);
                    title(['Label prediction for ' handles.group_names{these_groups_out(this_grNo)} ' mouse ' handles_out2.mouse_names{mouseNo} ' ' per_names{ii_pcorr}])
                end
            end
        end
    end


end

 
%Plot the overall predictions calculated from per mouse decoding
for grNo=these_groups_out

    ii_m_included=0;
    these_mean_moving_mean_sp=[];
    these_mean_moving_mean_sm=[];
    for mouseNo=1:length(handles_out2.mouse_names)
        if ~isempty(per_mouse_pred.group(grNo).mouse(mouseNo).mean_moving_mean_sp)
            ii_m_included=ii_m_included+1;
            these_mean_moving_mean_sp(ii_m_included,:)=per_mouse_pred.group(grNo).mouse(mouseNo).mean_moving_mean_sp;
            these_mean_moving_mean_sm(ii_m_included,:)=per_mouse_pred.group(grNo).mouse(mouseNo).mean_moving_mean_sm;
        end
    end

    if size(these_mean_moving_mean_sp,1)>2


        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        CIpv = bootci(1000, @mean, these_mean_moving_mean_sp);
        meanpv=mean(these_mean_moving_mean_sp,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;

        [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_moving_mean_sp,1)', CIpv', 'cmap',[0 114/255 178/255]);

        CIpv = bootci(1000, @mean, these_mean_moving_mean_sm);
        meanpv=mean(these_mean_moving_mean_sm,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(time_span',mean(these_mean_moving_mean_sm,1)', CIpv', 'cmap',[158/255 31/255 99/255]);

%         for ii_session=1:size(these_mean_accuracy,1)
%             plot(time_span',smoothdata(these_mean_accuracy(ii_session,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
%         end
% 
%         plot(time_span',mean(these_mean_accuracy,1)', 'k','LineWidth',1.5);

        ylim([-0.1 1.1])
        this_ylim=ylim;
       

        xlim([-10 20])

        this_grNo=floor((grNo-1)/4)+1;
        ii_pcorr=grNo-4*(this_grNo-1);

             %Odor on markers
        plot([0 0],this_ylim,'-k')
        odorhl=plot([0 delta_odor],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
        plot([delta_odor delta_odor],this_ylim,'-k')

        %Reinforcement markers
        plot([delta_odor_on_reinf_on delta_odor_on_reinf_on],this_ylim,'-r')
        reinfhl=plot([delta_odor_on_reinf_on delta_odor_on_reinf_on+delta_reinf],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
        plot([delta_odor_on_reinf_on+delta_reinf delta_odor_on_reinf_on+delta_reinf],this_ylim,'-r')



        title(['Prediction for ' handles.group_names{this_grNo} ' mouse ' per_names{ii_pcorr}])

    end

end
% 
% %Plot the posterior probabilities for Sp (scores(:,2) and Sm (scores(:,1))
% for iiMLalgo=handles.MLalgo_to_use
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .3 .6])
% 
%     ii_plot=0;
%     for grNo=1:length(these_groups)
% 
%         [ii_tspan,ii_file]=max(handles_out2.group_no(grNo).ii_time_span);
% 
%         if ~isempty(ii_tspan)
%             try
%                 %                 time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
% 
%                 ii_plot=ii_plot+1;
%                 subplot(6,2,ii_plot)
%                 hold on
% 
%                 %Extrapolate all points onto the longest ii_tspan
%                 time_span=handles_out2.group_no(grNo).time_span_euclid(ii_file,1:ii_tspan);
% 
%                 this_mean_per_trial_scores=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp;
%                 these_mean_per_trial_scores=zeros(size(this_mean_per_trial_scores,1),size(this_mean_per_trial_scores,2),ii_tspan);
%                 these_scores_sm=zeros(size(this_mean_per_trial_scores,1),ii_tspan);
%                 these_scores_sp=zeros(size(this_mean_per_trial_scores,1),ii_tspan);
% 
%                 for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%                     if ii_f==ii_file
%                         these_mean_per_trial_scores=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii_file,:,1:ii_tspan);
%                     else
%                         this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%                         this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
% 
%                         for ii_tsp=1:ii_tspan
%                             if time_span(ii_tsp)<this_time_span(1)
%                                 these_mean_per_trial_scores(ii_f,:,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii_f,:,1);
% 
%                             else
%                                 if time_span(ii_tsp)>this_time_span(end)
%                                     these_mean_per_trial_scores(ii_f,:,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii_f,:,end);
%                                 else
%                                     ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                                     ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                                     these_mean_per_trial_scores(ii_f,:,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii_f,:,ii_0)+...
%                                         (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii_f,:,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp(ii_f,:,ii_0))...
%                                         *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                                 end
%                             end
%                         end
% 
%                     end
% 
%                 end
% 
%                 %                 this_mean_per_trial_scores_sp=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sp;
%                 %                 this_mean_per_trial_scores_sp=this_mean_per_trial_scores_sp(:,:,1:length(time_span));
% 
% 
%                 these_scores_sm(:,:)=these_mean_per_trial_scores(:,1,:);
%                 these_scores_sp(:,:)=these_mean_per_trial_scores(:,2,:);
% 
%                 if size(these_scores_sm,1)>2
% 
%                     CIsm = bootci(1000, @mean, these_scores_sm);
% 
%                     meansm=mean(these_scores_sm,1);
%                     CIsm(1,:)=meansm-CIsm(1,:);
%                     CIsm(2,:)=CIsm(2,:)-meansm;
% 
% 
%                     [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
% 
%                 else
%                     plot(time_span',these_scores_sm', '-','Color',[158/255 31/255 99/255]);
%                 end
% 
%                 %                 these_scores_sp=zeros(size(this_mean_per_trial_scores_sp,1),size(this_mean_per_trial_scores_sp,3));
%                 %                 these_scores_sp(:,:)=this_mean_per_trial_scores_sp(:,2,:);
% 
%                 if size(these_scores_sp,1)>2
%                     CIsp = bootci(1000, @mean, these_scores_sp);
%                     meansp=mean(these_scores_sp,1);
%                     CIsp(1,:)=meansp-CIsp(1,:);
%                     CIsp(2,:)=CIsp(2,:)-meansp;
% 
%                     [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);
%                 else
%                     plot(time_span',these_scores_sp', '-','Color',[0 114/255 178/255]);
%                 end
% 
% 
%                 if (size(these_scores_sp,1)>2)&(size(these_scores_sm,1)>2)
%                     plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
%                     plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
%                 end
% 
%                 text(50,0.75,'S-','Color',[158/255 31/255 99/255])
%                 text(50,0.65,'S+','Color',[0 114/255 178/255])
% 
%                 ylim([0 1])
%                 title(['S+ ' handles.group_names{these_groups(grNo)}])
%                 xlabel('Time(sec)')
%                 ylabel('Posterior p')
% 
%                 ii_plot=ii_plot+1;
%                 subplot(6,2,ii_plot)
%                 hold on
% 
%                 %                 this_mean_per_trial_scores_sm=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm;
%                 %                 this_mean_per_trial_scores_sm=this_mean_per_trial_scores_sm(:,:,1:length(time_span));
%                 %
%                 %                 these_scores_sm=zeros(size(this_mean_per_trial_scores_sm,1),size(this_mean_per_trial_scores_sm,3));
%                 %                 these_scores_sm(:,:)=this_mean_per_trial_scores_sm(:,1,:);
% 
%                 this_mean_per_trial_scores=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm;
%                 these_mean_per_trial_scores=zeros(size(this_mean_per_trial_scores,1),size(this_mean_per_trial_scores,2),ii_tspan);
%                 these_scores_sm=zeros(size(this_mean_per_trial_scores,1),ii_tspan);
%                 these_scores_sp=zeros(size(this_mean_per_trial_scores,1),ii_tspan);
% 
%                 for ii_f=1:size(handles_out2.group_no(grNo).ii_time_span,1)
%                     if ii_f==ii_file
%                         these_mean_per_trial_scores=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii_file,:,1:ii_tspan);
%                     else
%                         this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
%                         this_time_span=handles_out2.group_no(grNo).time_span_euclid(ii_f,1:this_ii_tspan);
% 
%                         for ii_tsp=1:ii_tspan
%                             if time_span(ii_tsp)<this_time_span(1)
%                                 these_mean_per_trial_scores(ii_f,:,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii_f,:,1);
% 
%                             else
%                                 if time_span(ii_tsp)>this_time_span(end)
%                                     these_mean_per_trial_scores(ii_f,:,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii_f,:,end);
%                                 else
%                                     ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
%                                     ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
%                                     these_mean_per_trial_scores(ii_f,:,ii_tsp)=handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii_f,:,ii_0)+...
%                                         (handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii_f,:,ii_1)-handles_out2.group_no(grNo).ii_thr(ii_thr).MLalgo(iiMLalgo).mean_per_trial_scores_sm(ii_f,:,ii_0))...
%                                         *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
%                                 end
%                             end
%                         end
% 
%                     end
% 
%                 end
% 
%                 these_scores_sm(:,:)=these_mean_per_trial_scores(:,1,:);
%                 these_scores_sp(:,:)=these_mean_per_trial_scores(:,2,:);
% 
%                 if size(these_scores_sm,1)>2
%                     CIsm = bootci(1000, @mean, these_scores_sm);
%                     meansm=mean(these_scores_sm,1);
%                     CIsm(1,:)=meansm-CIsm(1,:);
%                     CIsm(2,:)=CIsm(2,:)-meansm;
% 
%                     [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
%                 else
%                     plot(time_span',these_scores_sm', '-','Color',[158/255 31/255 99/255]);
%                 end
% 
% %                 these_scores_sp=zeros(size(this_mean_per_trial_scores_sm,1),size(this_mean_per_trial_scores_sm,3));
% %                 these_scores_sp(:,:)=this_mean_per_trial_scores_sm(:,2,:);
% 
%                 if size(these_scores_sp,1)>2
%                     CIsp = bootci(1000, @mean, these_scores_sp);
%                     meansp=mean(these_scores_sp,1);
%                     CIsp(1,:)=meansp-CIsp(1,:);
%                     CIsp(2,:)=CIsp(2,:)-meansp;
% 
%                     [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);
%                 else
%                     plot(time_span',these_scores_sp', '-','Color',[0 114/255 178/255]);
%                 end
% 
% 
% 
%                 if (size(these_scores_sp,1)>2)&(size(these_scores_sm,1)>2)
%                     plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
%                     plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
%                 end
% 
%                 text(50,0.75,'S-','Color',[158/255 31/255 99/255])
%                 text(50,0.65,'S+','Color',[0 114/255 178/255])
% 
%                 ylim([0 1])
%                 title(['S- ' handles.group_names{these_groups(grNo)}])
%                 xlabel('Time(sec)')
%                 ylabel('Posterior p')
%             catch
%             end
%         end
%     end
%     sgtitle(['Posterior probability for ' handles_out2.classifier_names{iiMLalgo} ])
% end
% 

% 
% %Average decoding accuracy first and last session
% try
%     for iiMLalgo=handles.MLalgo_to_use
%         figureNo = figureNo + 1;
%         try
%             close(figureNo)
%         catch
%         end
%         hFig=figure(figureNo);
%         hold on
% 
%         ax=gca;ax.LineWidth=3;
%         set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%         length_corr=[];
%         for mouseNo=1:size(first_last,1)
%             length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).first_mean_correct_predict);
%         end
% 
%         first_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
%         for mouseNo=1:size(first_last,1)
%             first_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).first_mean_correct_predict;
%         end
% 
%         these_per_corr=first_mean_correct_predict;
%         this_time_spanf=zeros(1,size(these_per_corr,2));
%         this_time_spanf(1)=time_span(1);
%         dt=time_span(2)-time_span(1);
%         for ii=2:length(this_time_spanf)
%             this_time_spanf(ii)=this_time_spanf(ii-1)+dt;
%         end
% 
% 
%         if size(these_per_corr,1)>2
% 
%             CIpv = bootci(1000, @mean, these_per_corr);
%             meanpv=mean(these_per_corr,1);
%             CIpv(1,:)=meanpv-CIpv(1,:);
%             CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%             [hlpvl, hppvl] = boundedline(this_time_spanf',mean(these_per_corr,1)', CIpv', 'm');
%         else
%             if size(these_per_corr,1)>1
%                 plot(this_time_spanf',mean(these_per_corr,1)', 'm');
%             else
%                 plot(this_time_spanf',these_per_corr', 'm');
%             end
% 
%         end
% 
% 
% 
% 
%         length_corr=[];
%         for mouseNo=1:size(first_last,1)
%             length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).last_mean_correct_predict);
%         end
% 
%         last_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
%         for mouseNo=1:size(first_last,1)
%             last_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).last_mean_correct_predict;
%         end
% 
%         these_per_corr=last_mean_correct_predict;
%         this_time_span=zeros(1,size(these_per_corr,2));
%         this_time_span(1)=time_span(1);
%         for ii=2:length(this_time_span)
%             this_time_span(ii)=this_time_span(ii-1)+dt;
%         end
% 
%         if size(these_per_corr,1)>2
% 
%             CIpv = bootci(1000, @mean, these_per_corr);
%             meanpv=mean(these_per_corr,1);
%             CIpv(1,:)=meanpv-CIpv(1,:);
%             CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%             [hlpvl, hppvl] = boundedline(this_time_span',mean(these_per_corr,1)', CIpv', 'b');
%         else
%             if size(these_per_corr,1)>1
%                 plot(this_time_span',mean(these_per_corr,1)', 'b');
%             else
%                 plot(this_time_span',these_per_corr', 'b');
%             end
% 
%         end
% 
%         text(20,0.8,'First','Color','m')
%         text(20,0.9,'Last','Color','b')
%         plot(this_time_spanf',mean(first_mean_correct_predict,1)', 'm');
%         ylabel('Accuracy')
%         xlabel('Time(sec)')
%         title(['Accuracy ' handles_out2.classifier_names{iiMLalgo}])
%     end
% catch
% end
% 
% try
%     %Decoding accuracy first and last session per mouse
%     for iiMLalgo=handles.MLalgo_to_use
%         figureNo = figureNo + 1;
%         try
%             close(figureNo)
%         catch
%         end
%         hFig=figure(figureNo);
%         hold on
% 
%         ax=gca;ax.LineWidth=3;
%         set(hFig, 'units','normalized','position',[.2 .2 .3 .7])
% 
% 
%         length_corr=[];
%         for mouseNo=1:size(first_last,1)
%             length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).first_mean_correct_predict);
%         end
% 
%         first_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
%         for mouseNo=1:size(first_last,1)
%             first_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).first_mean_correct_predict;
%         end
% 
%         these_per_corr=first_mean_correct_predict;
%         this_time_spanf=zeros(1,size(these_per_corr,2));
%         this_time_spanf(1)=time_span(1);
%         dt=time_span(2)-time_span(1);
%         for ii=2:length(this_time_spanf)
%             this_time_spanf(ii)=this_time_spanf(ii-1)+dt;
%         end
% 
% 
%         length_corr=[];
%         for mouseNo=1:size(first_last,1)
%             length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).last_mean_correct_predict);
%         end
% 
%         last_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
%         for mouseNo=1:size(first_last,1)
%             last_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_thr).MLalgo(iiMLalgo).last_mean_correct_predict;
%         end
% 
%         these_per_corr=last_mean_correct_predict;
%         this_time_span=zeros(1,size(these_per_corr,2));
%         this_time_span(1)=time_span(1);
%         for ii=2:length(this_time_span)
%             this_time_span(ii)=this_time_span(ii-1)+dt;
%         end
% 
%         %Now plot them
%         for mouseNo=1:size(first_last,1)
%             subplot(size(first_last,1),1,mouseNo)
%             hold on
%             plot(this_time_spanf',first_mean_correct_predict(mouseNo,:),'-m')
%             plot(this_time_span',last_mean_correct_predict(mouseNo,:),'-b')
% 
% 
%             if mouseNo==1
%                 text(20,0.8,'First','Color','m')
%                 text(20,0.9,'Last','Color','b')
%             end
% 
%             ylabel('Accuracy')
%             xlabel('Time(sec)')
%             title(['Accuracy ' handles_out2.mouse_names{mouseNo}])
%         end
% 
%         sgtitle(['Accuracy ' handles_out2.classifier_names{iiMLalgo}])
%     end
% catch
% end
 
out_file=[choiceBatchPathName choiceFileName];
out_file=[out_file(1:end-2) '.mat'];
save(out_file,'handles_out2','handles','-v7.3')
fclose(fileID);

pfft=1;