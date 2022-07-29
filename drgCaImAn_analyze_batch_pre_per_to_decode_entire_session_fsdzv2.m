%drgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdzv2
close all
clear all



[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdzv2 run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

no_files=handles.no_files;
moving_mean_n=10;


%Now plot for each algorithm the prediction accuracy
handles_out2=[];

handles_out2.classifier_names{1}='Linear Discriminant';
handles_out2.classifier_names{2}='Support Vector Machine';
handles_out2.classifier_names{3}='Naive Bayes Classifier';
handles_out2.classifier_names{4}='Neural Network';
handles_out2.classifier_names{5}='Decision tree';
handles_out2.classifier_names{6}='Binomial glm';



%First and last sessions per mouse
first_last=[20 2;
    5 15;
    25 1;
    13 11;
    26 12];

handles_out2.first_last=first_last;

%Mouse names
handles_out2.mouse_names{1}='GRIN1';
handles_out2.mouse_names{2}='GRIN3';
handles_out2.mouse_names{3}='GRIN4';
handles_out2.mouse_names{4}='GRIN6';
handles_out2.mouse_names{5}='GRIN7';




if exist([choiceBatchPathName choiceFileName(1:end-2) '.mat'])==0
    %Process each file separately
    for grNo=1:max(handles.group)
        handles_out2.group_no(grNo).ii_euclid=0;
        for iiMLalgo=handles.MLalgo_to_use
            for ii_out=1:length(handles.p_threshold)
                handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=0;
            end
        end
    end

    for fileNo=1:no_files
        tic
        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};
        grNo=handles.group(fileNo);

        load([handles.PathName_out pre_per_FileName(1:end-4) handles.suffix_out])

        %Is this the first session?
        if sum(first_last(:,1)==fileNo)==1
            mouseNo=find(first_last(:,1)==fileNo);
            for iiMLalgo=handles.MLalgo_to_use
                for ii_out=1:length(handles_out.ii_out)
                    this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
                    handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict(1,1:length(this_mean_correct_predict))=this_mean_correct_predict;

                    this_mean_correct_predict_sh=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh);
                    handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict_sh(1,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;
                end
            end
        end

        %Is this the last session?
        if sum(first_last(:,2)==fileNo)
              mouseNo=find(first_last(:,2)==fileNo);
            for iiMLalgo=handles.MLalgo_to_use
                for ii_out=1:length(handles_out.ii_out)
                    this_mean_correct_predict=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict);
                    handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict(1,1:length(this_mean_correct_predict))=this_mean_correct_predict;

                    this_mean_correct_predict_sh=mean(handles_out.ii_out(ii_out).handles_out.MLalgo(iiMLalgo).this_correct_predict_sh);
                    handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict_sh(1,1:length(this_mean_correct_predict_sh))=this_mean_correct_predict_sh;
                end
            end
        end

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
                if (fileNo==3)&(ii_out==1)&(iiMLalgo==4)
                    pffft=1;
                end
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

            end
        end
        fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
    end
else
    load([choiceBatchPathName choiceFileName(1:end-2) '.mat'])
end

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

    for grNo=1:max(handles.group)

        for ii_thr=1:length(handles.p_threshold)
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

            end
        end
        bar_offset=bar_offset+1;
    end
    xticks([4.5:9:49.5])
    labels='xticklabels({';
    for ii_label=1:length(handles.group_names)
        labels=[labels '''' handles.group_names{ii_label} ''', '];
    end
    labels=[labels(1:end-2) '})'];
    eval(labels)
    title(['Predcition accuracy for ' handles_out2.classifier_names{iiMLalgo}])
    ylabel('Accuracy')
    ylim([0 1])
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

for grNo=1:max(handles.group)
    
    ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
    if ~isempty(ii_tspan)
        time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
        
        subplot(6,1,grNo)
        hold on
         
        these_meandFFsp=handles_out2.group_no(grNo).meandFFsp(1,1:ii_tspan);
        these_meandFFsm=handles_out2.group_no(grNo).meandFFsm(1,1:ii_tspan);
        
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
        xlim([-40 40])
        xlabel('Time(sec)')
        title(handles.group_names{grNo})
    end
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

for grNo=1:max(handles.group)
    
    ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
    if ~isempty(ii_tspan)
        time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
        
        subplot(6,1,grNo)
        hold on
        
        these_dist_euclid=handles_out2.group_no(grNo).dist_euclid(:,1:ii_tspan);
        
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
        xlim([-40 40])
        xlabel('Time(sec)')
        title(handles.group_names{grNo})
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

for grNo=1:max(handles.group)
    
    ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
    if ~isempty(ii_tspan)
        time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
        
        subplot(6,1,grNo)
        hold on
        
        these_KLdivergence=handles_out2.group_no(grNo).KLdivergence(:,1:ii_tspan);
        
        
        if size(these_KLdivergence,1)>2
            
            CIpv = bootci(1000, @mean, these_KLdivergence);
            meanpv=mean(these_KLdivergence,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;
            
            
            [hlpvl, hppvl] = boundedline(time_span',mean(these_KLdivergence,1)', CIpv', 'm');
        else
            if size(these_KLdivergence,1)>0
                plot(time_span',mean(these_KLdivergence,1)', 'm');
            end
            
        end
        xlim([-40 40])
        xlabel('Time(sec)')
        title(handles.group_names{grNo})
    end
end
sgtitle('KL divergence')

%Choose threshold
ii_out=1;

%Decoding accuracy
for iiMLalgo=handles.MLalgo_to_use
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .6])
    
    for grNo=1:max(handles.group)
        ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
        if ~isempty(ii_tspan)
            time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
            subplot(6,1,grNo)
            hold on
            
            these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh(:,1:ii_tspan);
            
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
            
            these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(:,1:ii_tspan);
      
            if size(these_per_corr,1)>2
                
                CIpv = bootci(1000, @mean, these_per_corr);
                meanpv=mean(these_per_corr,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;
                
                
                [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'm');
            else
                if size(these_per_corr,1)>1
                    plot(time_span',mean(these_per_corr,1)', 'm');
                else
                    plot(time_span',these_per_corr', 'm');
                end
                
            end
            xlim([-40 40])
            xlabel('Time(sec)')
            title(handles.group_names{grNo})
        end
    end
    sgtitle(['Percent correct ' handles_out2.classifier_names{iiMLalgo}])
end

%Decoding accuracy
for iiMLalgo=handles.MLalgo_to_use
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .6])
    
    for grNo=1:max(handles.group)
        ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
        if ~isempty(ii_tspan)
            time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
            subplot(6,1,grNo)
            hold on
            
            these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict_sh(:,1:ii_tspan);
            
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
            
            these_per_corr=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_correct_predict(:,1:ii_tspan);
      
            if size(these_per_corr,1)>2
                
                CIpv = bootci(1000, @mean, these_per_corr);
                meanpv=mean(these_per_corr,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;
                
                
                [hlpvl, hppvl] = boundedline(time_span',mean(these_per_corr,1)', CIpv', 'm');
            else
                if size(these_per_corr,1)>1
                    plot(time_span',mean(these_per_corr,1)', 'm');
                else
                    plot(time_span',these_per_corr', 'm');
                end
                
            end
            xlim([-40 40])
            xlabel('Time(sec)')
            title(handles.group_names{grNo})
        end
    end
    sgtitle(['Percent correct ' handles_out2.classifier_names{iiMLalgo}])
end

%Label prediction for S+ and S-
for iiMLalgo=handles.MLalgo_to_use
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .6])
    
    for grNo=1:max(handles.group)
        ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
        if ~isempty(ii_tspan)
            time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);
            subplot(6,1,grNo)
            hold on
            
            these_moving_mean_sm=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sm_timecourse(:,1:ii_tspan);
          
            
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
            
            these_moving_mean_sp=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_moving_mean_per_trial_sp_timecourse(:,1:ii_tspan);
          
            
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
            
            text(30,0.75,'S-','Color',[158/255 31/255 99/255])
            text(30,0.85,'S+','Color',[0 114/255 178/255])
            
            ylim([0 1])
            xlabel('Time(sec)')
            ylabel('Label prediction, S+=1, S-=0')
            
            title(handles.group_names{grNo})
        end
    end
    sgtitle(['Label prediction per trial for ' handles_out2.classifier_names{iiMLalgo} ])
    
end
 

%Plot the posterior probabilities for Sp (scores,:,2) and Sm (scores(:,1))
for iiMLalgo=handles.MLalgo_to_use
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end
    
    hFig = figure(figureNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .6])
    
    ii_plot=0;
    for grNo=1:max(handles.group)
        
        ii_tspan=min(handles_out2.group_no(grNo).ii_time_span);
        if ~isempty(ii_tspan)
            try
                time_span=handles_out2.group_no(grNo).time_span_euclid(1,1:ii_tspan);

                ii_plot=ii_plot+1;
                subplot(6,2,ii_plot)
                hold on

                this_mean_per_trial_scores_sp=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sp;
                this_mean_per_trial_scores_sp=this_mean_per_trial_scores_sp(:,:,1:length(time_span));

                these_scores_sm=zeros(size(this_mean_per_trial_scores_sp,1),size(this_mean_per_trial_scores_sp,3));
                these_scores_sm(:,:)=this_mean_per_trial_scores_sp(:,1,:);

                if size(these_scores_sm,1)>2

                    CIsm = bootci(1000, @mean, these_scores_sm);

                    meansm=mean(these_scores_sm,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;


                    [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);

                else
                    plot(time_span',these_scores_sm', '-','Color',[158/255 31/255 99/255]);
                end

                these_scores_sp=zeros(size(this_mean_per_trial_scores_sp,1),size(this_mean_per_trial_scores_sp,3));
                these_scores_sp(:,:)=this_mean_per_trial_scores_sp(:,2,:);

                if size(these_scores_sp,1)>2
                    CIsp = bootci(1000, @mean, these_scores_sp);
                    meansp=mean(these_scores_sp,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;

                    [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);
                else
                    plot(time_span',these_scores_sp', '-','Color',[158/255 31/255 99/255]);
                end


                if (size(these_scores_sp,1)>2)&(size(these_scores_sm,1)>2)
                    plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                    plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
                end

                text(50,0.75,'S-','Color',[158/255 31/255 99/255])
                text(50,0.65,'S+','Color',[0 114/255 178/255])

                ylim([0 1])
                title(['S+ ' handles.group_names{grNo}])
                xlabel('Time(sec)')
                ylabel('Posterior p')

                ii_plot=ii_plot+1;
                subplot(6,2,ii_plot)
                hold on

                this_mean_per_trial_scores_sm=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).mean_per_trial_scores_sm;
                this_mean_per_trial_scores_sm=this_mean_per_trial_scores_sm(:,:,1:length(time_span));

                these_scores_sm=zeros(size(this_mean_per_trial_scores_sm,1),size(this_mean_per_trial_scores_sm,3));
                these_scores_sm(:,:)=this_mean_per_trial_scores_sm(:,1,:);

                if size(these_scores_sm,1)>2
                    CIsm = bootci(1000, @mean, these_scores_sm);
                    meansm=mean(these_scores_sm,1);
                    CIsm(1,:)=meansm-CIsm(1,:);
                    CIsm(2,:)=CIsm(2,:)-meansm;

                    [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
                else
                    plot(time_span',these_scores_sm', '-','Color',[158/255 31/255 99/255]);
                end

                these_scores_sp=zeros(size(this_mean_per_trial_scores_sm,1),size(this_mean_per_trial_scores_sm,3));
                these_scores_sp(:,:)=this_mean_per_trial_scores_sm(:,2,:);

                if size(these_scores_sp,1)>2
                    CIsp = bootci(1000, @mean, these_scores_sp);
                    meansp=mean(these_scores_sp,1);
                    CIsp(1,:)=meansp-CIsp(1,:);
                    CIsp(2,:)=CIsp(2,:)-meansp;

                    [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);
                else
                    plot(time_span',these_scores_sp', '-','Color',[0 114/255 178/255]);
                end



                if (size(these_scores_sp,1)>2)&(size(these_scores_sm,1)>2)
                    plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                    plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
                end

                text(50,0.75,'S-','Color',[158/255 31/255 99/255])
                text(50,0.65,'S+','Color',[0 114/255 178/255])

                ylim([0 1])
                title(['S- ' handles.group_names{grNo}])
                xlabel('Time(sec)')
                ylabel('Posterior p')
            catch
            end
        end
    end
    sgtitle(['Posterior probability for ' handles_out2.classifier_names{iiMLalgo} ])
end



%Average decoding accuracy first and last session
for iiMLalgo=handles.MLalgo_to_use
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    length_corr=[];
    for mouseNo=1:size(first_last,1)
        length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict);
    end

    first_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
    for mouseNo=1:size(first_last,1)
        first_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict;
    end

    these_per_corr=first_mean_correct_predict;
    this_time_spanf=zeros(1,size(these_per_corr,2));
    this_time_spanf(1)=time_span(1);
    dt=time_span(2)-time_span(1);
    for ii=2:length(this_time_spanf)
        this_time_spanf(ii)=this_time_spanf(ii-1)+dt;
    end


    if size(these_per_corr,1)>2

        CIpv = bootci(1000, @mean, these_per_corr);
        meanpv=mean(these_per_corr,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(this_time_spanf',mean(these_per_corr,1)', CIpv', 'm');
    else
        if size(these_per_corr,1)>1
            plot(this_time_spanf',mean(these_per_corr,1)', 'm');
        else
            plot(this_time_spanf',these_per_corr', 'm');
        end

    end




    length_corr=[];
    for mouseNo=1:size(first_last,1)
        length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict);
    end

    last_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
    for mouseNo=1:size(first_last,1)
        last_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict;
    end

    these_per_corr=last_mean_correct_predict;
    this_time_span=zeros(1,size(these_per_corr,2));
    this_time_span(1)=time_span(1);
    for ii=2:length(this_time_span)
        this_time_span(ii)=this_time_span(ii-1)+dt;
    end

    if size(these_per_corr,1)>2

        CIpv = bootci(1000, @mean, these_per_corr);
        meanpv=mean(these_per_corr,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(this_time_span',mean(these_per_corr,1)', CIpv', 'b');
    else
        if size(these_per_corr,1)>1
            plot(this_time_span',mean(these_per_corr,1)', 'b');
        else
            plot(this_time_span',these_per_corr', 'b');
        end

    end

    text(20,0.8,'First','Color','m')
    text(20,0.9,'Last','Color','b')
    plot(this_time_spanf',mean(first_mean_correct_predict,1)', 'm');
    ylabel('Accuracy')
    xlabel('Time(sec)')
    title(['Accuracy ' handles_out2.classifier_names{iiMLalgo}])
end

%Decoding accuracy first and last session per mouse
for iiMLalgo=handles.MLalgo_to_use
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .7])


    length_corr=[];
    for mouseNo=1:size(first_last,1)
        length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict);
    end

    first_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
    for mouseNo=1:size(first_last,1)
        first_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).first_mean_correct_predict;
    end

    these_per_corr=first_mean_correct_predict;
    this_time_spanf=zeros(1,size(these_per_corr,2));
    this_time_spanf(1)=time_span(1);
    dt=time_span(2)-time_span(1);
    for ii=2:length(this_time_spanf)
        this_time_spanf(ii)=this_time_spanf(ii-1)+dt;
    end


    length_corr=[];
    for mouseNo=1:size(first_last,1)
        length_corr(mouseNo)=length(handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict);
    end

    last_mean_correct_predict=zeros(size(first_last,1),max(length_corr));
    for mouseNo=1:size(first_last,1)
        last_mean_correct_predict(mouseNo,1:length_corr(mouseNo))=handles_out2.mouse_no(mouseNo).ii_thr(ii_out).MLalgo(iiMLalgo).last_mean_correct_predict;
    end

    these_per_corr=last_mean_correct_predict;
    this_time_span=zeros(1,size(these_per_corr,2));
    this_time_span(1)=time_span(1);
    for ii=2:length(this_time_span)
        this_time_span(ii)=this_time_span(ii-1)+dt;
    end

    %Now plot them
    for mouseNo=1:size(first_last,1)
        subplot(size(first_last,1),1,mouseNo)
        hold on
        plot(this_time_spanf',first_mean_correct_predict(mouseNo,:),'-m')
        plot(this_time_span',last_mean_correct_predict(mouseNo,:),'-b')


        if mouseNo==1
            text(20,0.8,'First','Color','m')
            text(20,0.9,'Last','Color','b')
        end

        ylabel('Accuracy')
        xlabel('Time(sec)')
        title(['Accuracy ' handles_out2.mouse_names{mouseNo}])
    end
    
    sgtitle(['Accuracy ' handles_out2.classifier_names{iiMLalgo}])
end

out_file=[choiceBatchPathName choiceFileName];
out_file=[out_file(1:end-2) '.mat'];
save(out_file,'handles_out2','handles','-v7.3')
pfft=1;