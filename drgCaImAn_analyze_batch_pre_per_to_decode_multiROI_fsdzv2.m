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

%Process each file separately
 for grNo=1:max(handles.group)
%     handles_out2.group_no(grNo).ii_euclid=0;
   for iiMLalgo=handles.MLalgo_to_use
        for ii_out=1:length(handles.p_threshold)
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=0;
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_other_ROIs=[];
        end
    end
 end
 
for fileNo=1:no_files
    tic
    pre_per_PathName=handles.PathName_pre_per{fileNo};
    pre_per_FileName=handles.FileName_pre_per{fileNo};
    grNo=handles.group(fileNo);

    load([handles.PathName_out{fileNo} pre_per_FileName(1:end-4) handles.suffix_out])
    for iiMLalgo=handles.MLalgo_to_use
%         if iiMLalgo==handles.MLalgo_to_use(1)
%             handles_out2.group_no(grNo).ii_euclid=handles_out2.group_no(grNo).ii_euclid+1;
%             ii_euclid=handles_out2.group_no(grNo).ii_euclid;
%             handles_out2.group_no(grNo).dist_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.dist_euclid))=handles_out.ii_out(1).handles_out.dist_euclid-handles_out.ii_out(1).handles_out.dist_euclid_zero;
%             handles_out2.group_no(grNo).KLdivergence(ii_euclid,1:length(handles_out.ii_out(1).handles_out.KLdivergence))=handles_out.ii_out(1).handles_out.KLdivergence;
%             handles_out2.group_no(grNo).time_span_euclid(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.time_span;
%             handles_out2.group_no(grNo).ii_time_span(ii_euclid,1)=length(handles_out.ii_out(1).handles_out.time_span);
%             handles_out2.group_no(grNo).meandFFsp(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsp;
%             handles_out2.group_no(grNo).meandFFsm(ii_euclid,1:length(handles_out.ii_out(1).handles_out.time_span))=handles_out.ii_out(1).handles_out.meandFFsm;
%         end
        for ii_out=1:length(handles_out.ii_out)
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii+1;
            ii=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).ii;
            accuracy_tr=handles_out.ii_out(ii_out).handles_out.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr(ii)=accuracy_tr;
            accuracy_tr_other_ROIs=[];
            no_ROI_draws=handles_out.ii_out(ii_out).handles_choices.no_ROI_draws;
            for iiROI=2:no_ROI_draws+1
                accuracy_tr_other_ROIs=[accuracy_tr_other_ROIs handles_out.ii_out(ii_out).handles_out.ROI(iiROI).MLalgo(iiMLalgo).accuracy_tr];
            end
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_other_ROIs=[handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).accuracy_tr_other_ROIs accuracy_tr_other_ROIs];
            %Normalize the distribution
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).norm_accuracy_tr(ii)=(accuracy_tr-mean(accuracy_tr_other_ROIs))/std(accuracy_tr_other_ROIs);
            handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).norm_accuracy_tr_other_ROIs=(accuracy_tr_other_ROIs-mean(accuracy_tr_other_ROIs))/std(accuracy_tr_other_ROIs);
        end
    end
    fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
end

%Bar graph plot for accuracy
figureNo=0;
for iiMLalgo=handles.MLalgo_to_use


    all_accuracy_tr=[];
    all_accuracy_tr_other_ROIs=[];
    for grNo=1:max(handles.group)

        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.3 .3 .5 .25])

        hold on
        accuracy_tr=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).norm_accuracy_tr;
        all_accuracy_tr=[all_accuracy_tr accuracy_tr];
        accuracy_tr_other_ROIs=handles_out2.group_no(grNo).ii_thr(ii_out).MLalgo(iiMLalgo).norm_accuracy_tr_other_ROIs;
        all_accuracy_tr_other_ROIs=[all_accuracy_tr_other_ROIs accuracy_tr_other_ROIs];
        edges=[-3:0.1:3];
        histogram(accuracy_tr_other_ROIs,edges,'Normalization','probability')
        histogram(accuracy_tr,edges,'Normalization','probability')

        title(['Histogram of z scored accuracy for ' handles_out2.classifier_names{iiMLalgo} ' ' handles.group_names{grNo}])
        xlabel('z-score accuracy')

    end

    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .5 .25])

    hold on

   
    histogram(all_accuracy_tr_other_ROIs,edges,'Normalization','probability')
    histogram(all_accuracy_tr,edges,'Normalization','probability')

    title(['Histogram of z scored accuracy for ' handles_out2.classifier_names{iiMLalgo} ' all groups '])
    xlabel('z-scored accuracy')

end

pffft=1;
