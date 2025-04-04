%drgCaImAn_analyze_batch_pre_per_to_selectROI_decode_per_mouse

close all
clear all


[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_decode_entire_session_multi_ROI_fsdz for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

all_no_ROIs=[2000 2000];
all_no_ROI_draws=[1 1];


new_no_files=handles.no_files;

if isfield(handles,'processing_algo')
    processing_algo=handles.processing_algo;
else
    processing_algo=1;
end

if isfield(handles,'suffix_out')
    suffix_out=handles.suffix_out;
else
    suffix_out='_dec.mat';
end

if isfield(handles,'first_file')
    first_file=handles.first_file;
else
    first_file=1;
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

%Find the mouse numbers
mouseNo_per_file=[];
mouse_names=[];

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

no_mice=length(mouse_names);

figNo=0;
show_figures=1;
handles_out_all=[];
handles_out_all_rdec=[];

if all_files_present==1


    %Process each file separately
    fileNo_included=0;
    for fileNo=first_file:length(handles.FileName_pre_per)
        tic
        first_toc=toc;
        handles_out=[];
        ii_out=0;

        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};
        [percent_correct] = drgCaImAnFindPercentCorrect(pre_per_PathName, pre_per_FileName);

        if percent_correct>=80
            fileNo_included=fileNo_included+1;
            mouseNo_per_file_included(fileNo_included)=mouseNo_per_file(fileNo);

            %First bring in the rdec data for single ROI decoding
            handles_choices.ii_out=1;
            handles_choices.show_figures=0;
            handles_choices.pre_perFileName=[pre_per_FileName(1:end-4) '_rdec.mat'];
            handles_choices.pre_perPathName=pre_per_PathName;
            handles_choices.time_windows=[3.1 4.1];
            handles_choices.time_windows_pre=[-1 0];
            handles_choices.time_window_lat=[-1.5 10];
            handles_choices.pre_time_window=[-7 -1.5];
            handles_choices.acc_thr=[0.35 0.65];
            handles_choices.dt_span=15;
            handles_choices.MLalgo=6;

            handles_out_all_rdec.file(fileNo_included).handles_outd=drgCaImAnInspectDecodingMultiROI(handles_choices);

            %Do only for proficient
            for ii_process_low=1:2


                handles_choices.ii_out=ii_process_low;
                handles_choices.show_figures=0;
                handles_choices.pre_perFileName=[pre_per_FileName(1:end-4) '_sdec.mat'];
                handles_choices.pre_perPathName=pre_per_PathName;
                handles_choices.time_windows=[3.1 4.1];
                handles_choices.time_windows_pre=[-1 0];
                handles_choices.time_window_lat=[-1.5 10];
                handles_choices.pre_time_window=[-7 -1.5];
                handles_choices.acc_thr=[0.35 0.65];
                handles_choices.dt_span=15;
                handles_choices.MLalgo=6;



                handles_out_all.file(fileNo_included).ii_out(ii_process_low).handles_choices=handles_choices;

                start_toc=toc;

                handles_out_all.file(fileNo_included).ii_out(ii_process_low).handles_outd=drgCaImAnInspectDecodingSelectROI(handles_choices);

                fprintf(1, ['Data processed for file number %d, ii_out= %d\n'],fileNo,ii_process_low);
            end
            fprintf(1,'Processing time for file number %d is %d seconds\n',fileNo,(toc-first_toc));
        else
            fprintf(1,'File number %d was excluded because percent correct <80\n',fileNo);
        end
    end

    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));


    %Plot the histograms for accuracy for each number of ROIs in the odor
    %window
    figureNo=0;
    for ii_process_low=1:length(all_no_ROIs)
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);


        set(hFig, 'units','normalized','position',[.2 .3 .4 .4])

        edges=[0:0.05:1];

        %plot the accuracy determined using all ROIs per session
        subplot(2,1,2)
        ax=gca;ax.LineWidth=3;
        hold on

        all_accs=[];
        all_accs_sh=[];
        for fileNo=1:fileNo_included
            all_accs=[all_accs handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI];
            all_accs_sh=[all_accs_sh handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI_sh];
        end
        histogram(all_accs,edges,'Normalization','Probability')
        histogram(all_accs_sh,edges,'Normalization','Probability')


        if ii_process_low==1
            xlabel('Accuracy')
        end
        ylabel('P')

        title('Accuracy estimated with all ROIs per session')


        %plot the per ROI accuracy
        subplot(2,1,1)
        ax=gca;ax.LineWidth=3;
        hold on

        all_accs=[];
        all_accs_sh=[];
        for fileNo=1:fileNo_included
            these_accs=handles_out_all_rdec.file(fileNo).handles_outd.accuracy_per_ROI;
            these_accs_sh=handles_out_all_rdec.file(fileNo).handles_outd.accuracy_per_ROI_sh;
            if ii_process_low==1
                pruned_these_accs=these_accs(these_accs>=0.65);
                pruned_these_accs_sh=these_accs_sh(these_accs>=0.65);
            else
                pruned_these_accs=these_accs(these_accs<=0.35);
                pruned_these_accs_sh=these_accs_sh(these_accs<=0.35);
            end
            all_accs=[all_accs pruned_these_accs];
            all_accs_sh=[all_accs_sh pruned_these_accs_sh];
        end
        histogram(all_accs,edges,'Normalization','Probability')
        histogram(all_accs_sh,edges,'Normalization','Probability')


        if ii_process_low==1
            xlabel('Accuracy')
        end
        ylabel('P')

        title('Accuracy estimated per ROI')

        if ii_process_low==1
            sgtitle(['Histogram for odor decoding accuracy (single ROI accuracy >=0.65)'])
        else
            sgtitle(['Histogram for odor decoding accuracy (single ROI accuracy <=0.35)'])
        end
    end


    %Plot the histograms for accuracy for each number of ROIs in the pre
    %window
    for ii_process_low=1:length(all_no_ROIs)
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);


        set(hFig, 'units','normalized','position',[.3 .3 .4 .4])

        edges=[0:0.05:1];

        %plot the per ROI accuracy
        subplot(2,1,1)
        ax=gca;ax.LineWidth=3;
        hold on

        all_accs=[];
        all_accs_sh=[];
        for fileNo=1:fileNo_included
            these_accs=handles_out_all_rdec.file(fileNo).handles_outd.accuracy_per_ROI_pre;
            these_accs_sh=handles_out_all_rdec.file(fileNo).handles_outd.accuracy_per_ROI_sh_pre;
            if ii_process_low==1
                pruned_these_accs=these_accs(these_accs>=0.65);
                pruned_these_accs_sh=these_accs_sh(these_accs>=0.65);
            else
                pruned_these_accs=these_accs(these_accs<=0.35);
                pruned_these_accs_sh=these_accs_sh(these_accs<=0.35);
            end
            all_accs=[all_accs pruned_these_accs];
            all_accs_sh=[all_accs_sh pruned_these_accs_sh];
        end
        histogram(all_accs,edges,'Normalization','Probability')
        histogram(all_accs_sh,edges,'Normalization','Probability')


        if ii_process_low==1
            xlabel('Accuracy')
        end
        ylabel('P')

        title('Accuracy estimated per ROI')

        %Plot accuracy calculated for all ROIs per session
        subplot(2,1,2)
        ax=gca;ax.LineWidth=3;
        hold on

        all_accs=[];
        all_accs_sh=[];
        for fileNo=1:fileNo_included
            all_accs=[all_accs handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI_pre];
            all_accs_sh=[all_accs_sh handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI_sh_pre];
        end
        histogram(all_accs,edges,'Normalization','Probability')
        histogram(all_accs_sh,edges,'Normalization','Probability')


        if ii_process_low==1
            xlabel('Accuracy')
        end
        ylabel('P')

        title('Accuracy estimated with all ROIs per session')
        if ii_process_low==1
            sgtitle(['Histogram for pre decoding accuracy (single ROI accuracy >=0.65)'])
        else
            sgtitle(['Histogram for pre decoding accuracy (single ROI accuracy <=0.35)'])
        end

    end



    %Plot the histograms for latency
    for ii_process_low=1:length(all_no_ROIs)
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);


        set(hFig, 'units','normalized','position',[.2 .3 .4 .4])

        edges=[-2:0.5:10];


        %plot the per ROI latency
        subplot(2,1,1)
        ax=gca;ax.LineWidth=3;
        hold on

        all_lats=[];

        for fileNo=1:fileNo_included
            these_lats=handles_out_all_rdec.file(fileNo).handles_outd.latency_per_ROI;
            these_accs=handles_out_all_rdec.file(fileNo).handles_outd.accuracy_per_ROI;
            if ii_process_low==1
                pruned_these_lats=these_lats(these_accs>=0.65);
            else
                pruned_these_lats=these_lats(these_accs<=0.35);
            end
            all_lats=[all_lats pruned_these_lats];
        end
        histogram(all_lats,edges,'Normalization','Probability')

        if ii_process_low==1
            xlabel('Accuracy')
        end
        ylabel('P')

        title('Latency estimated per ROI')

        %Now plot latencies estimated for all ROIs per session
        subplot(2,1,2)
        ax=gca;ax.LineWidth=3;
        hold on

        all_lats=[];
        for fileNo=1:fileNo_included
            all_lats=[all_lats handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.latency_per_ROI];
        end
        histogram(all_lats,edges)



        if ii_process_low==1
            xlabel('Latency')
        end
        ylabel('P')

        title('Latency estimated for all ROIs per session')

        if ii_process_low==1
            sgtitle(['Histogram for latencies (single ROI accuracy >=0.65)'])
        else
            sgtitle(['Histogram for latencies (single ROI accuracy <=0.35)'])
        end
    end
    

    %Now generate per mouse histograms
    for ii_mouse=1:no_mice

        %Plot the histograms for accuracy for each number of ROIs in the odor
        %window
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);


        set(hFig, 'units','normalized','position',[.2 .3 .4 .4])

        edges=[0:0.05:1];
        for ii_process_low=1:length(all_no_ROIs)
            subplot(length(all_no_ROIs),1,ii_process_low)
            ax=gca;ax.LineWidth=3;
            hold on

            all_accs=[];
            all_accs_sh=[];
            for fileNo=1:fileNo_included
                if mouseNo_per_file_included(fileNo)==ii_mouse
                    all_accs=[all_accs handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI];
                    all_accs_sh=[all_accs_sh handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI_sh];
                end
            end
            histogram(all_accs,edges,'Normalization','Probability')
            histogram(all_accs_sh,edges,'Normalization','Probability')

            if ii_process_low==1
                title(['ROIs with decoding accuracy >=0.65'])
            else
                title(['ROIs with decoding accuracy <=0.35'])
            end
            if ii_process_low==1
                xlabel('Accuracy')
            end
            ylabel('P')
        end
        sgtitle(['Histograms for accuracy odor window for mouse No ' num2str(ii_mouse)])

        %Plot the histograms for accuracy for each number of ROIs in the pre
        %window
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);


        set(hFig, 'units','normalized','position',[.2 .3 .4 .4])

        edges=[0:0.05:1];
        for ii_process_low=1:length(all_no_ROIs)
            subplot(length(all_no_ROIs),1,ii_process_low)
            ax=gca;ax.LineWidth=3;
            hold on

            all_accs=[];
            all_accs_sh=[];
            for fileNo=1:fileNo_included
                if mouseNo_per_file_included(fileNo)==ii_mouse
                    all_accs=[all_accs handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI_pre];
                    all_accs_sh=[all_accs_sh handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.accuracy_per_ROI_sh_pre];
                end
            end
            histogram(all_accs,edges,'Normalization','Probability')
            histogram(all_accs_sh,edges,'Normalization','Probability')

            if ii_process_low==1
                title(['ROIs with decoding accuracy >=0.65'])
            else
                title(['ROIs with decoding accuracy <=0.35'])
            end
            if ii_process_low==1
                xlabel('Accuracy')
            end
            ylabel('P')
        end
        sgtitle(['Histograms for accuracy pre window for mouse No ' num2str(ii_mouse)])

        %Plot the histograms for latency for each number of ROIs
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);


        set(hFig, 'units','normalized','position',[.2 .3 .4 .4])

        edges=[-2:0.5:10];
        for ii_process_low=1:length(all_no_ROIs)
            subplot(length(all_no_ROIs),1,ii_process_low)
            ax=gca;ax.LineWidth=3;
            hold on

            all_lats=[];
            for fileNo=1:fileNo_included
                if mouseNo_per_file_included(fileNo)==ii_mouse
                all_lats=[all_lats handles_out_all.file(fileNo).ii_out(ii_process_low).handles_outd.latency_per_ROI];
                end
            end
            histogram(all_lats,edges)


            if ii_process_low==1
                title(['ROIs with decoding accuracy >=0.65'])
            else
                title(['ROIs with decoding accuracy <=0.35'])
            end
            if ii_process_low==1
                xlabel('Latency')
            end
            ylabel('P')
        end
        sgtitle(['Histograms for latency for mouse No ' num2str(ii_mouse)])

    end
end
pfft=1;