%drgCaImAn_analyze_batch_pre_per_to_mulitROI_decode_per_mouse

close all
clear all


[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_decode_entire_session_multi_ROI_fsdz for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

fileID = fopen([choiceBatchPathName 'drgCaImAn_analyze_batch_pre_per_to_mulitROI_decode_per_mouse.txt'],'w');

all_no_ROIs=[1 2 5 15 2000];
all_no_ROI_draws=[2000 40 40  40 1];

%separate the clusters according to diversion times
cluster_times=[-1.5 -1;-1 0; 0 1; 1 2; 2 3; 3 200];
cluster_labels{1}='-1.5 to -1 sec';
cluster_labels{2}='-1 to 0 sec';
cluster_labels{3}='0 to 1 sec';
cluster_labels{4}='1-2 sec';
cluster_labels{5}='2-3 sec';
cluster_labels{6}='>3sec';

epoch_labels{1}='Pre';
epoch_labels{2}='Odor';

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;
timeEvents = [-1.5 0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];


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

first_file=1;

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

if all_files_present==1


    %Process each file separately
    fileNo_included=0;
    mouseNo_per_file_included=[];
    for fileNo=first_file:length(handles.FileName_pre_per)
        tic
        first_toc=toc;
        handles_out=[];
        ii_out=0;

        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};
        [percent_correct] = drgCaImAnFindPercentCorrect(pre_per_PathName, pre_per_FileName);

        if percent_correct>=80
            %Do only for proficient
            fileNo_included=fileNo_included+1;
            mouseNo_per_file_included(fileNo_included)=mouseNo_per_file(fileNo);

            for ii_ROI_choices=1:length(all_no_ROIs)


                handles_choices.ii_out=ii_ROI_choices;
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



                handles_out_all.file(fileNo_included).ii_out(ii_ROI_choices).handles_choices=handles_choices;

                start_toc=toc;

                handles_out_all.file(fileNo_included).ii_out(ii_ROI_choices).handles_outd=drgCaImAnInspectDecodingMultiROI(handles_choices);

                fprintf(1, ['Data processed for file number %d, ii_out= %d\n'],fileNo,ii_ROI_choices);
            end
            fprintf(1,'Processing time for file number %d is %d seconds\n',fileNo,(toc-first_toc));

        else
            fprintf(1,'File number %d was excluded because percent correct <80\n',fileNo);
        end
    end


    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));

    fprintf(1,'Number of files included %d\n',fileNo_included);

    %Plot accurcy histograms and do glm
    glm_acc=[];
    glm_acc_ii=0;
    input_acc_data=[];
    id_acc_ii=0;

    %Plot the histograms for accuracy for each number of ROIs in the odor
    %window
    figureNo=0;
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);


    set(hFig, 'units','normalized','position',[.3 .1 .25 .7])

    fprintf(1, ['\nAccuracies for odor window\n'])
    fprintf(fileID, ['\nAccuracies for odor window\n']);

    edges=[0:0.05:1];
    
    for ii_ROI_choices=1:length(all_no_ROIs)
        
        subplot(length(all_no_ROIs),1,ii_ROI_choices)
        ax=gca;ax.LineWidth=3;
        hold on

        all_accs=[];
        all_accs_sh=[];
        for fileNo=first_file:fileNo_included
            this_time_span=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.time_span;
            these_accs=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI;
            these_accs_sh=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI_sh;

            vetted_accs=[];
            ii_v=0;
            for ii_repeats=1:length(these_accs)
                this_acc_time_course=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy;
                %NOte: some of the noisy traces yield artifactual
                %accuracy above 0.5

                if mean(this_acc_time_course(this_time_span<-2))<0.53
                    ii_v=ii_v+1;
                    vetted_accs(ii_v)=these_accs(ii_repeats);
                    vetted_accs_sh(ii_v)=these_accs_sh(ii_repeats);
                end
            end
            if ~isempty(vetted_accs)
                all_accs=[all_accs vetted_accs];
                all_accs_sh=[all_accs_sh vetted_accs_sh];

                glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(vetted_accs))=vetted_accs;
                glm_acc.ROI_group(glm_acc_ii+1:glm_acc_ii+length(vetted_accs))=ii_ROI_choices*ones(1,length(vetted_accs));
                glm_acc.epoch(glm_acc_ii+1:glm_acc_ii+length(vetted_accs))=0*ones(1,length(vetted_accs));
                glm_acc_ii=glm_acc_ii+length(vetted_accs);

                id_acc_ii=id_acc_ii+1;
                input_acc_data(id_acc_ii).data=vetted_accs;
                if ii_ROI_choices<length(all_no_ROIs)
                    input_acc_data(id_acc_ii).description=[epoch_labels{1} ' ' num2str(all_no_ROIs(ii_ROI_choices)) ' ROIs'];
                else
                    input_acc_data(id_acc_ii).description=[epoch_labels{1} ' all ROIs'];
                end

            end

        end

        fprintf(1, ['\nPercent ROIs with accuracy above 0.65 %d, n=%d\n'], 100*sum(all_accs>0.65)/length(all_accs),length(all_accs))
        fprintf(fileID, ['\nPercent ROIs with accuracy above 0.65 %d, n=%d\n'], 100*sum(all_accs>0.65)/length(all_accs),length(all_accs));

        histogram(all_accs,edges,'Normalization','Probability')
        histogram(all_accs_sh,edges,'Normalization','Probability')

        if ii_ROI_choices~=length(all_no_ROIs)
            title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
        else
            title(['All ROIs'])
        end
        if ii_ROI_choices==1
            xlabel('Accuracy')
        end
        ylabel('P')
    end
    sgtitle('Histograms for accuracy in odor window')

    %Plot the histograms for accuracy for each number of ROIs in the pre
    %window
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);


     set(hFig, 'units','normalized','position',[.3 .1 .25 .7])

    edges=[0:0.05:1];
    for ii_ROI_choices=1:length(all_no_ROIs)
        subplot(length(all_no_ROIs),1,ii_ROI_choices)
        ax=gca;ax.LineWidth=3;
        hold on

        all_accs=[];
        all_accs_sh=[];
        for fileNo=first_file:fileNo_included

            this_time_span=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.time_span;
            these_accs=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI_pre;
            these_accs_sh=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI_sh_pre;

            vetted_accs=[];
            ii_v=0;
            for ii_repeats=1:length(these_accs)
                this_acc_time_course=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy;
                %NOte: some of the noisy traces yield artifactual
                %accuracy above 0.5

                if mean(this_acc_time_course(this_time_span<-2))<0.53
                    ii_v=ii_v+1;
                    vetted_accs(ii_v)=these_accs(ii_repeats);
                    vetted_accs_sh(ii_v)=these_accs_sh(ii_repeats);
                end
            end
            if ~isempty(vetted_accs)
                all_accs=[all_accs vetted_accs];
                all_accs_sh=[all_accs_sh vetted_accs_sh];

                glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(vetted_accs))=vetted_accs;
                glm_acc.ROI_group(glm_acc_ii+1:glm_acc_ii+length(vetted_accs))=ii_ROI_choices*ones(1,length(vetted_accs));
                glm_acc.epoch(glm_acc_ii+1:glm_acc_ii+length(vetted_accs))=1*ones(1,length(vetted_accs));
                glm_acc_ii=glm_acc_ii+length(vetted_accs);

                id_acc_ii=id_acc_ii+1;
                input_acc_data(id_acc_ii).data=vetted_accs;
                  if ii_ROI_choices<length(all_no_ROIs)
                    input_acc_data(id_acc_ii).description=[epoch_labels{2} ' ' num2str(all_no_ROIs(ii_ROI_choices)) ' ROIs'];
                else
                    input_acc_data(id_acc_ii).description=[epoch_labels{2} ' all ROIs'];
                end
            end

        end
        histogram(all_accs,edges,'Normalization','Probability')
        histogram(all_accs_sh,edges,'Normalization','Probability')

        if ii_ROI_choices~=length(all_no_ROIs)
            title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
        else
            title(['All ROIs'])
        end
        if ii_ROI_choices==1
            xlabel('Accuracy')
        end
        ylabel('P')
    end
    sgtitle('Histograms for accuracy in pre window')

    %Perform the glm
    fprintf(1, ['\nglm for decoding accuracy\n'])
    fprintf(fileID, ['\nglm for decoding accuracy\n']);

    tbl = table(glm_acc.data',glm_acc.ROI_group',glm_acc.epoch',...
        'VariableNames',{'accuracy','number_of_ROIs','time_window'});
    mdl = fitglm(tbl,'accuracy~number_of_ROIs+time_window+number_of_ROIs*time_window'...
        ,'CategoricalVars',[2,3])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);


    %Plot the histograms for latency for each number of ROIs
    %and do glm
       glm_lat=[];
    glm_lat_ii=0;
    input_lat_data=[];
    id_lat_ii=0;

    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    set(hFig, 'units','normalized','position',[.3 .1 .25 .7])

    edges=[-2:0.5:10];
    all_lats_per_ROI_choice=[];
    for ii_ROI_choices=1:length(all_no_ROIs)
        subplot(length(all_no_ROIs),1,ii_ROI_choices)
        ax=gca;ax.LineWidth=3;
        hold on

        all_lats=[];
        for fileNo=first_file:fileNo_included
            these_lats=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.latency_per_ROI;


            this_time_span=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.time_span;


            vetted_lats=[];
            ii_v=0;
            for ii_repeats=1:length(these_lats)
                this_acc_time_course=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy;
                %NOte: some of the noisy traces yield artifactual
                %accuracy above 0.5
                if (mean(this_acc_time_course(this_time_span<-2))<0.53)&(~isnan(these_lats(ii_repeats)))
                    ii_v=ii_v+1;
                    vetted_lats(ii_v)=these_lats(ii_repeats);
                end
            end
            if ~isempty(vetted_lats)
                all_lats=[all_lats vetted_lats];

                  glm_lat.data(glm_lat_ii+1:glm_lat_ii+length(vetted_lats))=vetted_lats;
                glm_lat.ROI_group(glm_lat_ii+1:glm_lat_ii+length(vetted_lats))=ii_ROI_choices*ones(1,length(vetted_lats));
                % glm_lat.epoch(glm_lat_ii+1:glm_lat_ii+length(vetted_lats))=0*ones(1,length(vetted_lats));
                glm_lat_ii=glm_lat_ii+length(vetted_lats);

                id_lat_ii=id_lat_ii+1;
                input_lat_data(id_lat_ii).data=vetted_lats;
                if ii_ROI_choices<length(all_no_ROIs)
                    input_lat_data(id_lat_ii).description=[num2str(all_no_ROIs(ii_ROI_choices)) ' ROIs'];
                else
                    input_lat_data(id_lat_ii).description=['all ROIs'];
                end

            end
        end
      
        histogram(all_lats,edges)


        if ii_ROI_choices~=length(all_no_ROIs)
            title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
        else
            title(['All ROIs'])
        end
        if ii_ROI_choices==1
            xlabel('Latency')
        end
        ylabel('P')
    end
    sgtitle('Histograms for latency')

    %Perform the glm
    fprintf(1, ['\nglm for latency\n'])
    fprintf(fileID, ['\nglm for latency\n']);

    tbl = table(glm_acc.data',glm_acc.ROI_group',...
        'VariableNames',{'accuracy','number_of_ROIs'});
    mdl = fitglm(tbl,'accuracy~number_of_ROIs'...
        ,'CategoricalVars',[2])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);

    %Plot the pseudocolor time courses
    time_span=handles_out_all.file(1).ii_out(1).handles_outd.time_span;
    for ii_ROI_choices=1:length(all_no_ROIs)
        % figureNo = figureNo + 1;
        % try
        %     close(figureNo)
        % catch
        % end
        % hFig=figure(figureNo);
        % 
        % 
        % set(hFig, 'units','normalized','position',[.3 .1 .25 .7])
        % 
        % ax=gca;ax.LineWidth=3;
        % hold on

        all_acc_time_courses=zeros(5000,length(time_span));
        all_lats_for_acc_time_courses=[];
        ii_included=0;
        for fileNo=first_file:fileNo_included
            these_accs=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI;
            these_lats=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.latency_per_ROI;
            for ii_repeats=1:length(these_accs)
                if these_accs(ii_repeats)>=handles_choices.acc_thr(2)
                    this_time_span=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.time_span;
                    this_acc_time_course=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy;
                    %NOte: some of the noisy traces yield artifactual
                    %accuracy above 0.5
                    if mean(this_acc_time_course(this_time_span<-2))<0.53
                        ii_included=ii_included+1;
                        all_lats_for_acc_time_courses(ii_included)=these_lats(ii_repeats);

                        resample=0;
                        if length(this_time_span)~=length(time_span)
                            resample=1;
                        else
                            if sum(time_span==this_time_span)~=length(time_span)
                                resample=1;
                            end
                        end
                        if resample==0
                            all_acc_time_courses(ii_included,:)= handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy;
                        else
                            %Resample
                            for ii_t=1:length(time_span)
                                if time_span(ii_t)<this_time_span(1)
                                    all_acc_time_courses(ii_included,ii_t)= handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy(1);
                                else
                                    if time_span(ii_t)>this_time_span(end)
                                        all_acc_time_courses(ii_included,ii_t)= handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy(end);
                                    else
                                        if isempty(find(time_span(ii_t)==this_time_span))
                                            ii_before=find(this_time_span<time_span(ii_t),1,'last');
                                            ii_after=find(this_time_span>time_span(ii_t),1,'first');
                                            acc_before=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy(ii_before);
                                            acc_after=handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy(ii_after);
                                            all_acc_time_courses(ii_included,ii_t)= acc_before+(acc_after-acc_before)*(time_span(ii_t)-this_time_span(ii_before))/(this_time_span(ii_after)-this_time_span(ii_before));
                                        else
                                            ii_found=find(time_span(ii_t)==this_time_span,1,'first');
                                            all_acc_time_courses(ii_included,ii_t)= handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.ROI(ii_repeats).mean_accuracy(ii_found);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        all_acc_time_courses(ii_included+1:end,:)=[];


        %Now sort according to latencies or through cross correlation
        sort_latency=1;

        if sort_latency==0
            croscorr_traces=corrcoef(all_acc_time_courses');

            %Set autocorrelations to zero
            %     for ii=1:size(croscorr_traces,1)
            %         croscorr_traces(ii,ii)=0;
            %     end
            Z = linkage(croscorr_traces,'complete','correlation');
            no_clusters=4; %Note that this is different
            clusters = cluster(Z,'Maxclust',no_clusters);
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);
            %Do cutoff for 4 clusters
            cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
            [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
            set(H,'LineWidth',2)
            hFig=figure(figureNo);
            set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

            %re-sort the matrix
            for ii=1:size(croscorr_traces,1)
                for jj=1:size(croscorr_traces,1)
                    perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
                end
            end




            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);

            set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
            hold on
            pcolor(perm_croscorr_traces)
            colormap fire
            shading flat

            caxis([-1    1])
            title(['Cross correlations for all ROIs'])
        else
            %Sort by diverging time
            to_sort=[all_lats_for_acc_time_courses' [1:ii_included]'];
            sorted=sortrows(to_sort);
            outperm=sorted(:,2);
            
           

            clusters=[];
            for ii_tc=1:length(all_lats_for_acc_time_courses)
                for clus=1:size(cluster_times,1)
                    if (all_lats_for_acc_time_courses(ii_tc)>=cluster_times(clus,1))&(all_lats_for_acc_time_courses(ii_tc)<=cluster_times(clus,2))
                        clusters(ii_tc)=clus;
                    end
                end
            end

        end

        %Plot timecourses for all accuracy time courses
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
        hold on

        sorted_handles_out.all_acc_time_courses=[];
        sorted_handles_out.clusters=[];
        for ii=1:ii_included
            sorted_handles_out.all_acc_time_courses(ii_included-ii+1,:)=all_acc_time_courses(outperm(ii),:);
            sorted_handles_out.clusters(ii_included-ii+1)=clusters(outperm(ii));
        end

        %pcolor does not show the first row
        pseudo_acc=zeros(size(sorted_handles_out.all_acc_time_courses,1)+1,size(sorted_handles_out.all_acc_time_courses,2));
        pseudo_acc(1:end-1,:)=sorted_handles_out.all_acc_time_courses;

        time_span_mat=repmat(time_span,ii_included+1,1);
        ROI_mat=repmat(1:ii_included+1,length(time_span),1)';

        pcolor(time_span_mat,ROI_mat,pseudo_acc)


        %         time_span_mat=repmat(time_span,sorted_handles_out.all_spmresp_ii_dFF,1);
        %         ROI_mat=repmat(1:sorted_handles_out.all_spmresp_ii_dFF,length(time_span),1)';
        %
        %         pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_sporsm_dFFspm)

        colormap fire
        shading flat

        % caxis([prctile(sorted_handles_out.all_acc_time_courses(:),1) prctile(sorted_handles_out.all_acc_time_courses(:),99)])
        caxis([0.4,1])

        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],[0 ii_included],'-r')
        end

        xlim([-5 15])
        ylim([1 ii_included])
        if ii_ROI_choices~=length(all_no_ROIs)
            title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
        else
            title(['All ROIs'])
        end
        xlabel('Time (sec)')
        ylabel('ROI number')

        % %Plot rainbow
        % figureNo=figureNo+1;
        % try
        %     close(figureNo)
        % catch
        % end
        % 
        % hFig = figure(figureNo);
        % 
        % set(hFig, 'units','normalized','position',[.49 .1 .05 .3])
        % 
        % prain=[prctile(sorted_handles_out.all_acc_time_courses(:),1):(prctile(sorted_handles_out.all_acc_time_courses(:),99)-prctile(sorted_handles_out.all_acc_time_courses(:),1))/99:prctile(sorted_handles_out.all_acc_time_courses(:),99)];
        % pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
        % %             colormap jet
        % colormap fire
        % shading interp
        % ax=gca;
        % set(ax,'XTickLabel','')

        %Plot the average timecourses per cluster
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.25 .35 .25 .2])
        hold on

        for clus=max(clusters):-1:1

            % subplot(2,3, clus)
            % ax=gca;ax.LineWidth=3;
            % hold on



            %plot the accuracies
            this_cluster_acc=[];
            ii_inc=0;
            for ii=1:length(sorted_handles_out.clusters)
                if sorted_handles_out.clusters(ii)==clus
                    ii_inc=ii_inc+1;
                    this_cluster_acc(ii_inc,:)=sorted_handles_out.all_acc_time_courses(ii,:);
                end
            end

            if ii_inc>=5
                CIpv = bootci(1000, @mean, this_cluster_acc);
                meanpv=mean(this_cluster_acc,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;


                [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_acc), CIpv','cmap',[0 0 0]);


            end

        end

        ylim([0.4 0.9])
        this_ylim=ylim;
        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
        end

        xlim([-7 15])


        xlabel('Time(sec)')
        ylabel('Accuracy')


        % title(['Cluster from ' num2str(cluster_times(clus,1)) ' to ' num2str(cluster_times(clus,2)) ' sec'])
        %
        %
        %
        %
        if ii_ROI_choices~=length(all_no_ROIs)
            title(['Accuracy per cluster No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
        else
            title(['Accuracy per cluster All ROIs'])
        end

        %Plot average of all traces
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.25 .35 .25 .2])
        hold on



        % subplot(2,3, clus)
        % ax=gca;ax.LineWidth=3;
        % hold on



        % %plot the accuracies
        % this_cluster_acc=[];
        % ii_inc=0;
        % for ii=1:length(sorted_handles_out.clusters)
        %     if sorted_handles_out.clusters(ii)==clus
        %         ii_inc=ii_inc+1;
        %         this_cluster_acc(ii_inc,:)=sorted_handles_out.all_acc_time_courses(ii,:);
        %     end
        % end


        CIpv = bootci(1000, @mean, sorted_handles_out.all_acc_time_courses);
        meanpv=mean(sorted_handles_out.all_acc_time_courses,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(time_span,mean(sorted_handles_out.all_acc_time_courses), CIpv','cmap',[0 0 0]);



        ylim([0.4 0.9])
        this_ylim=ylim;
        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
        end

        xlim([-7 15])


        xlabel('Time(sec)')
        ylabel('Accuracy')


        % title(['Cluster from ' num2str(cluster_times(clus,1)) ' to ' num2str(cluster_times(clus,2)) ' sec'])
        %
        %
        %
        %
        if ii_ROI_choices~=length(all_no_ROIs)
            title(['Accuracy per cluster No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
        else
            title(['Accuracy per cluster All ROIs'])
        end
        pffft=1;
    end

    %Plot rainbow
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])

    % prain=[prctile(sorted_handles_out.all_acc_time_courses(:),1):(prctile(sorted_handles_out.all_acc_time_courses(:),99)-prctile(sorted_handles_out.all_acc_time_courses(:),1))/99:prctile(sorted_handles_out.all_acc_time_courses(:),99)];

    prain=[0.4:0.6/99:1];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    % %Now generate per mouse histograms
    % for ii_mouse=1:no_mice
    % 
    %     %Plot the histograms for accuracy for each number of ROIs in the odor
    %     %window
    %     figureNo = figureNo + 1;
    %     try
    %         close(figureNo)
    %     catch
    %     end
    %     hFig=figure(figureNo);
    % 
    % 
    %      set(hFig, 'units','normalized','position',[.3 .1 .25 .7])
    % 
    %     edges=[0:0.05:1];
    %     for ii_ROI_choices=1:length(all_no_ROIs)
    %         subplot(length(all_no_ROIs),1,ii_ROI_choices)
    %         ax=gca;ax.LineWidth=3;
    %         hold on
    % 
    %         all_accs=[];
    %         all_accs_sh=[];
    %         for fileNo=first_file:fileNo_included
    %             if mouseNo_per_file_included(fileNo)==ii_mouse
    %                 all_accs=[all_accs handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI];
    %                 all_accs_sh=[all_accs_sh handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI_sh];
    %             end
    %         end
    %         histogram(all_accs,edges,'Normalization','Probability')
    %         histogram(all_accs_sh,edges,'Normalization','Probability')
    % 
    %         if ii_ROI_choices~=length(all_no_ROIs)
    %             title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
    %         else
    %             title(['All ROIs'])
    %         end
    %         if ii_ROI_choices==1
    %             xlabel('Accuracy')
    %         end
    %         ylabel('P')
    %     end
    %     sgtitle(['Histograms for accuracy odor window for mouse No ' num2str(ii_mouse)])
    % 
    %     %Plot the histograms for accuracy for each number of ROIs in the pre
    %     %window
    %     figureNo = figureNo + 1;
    %     try
    %         close(figureNo)
    %     catch
    %     end
    %     hFig=figure(figureNo);
    % 
    % 
    %      set(hFig, 'units','normalized','position',[.3 .1 .25 .7])
    % 
    %     edges=[0:0.05:1];
    %     for ii_ROI_choices=1:length(all_no_ROIs)
    %         subplot(length(all_no_ROIs),1,ii_ROI_choices)
    %         ax=gca;ax.LineWidth=3;
    %         hold on
    % 
    %         all_accs=[];
    %         all_accs_sh=[];
    %         for fileNo=first_file:fileNo_included
    %             if mouseNo_per_file_included(fileNo)==ii_mouse
    %                 all_accs=[all_accs handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI_pre];
    %                 all_accs_sh=[all_accs_sh handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.accuracy_per_ROI_sh_pre];
    %             end
    %         end
    %         histogram(all_accs,edges,'Normalization','Probability')
    %         histogram(all_accs_sh,edges,'Normalization','Probability')
    % 
    %         if ii_ROI_choices~=length(all_no_ROIs)
    %             title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
    %         else
    %             title(['All ROIs'])
    %         end
    %         if ii_ROI_choices==1
    %             xlabel('Accuracy')
    %         end
    %         ylabel('P')
    %     end
    %     sgtitle(['Histograms for accuracy pre window for mouse No ' num2str(ii_mouse)])
    % 
    %     %Plot the histograms for latency for each number of ROIs
    %     figureNo = figureNo + 1;
    %     try
    %         close(figureNo)
    %     catch
    %     end
    %     hFig=figure(figureNo);
    % 
    % 
    %      set(hFig, 'units','normalized','position',[.3 .1 .25 .7])
    % 
    %     edges=[-2:0.5:10];
    %     for ii_ROI_choices=1:length(all_no_ROIs)
    %         subplot(length(all_no_ROIs),1,ii_ROI_choices)
    %         ax=gca;ax.LineWidth=3;
    %         hold on
    % 
    %         all_lats=[];
    %         for fileNo=first_file:fileNo_included
    %             if mouseNo_per_file_included(fileNo)==ii_mouse
    %             all_lats=[all_lats handles_out_all.file(fileNo).ii_out(ii_ROI_choices).handles_outd.latency_per_ROI];
    %             end
    %         end
    %         histogram(all_lats,edges)
    % 
    % 
    %         if ii_ROI_choices~=length(all_no_ROIs)
    %             title(['No ROIs =' num2str(all_no_ROIs(ii_ROI_choices))])
    %         else
    %             title(['All ROIs'])
    %         end
    %         if ii_ROI_choices==1
    %             xlabel('Latency')
    %         end
    %         ylabel('P')
    %     end
    %     sgtitle(['Histograms for latency for mouse No ' num2str(ii_mouse)])
    % 
    % end
end 

fclose(fileID)
pfft=1;