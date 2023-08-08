
%drgCaImAn_analyze_pval_batch.m
%Displays the data generated with drgCaImAn_pval_batch
clear all
close all


[outFileName,outhPathName] = uigetfile({'drgCaImAn_pval_choices*.mat'},'Select the .mat file with drgCaImAn_pval_batch');
load([outhPathName outFileName])

if ~isfield(handles_out,'use_pFDR')
    handles_out.use_pFDR=0;
end
time_span=handles_out.time_span;

which_task=0; %0=spm, 1=passive Ming
if which_task==1
    perCorr=50*ones(1,size(handles_out.perCorr,2));
    handles_out.perCorr=perCorr;
end

s=rng;

do_std=0;
no_clusters=3;

fr_per_names{1}='Fwd <40%';
fr_per_names{2}='Fwd 40-65%';
fr_per_names{3}='Fwd 65-80%';
fr_per_names{4}='Fwd >=80%';
fr_per_names{5}='Rev <40%';
fr_per_names{6}='Rev 40-65%';
fr_per_names{7}='Rev 65-80%';
fr_per_names{8}='Rev >=80%';

per_names{2}='40-65%';
per_names{3}='65-80%';
per_names{4}='>=80%';

naive_pro{1}='Naive';
naive_pro{3}='Proficient';

cluster_color{1}='m';
cluster_color{2}='b';
cluster_color{3}='y';
cluster_color{4}='r';
cluster_color{5}='c';
cluster_color{6}='g';

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

pre_t=[-3 -2];
fv_t=[-1 0];
odor_t=[2 4];

if isfield(handles_in,'min_tr_div')
    min_tr_div=handles_in.min_tr_div;
    min_tr_resp=handles_in.min_tr_resp;
else
    min_tr_div=12;
    min_tr_resp=6;
end

fileID = fopen([outhPathName 'drgCaImAn_analyze_pval_batch.txt'],'w');

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [-1.5 0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

%Time tinterval for shifting time base due to slow olfactometer computer
t_shift=0;
figureNo=0;


%Plot a bar graph showing percent signifcant divergence
glm_sig=[];
glm_sig_ii=0;
input_sig_data=[];
id_sig_ii=0;

%Find out which files were included (>=min_tr trials)
files_included_glm=zeros(1,length(handles_out.file));
for fNo=1:length(handles_out.file)
    if handles_out.file(fNo).output_data_odor.total_trials(1)>=min_tr_div
        files_included_glm(fNo)=1;
    end
end
files_included_glm=logical(files_included_glm);

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .2 .3])

hold on
bar_offset=0;
edges=[0:10:80];
rand_offset=0.5;

switch handles.group_algo
    case 1
        %Ming
        groups=[1 3];
    case 2
        %Fabio
        groups=unique(handles.group);
end
 

total_ROIs_per_group_per_mouse=zeros(max(groups),max((handles_out.mouseNo)));
for grNo=groups
    these_pre=[];
    these_FV=[];
    these_odor=[];
    
    for mouseNo=unique(handles_out.mouseNo)

        switch handles.group_algo
            case 1
                %Ming
                switch grNo
                    case 1
                        files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                    case 2
                        files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                    case 3
                        files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                end
            case 2
                %Fabio
                files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo)&files_included_glm;
        end
 
        switch handles_out.use_pFDR
            case {0,1}
                if sum(files_included)>0
                    n_pre=zeros(1,sum(files_included));
                    n_all=zeros(1,sum(files_included));
                    n_pre(1,:)=handles_out.output_data_pre(files_included);
                    n_all(1,:)=handles_out.output_data_N(files_included);
                    this_pre=n_pre./n_all;
                    these_pre=[these_pre 100*mean(this_pre)];

                    n_FV=zeros(1,sum(files_included));
                    n_FV(1,:)=handles_out.output_data_pre(files_included);
                    this_FV=n_FV./n_all;
                    these_FV=[these_FV 100*mean(this_FV)];

                    n_odor=zeros(1,sum(files_included));
                    n_odor(1,:)=handles_out.output_data_odor(files_included);
                    this_odor=n_odor./n_all;
                    these_odor=[these_odor 100*mean(this_odor)];
                end
            case {3,4}
                these_per_div=[];
                for fileNo=1:length(files_included)
                    if files_included(fileNo)==1
                        this_per_div=100*sum(handles_out.file(fileNo).output_data_odor.sig_div_glm)/length(handles_out.file(1).output_data_odor.sig_div_glm);

                        these_per_div=[these_per_div this_per_div];

                        total_ROIs_per_group_per_mouse(grNo,mouseNo)=total_ROIs_per_group_per_mouse(grNo,mouseNo)+sum(handles_out.file(fileNo).output_data_odor.sig_div_glm);

                    end
                end
                if ~isempty(these_per_div)
                    these_odor=[these_odor mean(these_per_div)];
                end
        end

    end 

    if ~isempty(these_odor)
        %             bar_offset=bar_offset+1;
        if grNo==1
            bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        else
            bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
        end
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',7);
        end

        glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
        glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
        glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
        glm_sig_ii=glm_sig_ii+length(these_odor);

        id_sig_ii=id_sig_ii+1;
        input_sig_data(id_sig_ii).data=these_odor;
        switch handles.group_algo
            case 1
                input_sig_data(id_sig_ii).description=naive_pro{grNo};
            case 2
                input_sig_data(id_sig_ii).description=naive_pro{grNo};
        end
    end

    bar_offset=bar_offset+1;
end

switch handles.group_algo
    case 1
        %Ming

        xticks([0 1 2])
        labels='xticklabels({';
        for ii_label=[1 3]
            labels=[labels '''' naive_pro{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
        xlim([-0.7 1.7])
        ylim([0 60])

    case 2
        %Fabio
        xticks(2*groups-1)
        labels='xticklabels({';
        for ii_label=1:length(groups)
            labels=[labels '''' handles.group_names{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
end


ylabel('% divergent')
title('Percent divergent ROIs per mouse')

% %Perform the glm
% fprintf(1, ['\nglm for percent divergent per mouse\n'])
% 
% 
% tbl = table(glm_sig.data',glm_sig.grNo',...
%     'VariableNames',{'percent_divergent','percent_correct'});
% mdl = fitglm(tbl,'percent_divergent~percent_correct'...
%     ,'CategoricalVars',[2])
% 
% 
% txt = evalc('mdl');
% txt=regexp(txt,'<strong>','split');
% txt=cell2mat(txt);
% txt=regexp(txt,'</strong>','split');
% txt=cell2mat(txt);


    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for percent divergent per mouse\n']);

    [output_data] = drgMutiRanksumorTtest(input_sig_data);

%Plot a bar graph showing percent responding to S+
glm_sig=[];
glm_sig_ii=0;
input_sig_data=[];
id_sig_ii=0;

%Find out which files were included (>=10 trials)
files_included_glm=zeros(1,length(handles_out.file));
for fNo=1:length(handles_out.file)
    if handles_out.file(fNo).output_data_odor_Sp.total_trials(1)>=min_tr_resp
        files_included_glm(fNo)=1;
    end
end
files_included_glm=logical(files_included_glm);

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .2 .3])

hold on
xticks([0 1 2])
bar_offset=0;

% rand_offset=0.8;

for grNo=groups
    these_pre=[];
    these_FV=[];
    these_odor=[];

    for mouseNo=unique(handles_out.mouseNo)

        switch handles.group_algo
            case 1
                %Ming
                switch grNo
                    case 1
                        files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                    case 2
                        files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                    case 3
                        files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                end
            case 2
                %Fabio
                files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo)&files_included_glm;
        end
        switch handles_out.use_pFDR
            case {0,1}
                if sum(files_included)>0
                    n_all=zeros(1,sum(files_included));
                    n_all(1,:)=handles_out.output_data_N(files_included);
                    n_odor=zeros(1,sum(files_included));
                    n_odor(1,:)=handles_out.output_data_Sp(files_included);
                    this_odor=n_odor./n_all;
                    these_odor=[these_odor 100*mean(this_odor)];
                end
            case {3,4}
                these_per_resp=[];
                for fileNo=1:length(files_included)
                    if files_included(fileNo)==1
                        this_per_resp=100*sum(handles_out.file(fileNo).output_data_odor_Sp.sig_resp_glm)/length(handles_out.file(1).output_data_odor_Sp.sig_resp_glm);

                        these_per_resp=[these_per_resp this_per_resp];

                    end
                end
                if ~isempty(these_per_resp)
                    these_odor=[these_odor mean(these_per_resp)];
                end
        end

    end

    if ~isempty(these_odor)
        %             bar_offset=bar_offset+1;
        if grNo==1
            bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        else
            bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
        end
        %         bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',7);
        end

        glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
        glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
        glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
        glm_sig_ii=glm_sig_ii+length(these_odor);

        id_sig_ii=id_sig_ii+1;
        input_sig_data(id_sig_ii).data=these_odor;

        switch handles.group_algo
            case 1
                input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];
            case 2
                input_sig_data(id_sig_ii).description=['Odor ' handles.group_names{grNo}];
        end
    end

    bar_offset=bar_offset+1;
end

switch handles.group_algo
    case 1
        %Ming

        xticks([0 1 2])
        labels='xticklabels({';
        for ii_label=[1 3]
            labels=[labels '''' naive_pro{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
        xlim([-0.7 1.7])
        ylim([0 50])
    case 2
        %Fabio
        xticks(2*groups-1)
        labels='xticklabels({';
        for ii_label=1:length(groups)
            labels=[labels '''' handles.group_names{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
end

ylabel('% responsive')
title('Percent ROIs responding to S+ per mouse')

%Perform the glm
fprintf(1, ['\nglm for percent responsive to S+ per mouse\n'])


tbl = table(glm_sig.data',glm_sig.grNo',...
    'VariableNames',{'percent_divergent','percent_correct'});
mdl = fitglm(tbl,'percent_divergent~percent_correct'...
    ,'CategoricalVars',[2])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);



%Plot a bar graph showing percent responding to S-
glm_sig=[];
glm_sig_ii=0;
input_sig_data=[];
id_sig_ii=0;

%Find out which files were included (>=10 trials)
files_included_glm=zeros(1,length(handles_out.file));
for fNo=1:length(handles_out.file)
    if handles_out.file(fNo).output_data_odor_Sp.total_trials(1)>=min_tr_resp
        files_included_glm(fNo)=1;
    end
end
files_included_glm=logical(files_included_glm);

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .2 .3])

hold on
xticks([0 1 2])
bar_offset=0;

% rand_offset=0.8;

for grNo=groups
    these_pre=[];
    these_FV=[];
    these_odor=[];

    for mouseNo=unique(handles_out.mouseNo)
        switch handles.group_algo
            case 1
                %Ming
                switch grNo
                    case 1
                        files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                    case 2
                        files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                    case 3
                        files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
                end
            case 2
                %Fabio
                files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo)&files_included_glm;
        end

        if sum(files_included)>0

            switch handles_out.use_pFDR
                case {0,1}
                    n_all=zeros(1,sum(files_included));
                    n_all(1,:)=handles_out.output_data_N(files_included);
                    n_odor=zeros(1,sum(files_included));
                    n_odor(1,:)=handles_out.output_data_Sm(files_included);
                    this_odor=n_odor./n_all;
                    these_odor=[these_odor 100*mean(this_odor)];

                case {3,4}
                    these_per_resp=[];
                    for fileNo=1:length(files_included)
                        if files_included(fileNo)==1
                            this_per_resp=100*sum(handles_out.file(fileNo).output_data_odor_Sm.sig_resp_glm)/length(handles_out.file(1).output_data_odor_Sm.sig_resp_glm);

                            these_per_resp=[these_per_resp this_per_resp];

                        end
                    end
                    if ~isempty(these_per_resp)
                        these_odor=[these_odor mean(these_per_resp)];
                    end
            end

        end

    end



    %     bar_offset=bar_offset+1;
    %     bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
    if grNo==1
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
    else
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    end
    if length(these_odor)>2
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_odor...
            ,edges,bar_offset,rand_offset,'k','k',7);
    end

    glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
    glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
    glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
    glm_sig_ii=glm_sig_ii+length(these_odor);

    id_sig_ii=id_sig_ii+1;
    input_sig_data(id_sig_ii).data=these_odor;
    switch handles.group_algo
        case 1
            input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];
        case 2
            input_sig_data(id_sig_ii).description=['Odor ' handles.group_names{grNo}];
    end

    bar_offset=bar_offset+1;
end

switch handles.group_algo
    case 1
        %Ming

        xticks([0 1 2])
        labels='xticklabels({';
        for ii_label=[1 3]
            labels=[labels '''' naive_pro{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
        xlim([-0.7 1.7])
        ylim([0 50])
    case 2
        %Fabio
        xticks(2*groups-1)
        labels='xticklabels({';
        for ii_label=1:length(groups)
            labels=[labels '''' handles.group_names{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
end

ylabel('% responsive')
title('Percent ROIs responding to S-')

%Perform the glm
fprintf(1, ['\nglm for percent responsive to S-\n'])


tbl = table(glm_sig.data',glm_sig.grNo',...
    'VariableNames',{'percent_divergent','percent_correct'});
mdl = fitglm(tbl,'percent_divergent~percent_correct'...
    ,'CategoricalVars',[2])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

pffft=1;



%Now generate the summary figures for the divergent responses

if handles_out.all_div_ii_dFF>5
%     
%     %Calculate the crosscorrelations
%     croscorr_traces=corrcoef(handles_out.all_div_dFFspm);
% 
%     Z = linkage(croscorr_traces,'complete','correlation');
% 
%     clusters = cluster(Z,'Maxclust',no_clusters);
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
%     %Do cutoff for 4 clusters
%     cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
%     [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
%     set(H,'LineWidth',2)
%     hFig=figure(figureNo);
%     set(hFig, 'units','normalized','position',[.05 .1 .14 .8])
% 
%     %re-sort the matrix
%     for ii=1:size(croscorr_traces,1)
%         for jj=1:size(croscorr_traces,1)
%             perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
%         end
%     end
% 
% 
% 
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
%     hold on
%     pcolor(perm_croscorr_traces)
%     colormap fire
%     shading flat
% 
%     caxis([-1  1])
%     title(['Cross correlations for all ROIs'])
%     xlim([1 handles_out.all_div_ii_dFF])
%     ylim([1 handles_out.all_div_ii_dFF])
% 
%     %Plot rainbow
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.49 .1 .03 .3])
% 
% 
%     prain=[-1:2/99:1];
%     pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%     %             colormap jet
%     colormap fire
%     shading interp
%     ax=gca;
%     set(ax,'XTickLabel','')
% 
%     %Plot timecourses for all ROIs
%  
% 
%     %S+
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
%     hold on
% 
% 
% 
%     sorted_handles_out.all_div_dFFsplus=[];
%     for ii=1:handles_out.all_div_ii_dFF
%         sorted_handles_out.all_div_dFFsplus(ii,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
%     end
% 
%     time_span_mat=repmat(time_span,ii,1);
%     ROI_mat=repmat(1:ii,length(time_span),1)';
% 
%     pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
%     colormap fire
%     shading flat
% 
%     caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])
% 
%     for ii_te=1:length(timeEvents)
%         plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
%     end
% 
%     xlim([-7 15])
%     ylim([1 ii])
%     title(['S+ for all ROIs'])
%     xlabel('Time (sec)')
%     ylabel('ROI number')
% 
%     %S-
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
%     hold on
% 
% 
%     sorted_handles_out.all_div_dFFsminus=[];
%     for ii=1:handles_out.all_div_ii_dFF
%         sorted_handles_out.all_div_dFFsminus(ii,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
%     end
% 
%     time_span_mat=repmat(time_span,ii,1);
%     ROI_mat=repmat(1:ii,length(time_span),1)';
% 
%     pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
%     colormap fire
%     shading flat
% 
%     caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])
% 
%     for ii_te=1:length(timeEvents)
%         plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
%     end
% 
%     xlim([-7 15])
%     ylim([1 ii])
%     title(['S- for all ROIs'])
%     xlabel('Time (sec)')
%     ylabel('ROI number')
% 
%     %Plot the average timecourses per cluster
%     for clus=1:max(clusters)
%         figureNo = figureNo + 1;
%         try
%             close(figureNo)
%         catch
%         end
%         hFig=figure(figureNo);
% 
%         ax=gca;ax.LineWidth=3;
%         set(hFig, 'units','normalized','position',[.3 .2 .3 .3])
% 
% 
%         hold on
% 
%         %get the dF/F
% 
% 
%         try
% 
%             %S-
%             this_cluster_dFFsminus=[];
%             ii_included=0;
%             for ii=1:handles_out.all_div_ii_dFF
%                 if clusters(ii)==clus
%                     ii_included=ii_included+1;
%                     this_cluster_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(ii,:);
%                 end
%             end
% 
% 
%             CIpv = bootci(1000, @mean, this_cluster_dFFsminus);
%             meanpv=mean(this_cluster_dFFsminus,1);
%             CIpv(1,:)=meanpv-CIpv(1,:);
%             CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%             [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_dFFsminus), CIpv','cmap',[158/255 31/255 99/255]);
% 
%             %S+
%             this_cluster_dFFsplus=[];
%             ii_included=0;
%             for ii=1:handles_out.all_div_ii_dFF
%                 if clusters(ii)==clus
%                     ii_included=ii_included+1;
%                     this_cluster_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(ii,:);
%                 end
%             end
% 
%             CIpv = bootci(1000, @mean, this_cluster_dFFsplus);
%             meanpv=mean(this_cluster_dFFsplus,1);
%             CIpv(1,:)=meanpv-CIpv(1,:);
%             CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%             [hlpvl, hppvl] = boundedline(time_span, mean(this_cluster_dFFsplus), CIpv','cmap',[0 114/255 178/255]);
% 
% 
%             plot(time_span',mean(this_cluster_dFFsminus)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
%             plot(time_span',mean(this_cluster_dFFsplus)','Color',[0 114/255 178/255],'LineWidth',1.5);
% 
% 
% 
%         catch
%         end
% 
%         ylim([-0.5 1.5])
%         xlim([-7 15])
%         this_ylim=ylim;
%         
%         %FV
%          fvhl=plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-'...
%              ,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',5);
%         plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])
% 
%         %Odor on markers
%         plot([0 0],this_ylim,'-k')
%         odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
%         plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')
% 
%         %Reinforcement markers
%         plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
%         reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
%         plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')
% 
% 
%         xlabel('Time(sec)')
%         ylabel('dFF')
% 
% 
%         title(['dFF divergent for cluster no ' num2str(clus)])
% 
%     end


    %Now plot the timecourses per group
    divergence_times_per_group=[];
    dFF_timecourse_per_clus=[]
    for grNo=[1 3]

        switch grNo
            case 1
                ROIs_included=(handles_out.all_div_dFFspm_pcorr>=45)&(handles_out.all_div_dFFspm_pcorr<=65);
            case 2
                ROIs_included=(handles_out.all_div_dFFspm_pcorr>65)&(handles_out.all_div_dFFspm_pcorr<80);
            case 3
                ROIs_included=(handles_out.all_div_dFFspm_pcorr>=80);
        end

        if sum(ROIs_included)>0

            %Calculate the crosscorrelations
            these_all_div_dFFspm=zeros(size(handles_out.all_div_dFFspm,1),sum(ROIs_included));
            these_all_div_dFFspm(:,:)=handles_out.all_div_dFFspm(:,logical(ROIs_included));

            these_all_div_dFFsplus=zeros(size(handles_out.all_div_dFFsplus,2),sum(ROIs_included));
            these_all_div_dFFsplus(:,:)=handles_out.all_div_dFFsplus(logical(ROIs_included),:)';

            these_all_div_dFFsminus=zeros(size(handles_out.all_div_dFFsminus,2),sum(ROIs_included));
            these_all_div_dFFsminus(:,:)=handles_out.all_div_dFFsminus(logical(ROIs_included),:)';

            these_all_delta_dFFsplus=zeros(1,sum(ROIs_included));
            these_all_delta_dFFsplus(1,:)=handles_out.all_div_delta_dFFsplus(1,logical(ROIs_included));

            these_all_delta_dFFsminus=zeros(1,sum(ROIs_included));
            these_all_delta_dFFsminus(1,:)=handles_out.all_div_delta_dFFsminus(1,logical(ROIs_included));

            these_all_div_t=zeros(1,sum(ROIs_included));
            these_all_div_t(1,:)=handles_out.all_div_t(1,logical(ROIs_included));


%             croscorr_traces=corrcoef(these_all_div_dFFspm);
             croscorr_traces=corrcoef(these_all_div_dFFsminus);

            Z = linkage(croscorr_traces,'complete','correlation');

            clusters = cluster(Z,'Maxclust',no_clusters);
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);

            %Do cutoff for no_clusters
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

            caxis([-1  1])
            title(['Cross correlations for all ROIs'])
            xlim([1 sum(ROIs_included)])
            ylim([1 sum(ROIs_included)])

            %S+
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);

            set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
            hold on



            sorted_handles_out.all_div_dFFsplus=[];
            sorted_handles_out.all_div_dFFsplus=these_all_div_dFFsplus(:,outperm)';

            ii_included=size(sorted_handles_out.all_div_dFFsplus,1);

            time_span_mat=repmat(time_span,ii_included,1);
            ROI_mat=repmat(1:ii_included,length(time_span),1)';

            pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
            colormap fire
            shading flat

            caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99.9)])

            for ii_te=1:length(timeEvents)
                plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
            end

            xlim([-7 15])
            ylim([1 ii_included])
            title(['S+ ' per_names{grNo+1}])
            xlabel('Time (sec)')
            ylabel('ROI number')

            %S-
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);

            set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
            hold on



            sorted_handles_out.all_div_dFFsminus=[];
            sorted_handles_out.all_div_dFFsminus=these_all_div_dFFsminus(:,outperm)';

            ii_included=size(sorted_handles_out.all_div_dFFsplus,1);

            time_span_mat=repmat(time_span,ii_included,1);
            ROI_mat=repmat(1:ii_included,length(time_span),1)';

            pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
            colormap fire
            shading flat

            caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

            for ii_te=1:length(timeEvents)
                plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
            end

            xlim([-7 15])
            ylim([1 ii_included])
            title(['S- ' per_names{grNo+1}])
            xlabel('Time (sec)')
            ylabel('ROI number')

            %Plot the average timecourses per cluster
            for clus=1:max(clusters)
                figureNo = figureNo + 1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);

                ax=gca;ax.LineWidth=3;
                set(hFig, 'units','normalized','position',[.3 .2 .3 .3])


                hold on

                %get the dF/F


                try

                    %S-
                    this_cluster_dFFsminus=[];
                    ii_included=0;
                    for ii=1:size(these_all_div_dFFsminus,2)
                        if clusters(ii)==clus
                            ii_included=ii_included+1;
                            this_cluster_dFFsminus(ii_included,:)=these_all_div_dFFsminus(:,ii)';
                        end
                    end
                    
                    dFF_timecourse_per_clus.group(grNo).cluster(clus).sminus_timecourses=this_cluster_dFFsminus;

                    if do_std==0
                        CIpv = bootci(1000, @mean, this_cluster_dFFsminus);
                        meanpv=mean(this_cluster_dFFsminus,1);
                        CIpv(1,:)=meanpv-CIpv(1,:);
                        CIpv(2,:)=CIpv(2,:)-meanpv;
                    else
                        STDpv = std(this_cluster_dFFsminus);
                        meanpv=mean(this_cluster_dFFsminus,1);
                        CIpv(1,:)=STDpv;
                        CIpv(2,:)=STDpv;
                    end


                    [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_dFFsminus), CIpv','cmap',[158/255 31/255 99/255]);

                    %S+
                    this_cluster_dFFsplus=[];
                    ii_included=0;
                    for ii=1:size(these_all_div_dFFsplus,2)
                        if clusters(ii)==clus
                            ii_included=ii_included+1;
                            this_cluster_dFFsplus(ii_included,:)=these_all_div_dFFsplus(:,ii)';
                        end
                    end

                    dFF_timecourse_per_clus.group(grNo).cluster(clus).splus_timecourses=this_cluster_dFFsplus;

                    if do_std==0
                        CIpv = bootci(1000, @mean, this_cluster_dFFsplus);
                        meanpv=mean(this_cluster_dFFsplus,1);
                        CIpv(1,:)=meanpv-CIpv(1,:);
                        CIpv(2,:)=CIpv(2,:)-meanpv;
                    else
                        STDpv = std(this_cluster_dFFsplus);
                        meanpv=mean(this_cluster_dFFsplus,1);
                        CIpv(1,:)=STDpv;
                        CIpv(2,:)=STDpv;
                    end


                    [hlpvl, hppvl] = boundedline(time_span, mean(this_cluster_dFFsplus), CIpv','cmap',[0 114/255 178/255]);


                    plot(time_span',mean(this_cluster_dFFsminus)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
                    plot(time_span',mean(this_cluster_dFFsplus)','Color',[0 114/255 178/255],'LineWidth',1.5);



                catch
                end

                ylim([-0.5 1.5])
                xlim([-7 15])
                this_ylim=ylim;

                %FV
                fvhl=plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-'...
                    ,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',5);
                plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

                %Odor on markers
                plot([0 0],this_ylim,'-k')
                odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
                plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')


                xlabel('Time(sec)')
                ylabel('dFF')


                title(['dFF divergent for ' per_names{grNo+1} ' cluster no ' num2str(clus)])

            end

            %Plot the relationship between delta dFFs with a different
            %color per cluster
%             these_all_delta_dFFsplus=zeros(1,sum(ROIs_included));
%             these_all_delta_dFFsplus(1,:)=handles_out.all_div_delta_dFFsplus(1,logical(ROIs_included));
% 
%             these_all_delta_dFFsminus=zeros(1,sum(ROIs_included));
%             these_all_delta_dFFsminus(1,:)=handles_out.all_div_delta_dFFsminus(1,logical(ROIs_included));
% 
%             these_all_div_t=zeros(1,sum(ROIs_included));
%             these_all_div_t(1,:)=handles_out.all_div_t(1,logical(ROIs_included));

            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);

            set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

            hold on

            %Get the color for each cluster
            linesColor = cell2mat(get(H,'Color')); % get lines color;
            colorList = unique(linesColor, 'rows');

            X_color     = zeros(sum(ROIs_included),3);
            X_cluster   = zeros(sum(ROIs_included),1);
            for iLeaf = 1:sum(ROIs_included)
                [iRow, ~] = find(Z==iLeaf);
                color = linesColor(iRow,:); % !

                % assign color to each observation
                X_color(iLeaf,:) = color;
                % assign cluster number to each observation
                X_cluster(iLeaf,:) = find(ismember(colorList, color, 'rows'));
            end
            cluster_colors=zeros(max(clusters),3);
            for clus= 1:max(clusters)
                this_ii=find(clusters==clus,1,'first');
                cluster_colors(clus,:)=X_color(this_ii,:);
            end

            for clus=1:max(clusters)


                %S-
                this_cluster_delta_dFFsminus=[];
                ii_included=0;
                for ii=1:size(these_all_delta_dFFsminus,2)
                    if clusters(ii)==clus
                        ii_included=ii_included+1;
                        this_cluster_delta_dFFsminus(1,ii_included)=these_all_delta_dFFsminus(1,ii)';
                    end
                end

                %S+
                this_cluster_delta_dFFsplus=[];
                ii_included=0;
                for ii=1:size(these_all_delta_dFFsplus,2)
                    if clusters(ii)==clus
                        ii_included=ii_included+1;
                        this_cluster_delta_dFFsplus(1,ii_included)=these_all_delta_dFFsplus(1,ii)';
                    end
                end
                this_color=zeros(1,3);
                this_color(1,:)=cluster_colors(clus,:);
                plot(this_cluster_delta_dFFsminus,this_cluster_delta_dFFsplus,'.','Color',this_color)
            end


            ylabel('Delta zdFF S+')
            xlabel('Delta zdFF S-')
            title(['Delta zdFF S+ vs S- for' per_names{grNo+1}])




            %Divergence time
            figureNo=figureNo+1;
            try
                close(figureNo)
            catch
            end

            hFig = figure(figureNo);

            set(hFig, 'units','normalized','position',[.2 .2 .4 .3])

            edges=[-1.5:0.5:6.5];

            histogram(these_all_div_t,edges,'Normalization','probability')
            divergence_times_per_group.group(grNo).div_t=these_all_div_t;

            xlim([-3 7])
            %             ylim([0 120])
            ylabel('No of ROIs')
            xlabel('Divergence time')
            title(['Histogram for divergence time for ' per_names{grNo+1}])
        end
    end

    %Quantify changes in dFF
    pffft=1;

    cluster_pairing.group(1).clus(1)=1;
    cluster_pairing.group(1).clus(2)=2;
    cluster_pairing.group(1).clus(3)=2;

    cluster_pairing.group(3).clus(1)=2;
    cluster_pairing.group(3).clus(2)=1;
    cluster_pairing.group(3).clus(3)=2;

    clus_per_clus_ii=[1 1 2 3];
    time_window_clus=[1.5 10 1.3 0];

    clus_ii_legends{1}='cluster 1 peak 1';
    clus_ii_legends{2}='cluster 1 peak 2';
    clus_ii_legends{3}='cluster 2';
    clus_ii_legends{4}='cluster 3';

    for clus_ii=1:4

        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        hold on

        clus=clus_per_clus_ii(clus_ii);
        
        this_time_window=[time_window_clus(clus_ii)-0.5 time_window_clus(clus_ii)+0.5];

        glm_clus=[];
        glm_clus_ii=0;
        input_clus_data=[];
        id_clus_ii=0;
        bar_offset=0;

        edges=[-1:0.05:5];
        rand_offset=0.4;

        for grNo =[1 3]
 
            actual_clus=cluster_pairing.group(grNo).clus(clus);

            %do S-
            these_dFF_minus=mean(dFF_timecourse_per_clus.group(grNo).cluster(actual_clus).sminus_timecourses...
                (:,(time_span>=this_time_window(1))&(time_span<=this_time_window(2))),2);

            bar_offset=bar_offset+1;
            bar(bar_offset,mean(these_dFF_minus),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_dFF_minus...
                ,edges,bar_offset,rand_offset,'k','k',1);

            glm_clus.data(glm_clus_ii+1:glm_clus_ii+length(these_dFF_minus))=these_dFF_minus;
            glm_clus.pcorr(glm_clus_ii+1:glm_clus_ii+length(these_dFF_minus))=grNo*ones(1,length(these_dFF_minus));
            glm_clus.spm(glm_clus_ii+1:glm_clus_ii+length(these_dFF_minus))=0*ones(1,length(these_dFF_minus));
            glm_clus_ii=glm_clus_ii+length(these_dFF_minus);

            id_clus_ii=id_clus_ii+1;
            input_clus_data(id_clus_ii).data=these_dFF_minus;
            input_clus_data(id_clus_ii).description=[naive_pro{grNo} ' S-'];

            %do S+
            these_dFF_plus=mean(dFF_timecourse_per_clus.group(grNo).cluster(actual_clus).splus_timecourses...
                (:,(time_span>=this_time_window(1))&(time_span<=this_time_window(2))),2);

            bar_offset=bar_offset+1;
            bar(bar_offset,mean(these_dFF_plus),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_dFF_plus...
                ,edges,bar_offset,rand_offset,'k','k',1);

            bar_offset=bar_offset+1;

            glm_clus.data(glm_clus_ii+1:glm_clus_ii+length(these_dFF_plus))=these_dFF_plus;
            glm_clus.pcorr(glm_clus_ii+1:glm_clus_ii+length(these_dFF_plus))=grNo*ones(1,length(these_dFF_plus));
            glm_clus.spm(glm_clus_ii+1:glm_clus_ii+length(these_dFF_plus))=1*ones(1,length(these_dFF_plus));
            glm_clus_ii=glm_clus_ii+length(these_dFF_plus);

            id_clus_ii=id_clus_ii+1;
            input_clus_data(id_clus_ii).data=these_dFF_plus;
            input_clus_data(id_clus_ii).description=[naive_pro{grNo} ' S-'];
        end

        xticks([1 2 4 5])
        xticklabels({'S-','S+','S-','S+'})
      
        ylabel('zdFF')

        title(['zdFF for ' clus_ii_legends{clus_ii}])



        %Perform the glm for this cluster
        fprintf(1, ['\nglm for zdFF for '  clus_ii_legends{clus_ii} '\n'])
        fprintf(fileID, ['\nglm for zdFF for '  clus_ii_legends{clus_ii} '\n']);

        tbl = table(glm_clus.data',glm_clus.pcorr',glm_clus.spm',...
            'VariableNames',{'zdFF','p_correct','spm'});
        mdl = fitglm(tbl,'zdFF~p_correct+spm+p_correct*spm'...
            ,'CategoricalVars',[2,3])

        txt = evalc('mdl');
        txt=regexp(txt,'<strong>','split');
        txt=cell2mat(txt);
        txt=regexp(txt,'</strong>','split');
        txt=cell2mat(txt);

        fprintf(fileID,'%s\n', txt);



        %Do the ranksum/t-test
        fprintf(1, ['\n\nRanksum or t-test p values for zdFF for ' clus_ii_legends{clus_ii} '\n'])
        fprintf(fileID, ['\n\nRanksum or t-test p values for zdFF for ' clus_ii_legends{clus_ii} '\n']);

        [output_data] = drgMutiRanksumorTtest(input_clus_data, fileID,0);

    end

    %Divergence time
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    hold on

    [f_aic,x_aic] = drg_ecdf(divergence_times_per_group.group(1).div_t(divergence_times_per_group.group(1).div_t<=7)');
    plot(x_aic,f_aic,'Color',[80/255 194/255 255/255],'LineWidth',2)
    plot([-3 -1.5],[0 0],'Color',[80/255 194/255 255/255],'LineWidth',2)

    try
        [f_aic,x_aic] = drg_ecdf(divergence_times_per_group.group(3).div_t(divergence_times_per_group.group(3).div_t<=7)');
        plot(x_aic,f_aic,'Color',[0 114/255 178/255],'LineWidth',2)
        plot([-3 -1.5],[0 0],'Color',[0 114/255 178/255],'LineWidth',2)
    catch
    end

    xlim([-3 7])
    ylim([-0.1 1.2])

    this_ylim=ylim;

    this_ylim=ylim;

    %FV
    fvhl=plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-'...
        ,'Color',[0.9 0.9 0.9],'LineWidth',5);
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

    ylabel('Probability')
    xlabel('Divegence time (sec)')
    title(['Divergence time cumulative probability'])


     %Divergence time histogram
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
    hold on

    edges=[-1.5:0.5:6.5];


    histogram(divergence_times_per_group.group(1).div_t(divergence_times_per_group.group(1).div_t<=7)',edges,...
        'Normalization','probability','FaceColor',[80/255 194/255 255/255])

    histogram(divergence_times_per_group.group(3).div_t(divergence_times_per_group.group(1).div_t<=7)',edges,...
        'Normalization','probability','FaceColor',[0 114/255 178/255])

    plot([-3 7],[0 0],'-k')
    xlim([-2 7])
    ylim([-0.04 0.18])

    this_ylim=ylim;

    %FV
    fvhl=plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-'...
        ,'Color',[0.9 0.9 0.9],'LineWidth',5);
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

    ylabel('No of ROIs')
    xlabel('Divergence time')
    title(['Histogram for divergence time'])

    try
        input_data=[];

        input_data(1).data=divergence_times_per_group.group(1).div_t(divergence_times_per_group.group(1).div_t<=7)';
        input_data(1).description='Naive';
        input_data(2).data=divergence_times_per_group.group(3).div_t(divergence_times_per_group.group(1).div_t<=7)';
        input_data(2).description='Proficient';

        fprintf(1, ['\n\nRanksum or t-test p values for discrimination times\n'])

        [output_data] = drgMutiRanksumorTtest(input_data);
    catch
    end

    %Plot rainbow
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


    prain=[prctile(handles_out.all_div_dFFspm(:),1):(prctile(handles_out.all_div_dFFspm(:),99)-prctile(handles_out.all_div_dFFspm(:),1))/99:prctile(handles_out.all_div_dFFspm(:),99)];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Generate a histogram for the magnitude of the responses

    %     handles_out.all_div_t=[];
    %     handles_out.all_div_delta_dFFsplus=[];
    %     handles_out.all_div_delta_dFFsminus=[];

%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
% 
%     edges=-0.5:0.5:4;
%     histogram(handles_out.all_div_delta_dFFsplus,edges,'FaceColor',[0 114/255 178/255])
%     hold on
%     histogram(handles_out.all_div_delta_dFFsminus,edges,'FaceColor',[158/255 31/255 99/255])
% 
%     ylabel('No of ROIs')
%     xlabel('Delta zdFF')
%     title('Histogram for Delta zdFF')
% 
%     this_ylim=ylim;
% 
%     text(2.5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
%     text(2.5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)
% 
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%     plot(handles_out.all_div_delta_dFFsminus,handles_out.all_div_delta_dFFsplus,'ok')
% 
%     ylabel('S+')
%     xlabel('S-')
%     title('Delta zdFF')
% 
%     %Divergence time
%     figureNo=figureNo+1;
%     try
%         close(figureNo)
%     catch
%     end
% 
%     hFig = figure(figureNo);
% 
%     set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
% 
%     edges=[-1.5:0.5:6.5];
%     histogram(handles_out.all_div_t(~isnan(handles_out.all_div_t)),edges)
% 
%     xlim([-3 7])
%     ylim([0 120])
%     ylabel('No of ROIs')
%     xlabel('Divegence time')
%     title('Histogram for divergence time')

end


%Now generate the summary figures for the divergent responses

if handles_out.all_spresp_ii_dFF+handles_out.all_smresp_ii_dFF>5
    %Calculate the crosscorrelations
    %     croscorr_traces=abs(corrcoef(handles_out.all_div_dFFspm)); %please note that I am using the absolute value
    all_spandmresp_dFFspm=[handles_out.all_spresp_dFFspm handles_out.all_smresp_dFFspm];
    all_spandmresp_iidFF=[1:handles_out.all_spresp_ii_dFF 1:handles_out.all_smresp_ii_dFF];
    all_spandmresp_sporsm=[ones(1,handles_out.all_spresp_ii_dFF) zeros(1,handles_out.all_smresp_ii_dFF)];
    croscorr_traces=corrcoef(all_spandmresp_dFFspm);

    %Set autocorrelations to zero
%     for ii=1:size(croscorr_traces,1)
%         croscorr_traces(ii,ii)=0;
%     end
    Z = linkage(croscorr_traces,'complete','correlation');
    no_clusters=2; %Note that this is different
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
    xlim([1 handles_out.all_spresp_ii_dFF+handles_out.all_smresp_ii_dFF])
    ylim([1 handles_out.all_spresp_ii_dFF+handles_out.all_smresp_ii_dFF])

    %Plot rainbow
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


    prain=[0:0.6/99:0.6];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Plot timecourses for all ROIs for S+ and S-

    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

     set(hFig, 'units','normalized','position',[.05 .1 .18 .8])
    hold on



    sorted_handles_out.all_spmresp_ii_dFF=0;
    sorted_handles_out.all_spmresp_sporsm=[];
    sorted_handles_out.all_sporsm_dFFspm=[];
    sorted_handles_out.clusters=[];
    for ii=1:length(all_spandmresp_iidFF)
        this_ii_dFF=all_spandmresp_iidFF(outperm(ii));
        this_sporsm=all_spandmresp_sporsm(outperm(ii));
        sorted_handles_out.all_spmresp_ii_dFF=sorted_handles_out.all_spmresp_ii_dFF+1;
        if this_sporsm==0
            %This is S-
            sorted_handles_out.all_sporsm_dFFspm(sorted_handles_out.all_spmresp_ii_dFF,:)=handles_out.all_smresp_dFFsminus(this_ii_dFF,:);
        else
            %This is S+
            sorted_handles_out.all_sporsm_dFFspm(sorted_handles_out.all_spmresp_ii_dFF,:)=handles_out.all_spresp_dFFsplus(this_ii_dFF,:);
        end
        sorted_handles_out.clusters(sorted_handles_out.all_spmresp_ii_dFF)=clusters(outperm(ii));
    end

    %pcolor does not show the first row
    pseudo_dFFspm=zeros(size(sorted_handles_out.all_sporsm_dFFspm,1)+1,size(sorted_handles_out.all_sporsm_dFFspm,2));
    pseudo_dFFspm(1:end-1,:)=sorted_handles_out.all_sporsm_dFFspm;

    time_span_mat=repmat(time_span,sorted_handles_out.all_spmresp_ii_dFF+1,1);
    ROI_mat=repmat(1:sorted_handles_out.all_spmresp_ii_dFF+1,length(time_span),1)';

    pcolor(time_span_mat,ROI_mat,pseudo_dFFspm)


    %         time_span_mat=repmat(time_span,sorted_handles_out.all_spmresp_ii_dFF,1);
    %         ROI_mat=repmat(1:sorted_handles_out.all_spmresp_ii_dFF,length(time_span),1)';
    %
    %         pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_sporsm_dFFspm)

    colormap fire
    shading flat

    caxis([prctile(sorted_handles_out.all_sporsm_dFFspm(:),1) prctile(sorted_handles_out.all_sporsm_dFFspm(:),99)])

    for ii_te=1:length(timeEvents)
        plot([timeEvents(ii_te) timeEvents(ii_te)],[0 sorted_handles_out.all_spmresp_ii_dFF],'-r')
    end

    xlim([-7 15])
    ylim([1 ii])
    title(['Timecourses for significant dFF changes, all ROIs'])
    xlabel('Time (sec)')
    ylabel('ROI number')



    %Plot the average timecourses per cluster
    for clus=1:max(clusters)
        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.3 .2 .3 .3])


        hold on

        %get the dF/F


        try


            this_cluster_dFF=[];
            ii_included=0;
            for ii=1:length(sorted_handles_out.clusters)
                if sorted_handles_out.clusters(ii)==clus
                    ii_included=ii_included+1;
                    this_cluster_dFF(ii_included,:)=sorted_handles_out.all_sporsm_dFFspm(ii,:);
                end
            end


            CIpv = bootci(1000, @mean, this_cluster_dFF);
            meanpv=mean(this_cluster_dFF,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;


            [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_dFF), CIpv','cmap',[0 0 0]);



        catch
        end

        this_ylim=ylim;
        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
        end

        xlim([-7 15])


        xlabel('Time(sec)')
        ylabel('dFF')


        title(['dFF for cluster no ' num2str(clus)])

    end


    %         %Now plot the timecourses per group
    %         for grNo=1:3
    %
    %             switch grNo
    %                 case 1
    %                     ROIs_included=(handles_out.all_div_dFFspm_pcorr>=45)&(handles_out.all_div_dFFspm_pcorr<=65);
    %                 case 2
    %                     ROIs_included=(handles_out.all_div_dFFspm_pcorr>65)&(handles_out.all_div_dFFspm_pcorr<80);
    %                 case 3
    %                     ROIs_included=(handles_out.all_div_dFFspm_pcorr>=80);
    %             end
    %
    %             if sum(ROIs_included)>0
    %
    %                 %S+
    %                 figureNo=figureNo+1;
    %                 try
    %                     close(figureNo)
    %                 catch
    %                 end
    %
    %                 hFig = figure(figureNo);
    %
    %                 set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
    %                 hold on
    %
    %
    %
    %                 sorted_handles_out.all_div_dFFsplus=[];
    %                 ii_included=0;
    %                 for ii=1:handles_out.all_div_ii_dFF
    %                     if ROIs_included(ii)
    %                         ii_included=ii_included+1;
    %                         sorted_handles_out.all_div_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
    %                     end
    %                 end
    %
    %                 time_span_mat=repmat(time_span,ii_included,1);
    %                 ROI_mat=repmat(1:ii_included,length(time_span),1)';
    %
    %                 pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
    %                 colormap fire
    %                 shading flat
    %
    %                 caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])
    %
    %                 for ii_te=1:length(timeEvents)
    %                     plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
    %                 end
    %
    %                 xlim([-7 15])
    %                 ylim([1 ii_included])
    %                 title(['S+ ' per_names{grNo+1}])
    %                 xlabel('Time (sec)')
    %                 ylabel('ROI number')
    %
    %                 %S-
    %                 figureNo=figureNo+1;
    %                 try
    %                     close(figureNo)
    %                 catch
    %                 end
    %
    %                 hFig = figure(figureNo);
    %
    %                 set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
    %                 hold on
    %
    %
    %
    %                 sorted_handles_out.all_div_dFFsminus=[];
    %                 ii_included=0;
    %                 for ii=1:handles_out.all_div_ii_dFF
    %                     if ROIs_included(ii)
    %                         ii_included=ii_included+1;
    %                         sorted_handles_out.all_div_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
    %                     end
    %                 end
    %
    %                 time_span_mat=repmat(time_span,ii_included,1);
    %                 ROI_mat=repmat(1:ii_included,length(time_span),1)';
    %
    %                 pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
    %                 colormap fire
    %                 shading flat
    %
    %                 caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])
    %
    %                 for ii_te=1:length(timeEvents)
    %                     plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
    %                 end
    %
    %                 xlim([-7 15])
    %                 ylim([1 ii_included])
    %                 title(['S- ' per_names{grNo+1}])
    %                 xlabel('Time (sec)')
    %                 ylabel('ROI number')
    %             end
    %         end
    %
    %         %Plot rainbow
    %         figureNo=figureNo+1;
    %         try
    %             close(figureNo)
    %         catch
    %         end
    %
    %         hFig = figure(figureNo);
    %
    %         set(hFig, 'units','normalized','position',[.49 .1 .05 .3])
    %
    %
    %         prain=[prctile(handles_out.all_div_dFFspm(:),1):(prctile(handles_out.all_div_dFFspm(:),99)-prctile(handles_out.all_div_dFFspm(:),1))/99:prctile(handles_out.all_div_dFFspm(:),99)];
    %         pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %         %             colormap jet
    %         colormap fire
    %         shading interp
    %         ax=gca;
    %         set(ax,'XTickLabel','')
    %
    %         %Generate a histogram for the magnitude of the responses
    %
    %         %     handles_out.all_div_t=[];
    %         %     handles_out.all_div_delta_dFFsplus=[];
    %         %     handles_out.all_div_delta_dFFsminus=[];
    %
    %         figureNo=figureNo+1;
    %         try
    %             close(figureNo)
    %         catch
    %         end
    %
    %         hFig = figure(figureNo);
    %
    %         set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
    %
    %         edges=-0.5:0.5:4;
    %         histogram(handles_out.all_div_delta_dFFsplus,edges,'FaceColor',[0 114/255 178/255])
    %         hold on
    %         histogram(handles_out.all_div_delta_dFFsminus,edges,'FaceColor',[158/255 31/255 99/255])
    %
    %         ylabel('No of ROIs')
    %         xlabel('Delta zdFF')
    %         title('Histogram for Delta zdFF')
    %
    %         this_ylim=ylim;
    %
    %         text(2.5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
    %         text(2.5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)
    %
    %         figureNo=figureNo+1;
    %         try
    %             close(figureNo)
    %         catch
    %         end
    %
    %         hFig = figure(figureNo);
    %
    %         set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
    %
    %         plot(handles_out.all_div_delta_dFFsminus,handles_out.all_div_delta_dFFsplus,'ok')
    %
    %         ylabel('S+')
    %         xlabel('S-')
    %         title('Delta zdFF')
    %
    %         %Divergence time
    %         figureNo=figureNo+1;
    %         try
    %             close(figureNo)
    %         catch
    %         end
    %
    %         hFig = figure(figureNo);
    %
    %         set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
    %
    %         edges=[-1.75:0.5:5.25];
    %         histogram(handles_out.all_div_t(~isnan(handles_out.all_div_t)),edges)
    %
    %         ylabel('No of ROIs')
    %         xlabel('Divegence time')
    %         title('Histogram for divergence time')

end
%Uncomment this if you want to browse through the figures
%     tic
%     for figNo=1:figureNo
%         figure(figNo)
%         this_toc=toc;
%         while toc-this_toc<4
%         end
%     end

fclose(fileID)
pffft=1;
