
clear all
close all

 
[outFileName,outhPathName] = uigetfile({'drgCaImAn_pval_choices*.mat'},'Select the .mat file with drgCaImAn_pval_batch');
load([outhPathName outFileName])

time_span=handles_out.time_span;

s=rng;

no_clusters=6;

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

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

pre_t=[-3 -2];
fv_t=[-1 0];
odor_t=[2 4];

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [-1.5 0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

%Time tinterval for shifting time base due to slow olfactometer computer
t_shift=0;

%Plot a bar graph showing percent signifcant divergence
glm_sig=[];
glm_sig_ii=0;
input_sig_data=[];
id_sig_ii=0;

figureNo=0;
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

hold on
bar_offset=0;
edges=[0:10:20];
rand_offset=0.5;

switch handles.group_algo
    case 1
        %Ming
        groups=[1:3];
    case 2
        %Fabio
        groups=unique(handles.group);
end

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
                        files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                    case 2
                        files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                    case 3
                        files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
                end
            case 2
                %Fabio
                files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo);
        end

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

    end

    if ~isempty(these_odor)
        %             bar_offset=bar_offset+1;
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',5);
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
        for ii_label=2:4
            labels=[labels '''' per_names{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
        xlim([-0.7 2.7])
        ylim([0 20])
        
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

%Perform the glm
fprintf(1, ['\nglm for percent divergent per mouse\n'])


tbl = table(glm_sig.data',glm_sig.grNo',...
    'VariableNames',{'percent_divergent','percent_correct'});
mdl = fitglm(tbl,'percent_divergent~percent_correct'...
    ,'CategoricalVars',[2])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);


%Plot a bar graph showing percent responding to S+
glm_sig=[];
glm_sig_ii=0;
input_sig_data=[];
id_sig_ii=0;

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

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
                        files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                    case 2
                        files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                    case 3
                        files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
                end
            case 2
                %Fabio
                files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo);
        end

        if sum(files_included)>0
            n_all=zeros(1,sum(files_included));
            n_all(1,:)=handles_out.output_data_N(files_included);
            n_odor=zeros(1,sum(files_included));
            n_odor(1,:)=handles_out.output_data_Sp(files_included);
            this_odor=n_odor./n_all;
            these_odor=[these_odor 100*mean(this_odor)];
        end

    end

    if ~isempty(these_odor)
        %             bar_offset=bar_offset+1;
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',5);
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
        for ii_label=2:4
            labels=[labels '''' per_names{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
        xlim([-0.7 2.7])
        ylim([0 20])
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

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

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
                        files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                    case 2
                        files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                    case 3
                        files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
                end
            case 2
                %Fabio
                files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo);
        end

        if sum(files_included)>0
            n_all=zeros(1,sum(files_included));
            n_all(1,:)=handles_out.output_data_N(files_included);
            n_odor=zeros(1,sum(files_included));
            n_odor(1,:)=handles_out.output_data_Sm(files_included);
            this_odor=n_odor./n_all;
            these_odor=[these_odor 100*mean(this_odor)];
        end

    end



%     bar_offset=bar_offset+1;
    bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
    if length(these_odor)>2
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_odor...
            ,edges,bar_offset,rand_offset,'k','k',5);
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
        for ii_label=2:4
            labels=[labels '''' per_names{ii_label} ''', '];
        end
        labels=[labels(1:end-2) '})'];
        eval(labels)
        xlim([-0.7 2.7])
        ylim([0 20])
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
    %Calculate the crosscorrelations
    %     croscorr_traces=abs(corrcoef(handles_out.all_div_dFFspm)); %please note that I am using the absolute value
    croscorr_traces=corrcoef(handles_out.all_div_dFFspm);

    %Set autocorrelations to zero
%     for ii=1:size(croscorr_traces,1)
%         croscorr_traces(ii,ii)=0;
%     end
    Z = linkage(croscorr_traces,'complete','correlation');
    
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

    caxis([-1  1])
    title(['Cross correlations for all ROIs'])
    xlim([1 handles_out.all_div_ii_dFF])
    ylim([1 handles_out.all_div_ii_dFF])

    %Plot rainbow
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.49 .1 .03 .3])


    prain=[0:0.6/99:0.6];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Plot timecourses for all ROIs
 

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
    for ii=1:handles_out.all_div_ii_dFF
        sorted_handles_out.all_div_dFFsplus(ii,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
    end

    time_span_mat=repmat(time_span,ii,1);
    ROI_mat=repmat(1:ii,length(time_span),1)';

    pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
    colormap fire
    shading flat

    caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

    for ii_te=1:length(timeEvents)
        plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
    end

    xlim([-7 15])
    ylim([1 ii])
    title(['S+ for all ROIs'])
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
    for ii=1:handles_out.all_div_ii_dFF
        sorted_handles_out.all_div_dFFsminus(ii,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
    end

    time_span_mat=repmat(time_span,ii,1);
    ROI_mat=repmat(1:ii,length(time_span),1)';

    pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
    colormap fire
    shading flat

    caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

    for ii_te=1:length(timeEvents)
        plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
    end

    xlim([-7 15])
    ylim([1 ii])
    title(['S- for all ROIs'])
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
            for ii=1:handles_out.all_div_ii_dFF
                if clusters(ii)==clus
                    ii_included=ii_included+1;
                    this_cluster_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(ii,:);
                end
            end


            CIpv = bootci(1000, @mean, this_cluster_dFFsminus);
            meanpv=mean(this_cluster_dFFsminus,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;


            [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_dFFsminus), CIpv','cmap',[158/255 31/255 99/255]);

            %S+
            this_cluster_dFFsplus=[];
            ii_included=0;
            for ii=1:handles_out.all_div_ii_dFF
                if clusters(ii)==clus
                    ii_included=ii_included+1;
                    this_cluster_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(ii,:);
                end
            end

            CIpv = bootci(1000, @mean, this_cluster_dFFsplus);
            meanpv=mean(this_cluster_dFFsplus,1);
            CIpv(1,:)=meanpv-CIpv(1,:);
            CIpv(2,:)=CIpv(2,:)-meanpv;


            [hlpvl, hppvl] = boundedline(time_span, mean(this_cluster_dFFsplus), CIpv','cmap',[0 114/255 178/255]);


            plot(time_span',mean(this_cluster_dFFsminus)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
            plot(time_span',mean(this_cluster_dFFsplus)','Color',[0 114/255 178/255],'LineWidth',1.5);



        catch
        end

        ylim([-0.5 2.5])
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


        title(['dFF divergent for cluster no ' num2str(clus)])

    end


    %Now plot the timecourses per group
    for grNo=1:3

        switch grNo
            case 1
                ROIs_included=(handles_out.all_div_dFFspm_pcorr>=45)&(handles_out.all_div_dFFspm_pcorr<=65);
            case 2
                ROIs_included=(handles_out.all_div_dFFspm_pcorr>65)&(handles_out.all_div_dFFspm_pcorr<80);
            case 3
                ROIs_included=(handles_out.all_div_dFFspm_pcorr>=80);
        end

        if sum(ROIs_included)>0
            try
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
                ii_included=0;
                for ii=1:handles_out.all_div_ii_dFF
                    if ROIs_included(ii)
                        ii_included=ii_included+1;
                        sorted_handles_out.all_div_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
                    end
                end

                time_span_mat=repmat(time_span,ii_included,1);
                ROI_mat=repmat(1:ii_included,length(time_span),1)';

                pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
                colormap fire
                shading flat

                caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

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
                ii_included=0;
                for ii=1:handles_out.all_div_ii_dFF
                    if ROIs_included(ii)
                        ii_included=ii_included+1;
                        sorted_handles_out.all_div_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
                    end
                end

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

            catch
            end
        end
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

    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.2 .2 .4 .3])

    edges=-0.5:0.5:4;
    histogram(handles_out.all_div_delta_dFFsplus,edges,'FaceColor',[0 114/255 178/255])
    hold on
    histogram(handles_out.all_div_delta_dFFsminus,edges,'FaceColor',[158/255 31/255 99/255])

    ylabel('No of ROIs')
    xlabel('Delta zdFF')
    title('Histogram for Delta zdFF')

    this_ylim=ylim;

    text(2.5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
    text(2.5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)

    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    plot(handles_out.all_div_delta_dFFsminus,handles_out.all_div_delta_dFFsplus,'ok')

    ylabel('S+')
    xlabel('S-')
    title('Delta zdFF')

    %Divergence time
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.2 .2 .4 .3])

    edges=[-1.5:0.5:6.5];
    histogram(handles_out.all_div_t(~isnan(handles_out.all_div_t)),edges)

    xlim([-3 7])
    ylim([0 120])
    ylabel('No of ROIs')
    xlabel('Divegence time')
    title('Histogram for divergence time')

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
pffft=1;
