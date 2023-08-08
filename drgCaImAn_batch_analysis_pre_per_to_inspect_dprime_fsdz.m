function drgCaImAn_batch_analysis_pre_per_to_inspect_dprime_fsdz
%Note: fitcnet will not work in Matlab versions earlier than 2021a

close all
clear all

first_file=1; %Do not change

prof_labels{1}='naive';
prof_labels{2}='intermediate';
prof_labels{3}='proficient';

% time_periods_eu=[-5 -3;
%     -2.5 -1.5;
%     -1 0;
%     2 4.1];

time_periods_eu=[
    -4.5 -3.5
    -1 0;
    3.1 4.1;
    4.5 5.5];

period_labels{1}='Baseline';
period_labels{2}='Pre';
period_labels{3}='Odor';
period_labels{4}='Reinf';



figNo=0;

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_inspect_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_inspect_fsdz run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

fileID = fopen([choiceBatchPathName 'drgCaImAn_batch_analysis_pre_per_to_inspect_fsdz.txt'],'w');

%Parallel batch processing for each file
all_files_present=1;
for filNum=1:handles.no_files


    %Make sure that all the files exist
    pre_per_FileName=handles.FileName_pre_per{filNum};
    if iscell(handles.PathName_pre_per)
        pre_per_PathName=handles.PathName_pre_per{filNum};
    else
        pre_per_PathName=handles.PathName_pre_per;
    end

    if exist([pre_per_PathName pre_per_FileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' pre_per_PathName pre_per_FileName ' does not exist\n'],filNum);
        all_files_present=0;
    end

end


no_mice=1;
these_mice=[];
these_mice{1}=handles.mouse{1};

for fileNo=2:length(handles.FileName_pre_per)
    mouse_found=0;
    for ii_mouse=1:length(these_mice)
        if strcmp(these_mice{ii_mouse},handles.mouse{fileNo})
            mouse_found=1;
        end
    end
    if mouse_found==0
        no_mice=no_mice+1;
        these_mice{no_mice}=handles.mouse{fileNo};
    end
end


if all_files_present==1


    %Process each file separately
    all_data=[];
    for ii_pc_group=1:3
        for ii_mouse=1:length(these_mice)
            all_data.mouse(ii_mouse).p_corr(ii_pc_group).ii_movies=0;
        end
    end

    first_time=1;
    for fileNo=1:length(handles.FileName_pre_per)
        tic
        first_toc=toc;

        for ii_mouse=1:length(these_mice)
            if strcmp(handles.mouse{fileNo},these_mice{ii_mouse})
                this_ii_mouse=ii_mouse;
            end
        end

        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};

        old_handles=handles;
        if fileNo==1
            this_pre_per_name=[pre_per_PathName pre_per_FileName];
            load(this_pre_per_name)
        end
        handles=old_handles;

        for ii_pc_group=1:3
            this_full_name=[pre_per_PathName pre_per_FileName(1:end-11) 'insp' num2str(ii_pc_group) '.mat'];
            if exist(this_full_name)~=0
                handles_out=[];
                load(this_full_name)
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).ii_movies=all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).ii_movies+1;
                ii_movies=all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).ii_movies;
                if first_time==1
                    time_span=handles_out.time_span;
                    trimmed_time_span=time_span((time_span>=-7)&(time_span<=15));
                    all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).trimmed_time_span=trimmed_time_span;
                    first_time=0;
                end
                %I do this because images are not always taken with the
                %same fps
                this_time_span=handles_out.time_span;

                if isequal(this_time_span,time_span)
                    zdist_euclid=handles_out.zdist_euclid((time_span>=-7)&(time_span<=15));
                    zdist_euclid_within=handles_out.zdist_euclid_within((time_span>=-7)&(time_span<=15));
                    meandFF_sp=mean(handles_out.meandFF_per_trial_sp(:,(time_span>=-7)&(time_span<=15)),1);
                    meandFF_sm=mean(handles_out.meandFF_per_trial_sm(:,(time_span>=-7)&(time_span<=15)),1);
                    KLdivergence=handles_out.KLdivergence((time_span>=-7)&(time_span<=15));
                    KLdivergence_mix=handles_out.KLdivergence_mix((time_span>=-7)&(time_span<=15));
                    dprime=handles_out.dprime((time_span>=-7)&(time_span<=15));
                    dprime_mix=handles_out.dprime_mix((time_span>=-7)&(time_span<=15));
                else
                    %time_span has a different time base
                    this_meandFF_sp=mean(handles_out.meandFF_per_trial_sp,1);
                    this_meandFF_sm=mean(handles_out.meandFF_per_trial_sm,1);
                    meandFF_sp=zeros(1,length(trimmed_time_span));
                    meandFF_sm=zeros(1,length(trimmed_time_span));
                    zdist_euclid=zeros(1,length(trimmed_time_span));
                    zdist_euclid_within=zeros(1,length(trimmed_time_span));
                    KLdivergence=zeros(1,length(trimmed_time_span));
                    KLdivergence_mix=zeros(1,length(trimmed_time_span));
                    dprime=zeros(1,length(trimmed_time_span));
                    dprime_mix=zeros(1,length(trimmed_time_span));

                    for ii_t=1:length(trimmed_time_span)
                        if ii_t<length(trimmed_time_span)
                            these_jj=find((this_time_span>=trimmed_time_span(ii_t))&(this_time_span<trimmed_time_span(ii_t+1)));
                            if ~isempty(these_jj)
                                zdist_euclid(ii_t)=mean(handles_out.zdist_euclid(these_jj));
                                zdist_euclid_within(ii_t)=mean(handles_out.zdist_euclid_within(these_jj));
                                KLdivergence(ii_t)=mean(handles_out.KLdivergence(these_jj));
                                KLdivergence_mix(ii_t)=mean(handles_out.KLdivergence_mix(these_jj));
                                dprime(ii_t)=mean(handles_out.dprime(these_jj));
                                dprime_mix(ii_t)=mean(handles_out.dprime_mix(these_jj));
                                meandFF_sp(ii_t)=mean(this_meandFF_sp(these_jj));
                                meandFF_sm(ii_t)=mean(this_meandFF_sm(these_jj));
                            else
                                this_jj=find((this_time_span<=trimmed_time_span(ii_t)),1,"first");
                                zdist_euclid(ii_t)=zdist_euclid(this_jj);
                                zdist_euclid_within(ii_t)=zdist_euclid_within(this_jj);
                                KLdivergence(ii_t)=KLdivergence(this_jj);
                                KLdivergence_mix(ii_t)=KLdivergence_mix(this_jj);
                                dprime(ii_t)=dprime(this_jj);
                                dprime_mix(ii_t)=dprime_mix(this_jj);
                                meandFF_sp(ii_t)=this_meandFF_sp(this_jj);
                                meandFF_sm(ii_t)=this_meandFF_sm(this_jj);
                            end
                        else
                            this_jj=find((this_time_span>=trimmed_time_span(ii_t)),1,'first');
                            zdist_euclid(ii_t)=handles_out.zdist_euclid(this_jj);
                            zdist_euclid_within(ii_t)=handles_out.zdist_euclid_within(this_jj);
                            KLdivergence(ii_t)=handles_out.KLdivergence(this_jj);
                            KLdivergence_mix(ii_t)=handles_out.KLdivergence_mix(this_jj);
                            dprime(ii_t)=handles_out.dprime(this_jj);
                            dprime_mix(ii_t)=handles_out.dprime_mix(this_jj);
                            meandFF_sp(ii_t)=this_meandFF_sp(this_jj);
                            meandFF_sm(ii_t)=this_meandFF_sm(this_jj);
                        end
                    end
                end
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).delta_zdist_euclid=(zdist_euclid_within-mean(zdist_euclid_within((trimmed_time_span>-7)&(trimmed_time_span<=-2))))';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).delta_zdist_euclid_within=(zdist_euclid-mean(zdist_euclid((trimmed_time_span>-7)&(trimmed_time_span<=-2))))';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).meandFF_sp=meandFF_sp';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).meandFF_sm=meandFF_sm';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).delta_KLdivergence=(KLdivergence-mean(KLdivergence((trimmed_time_span>-7)&(trimmed_time_span<=-2))))';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).delta_KLdivergence_mix=(KLdivergence_mix-mean(KLdivergence_mix((trimmed_time_span>-7)&(trimmed_time_span<=-2))))';
                delta_dprime=[];
                delta_dprime=(dprime-mean(dprime((trimmed_time_span>-7)&(trimmed_time_span<=-2))))';
                delta_dprime_mix=[];
                delta_dprime_mix=(dprime_mix-mean(dprime_mix((trimmed_time_span>-7)&(trimmed_time_span<=-2))))';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).delta_dprime=delta_dprime';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).delta_dprime_mix=delta_dprime_mix';
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).time_span=handles_out.time_span;
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).lick_sp=handles_out.lick_sp;
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).lick_sm=handles_out.lick_sm;
                all_data.mouse(this_ii_mouse).p_corr(ii_pc_group).movie(ii_movies).time_p_lick=handles_out.time_p_lick;
            end
        end


    end

    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));

    %Plot dprime
    for ii_pc_group=[1 3]
        %Calculate the mean time course for each mouse
        per_mouse_delta_dprime=[];
        per_mouse_delta_dprime_mix=[];
        for ii_mouse=1:length(these_mice)
            ii_movies=all_data.mouse(ii_mouse).p_corr(ii_pc_group).ii_movies;
            if ii_movies>0
                t_points=length(all_data.mouse(ii_mouse).p_corr(ii_pc_group).movie(1).delta_dprime);
                these_delta_dprime=zeros(t_points,ii_movies);
                these_delta_dprime_mix=zeros(t_points,ii_movies);
                for ii_m=1:ii_movies
                    these_delta_dprime(:,ii_m)=all_data.mouse(ii_mouse).p_corr(ii_pc_group).movie(ii_m).delta_dprime;
                    these_delta_dprime_mix(:,ii_m)=all_data.mouse(ii_mouse).p_corr(ii_pc_group).movie(ii_m).delta_dprime_mix;
                end
            end
            per_mouse_delta_dprime(:,ii_mouse)=mean(these_delta_dprime,2);
            per_mouse_delta_dprime_mix(:,ii_mouse)=mean(these_delta_dprime_mix,2);
        end
      
        %delta dprime all in one plot
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        ax=gca;

        set(hFig, 'units','normalized','position',[.2 .2 .4 .25])
        hold on

        try
            CIeu = bootci(1000, @mean, per_mouse_delta_dprime');
            meaneu=mean(per_mouse_delta_dprime',1);
            CIeu(1,:)=meaneu-CIeu(1,:);
            CIeu(2,:)=CIeu(2,:)-meaneu;

            [hlsp, hpsp] = boundedline(trimmed_time_span',meaneu', CIeu', 'cmap',[0/255 0/255 0/255]);
        catch
            plot(trimmed_time_span',meaneu','Color',[150/255 150/255 150/255]);
        end

        try
            CIeu = bootci(1000, @mean, per_mouse_delta_dprime_mix');
            meaneu=mean(per_mouse_delta_dprime_mix',1);
            CIeu(1,:)=meaneu-CIeu(1,:);
            CIeu(2,:)=CIeu(2,:)-meaneu;

            [hlsp, hpsp] = boundedline(trimmed_time_span',meaneu', CIeu', 'cmap',[150/255 150/255 150/255]);
        catch
            plot(trimmed_time_span',meaneu', 'Color',[0/255 0/255 0/255]);
        end

        try
            CIeu = bootci(1000, @mean, per_mouse_delta_dprime'-per_mouse_delta_dprime_mix');
            meaneu=mean(per_mouse_delta_dprime'-per_mouse_delta_dprime_mix',1);
            CIeu(1,:)=meaneu-CIeu(1,:);
            CIeu(2,:)=CIeu(2,:)-meaneu;

            [hlsp, hpsp] = boundedline(trimmed_time_span',meaneu', CIeu', 'cmap',[244/255 165/255 30/255]);
        catch
            plot(trimmed_time_span',meaneu', 'Color',[244/255 165/255 30/255]);
        end

        this_ylim=ylim;

        %Place epoch markers

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

        % legend('Within S+','Within S-', 'Between')

        xlim([-7 15])
        ax.LineWidth=3;
        title(['Delta dprime ' prof_labels{ii_pc_group}])
        xlabel('Time(sec)')
        ylabel('delta dprime, black=between, gray=within, tan=bth-wth')
    end

    %Perform glm and generate bar graphs for delta prime
    %This try is here because of an error in drgCaImAn_batch_traces_pre_pre

    glm_ddp=[];
    glm_ddpbmw=[];
    glm_ii=0;
    glm_iibmw=0;

    id_ii=0;
    id_ii_bmw=0;
    input_data=[];
    input_data_bmw=[];

    %Do the bar graph for dprime
    %I will use the following time bins
    %-5 to -3
    %-2.5 to -1.5
    %-1 to 0
    %2 to 4.1

    figNo = figNo + 1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .5 .25])

    hold on

    figNo = figNo + 1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.3 .3 .5 .25])

    hold on

    edges=[-1:0.5:5];
    rand_offset_s=0.5;
    rand_offset=0.3;


    bar_offset=0;
    bar_offset_bmw=0;

    for ii_pc_group=[1 3]

        %Calculate the mean for each mouse for each time period
        for ii_t_period=1:size(time_periods_eu,1)
            per_mouse_delta_dprime=[];
            per_mouse_delta_dprime_mix=[];
            all_delta_dprime=[];
            all_delta_dprime_mix=[];
            for ii_mouse=1:length(these_mice)
                ii_movies=all_data.mouse(ii_mouse).p_corr(ii_pc_group).ii_movies;
                if ii_movies>0
                    these_delta_dprime=zeros(1,ii_movies);
                    these_delta_dprime_mix=zeros(1,ii_movies);
                    for ii_m=1:ii_movies
                        these_delta_dprime(1,ii_m)=mean(all_data.mouse(ii_mouse).p_corr(ii_pc_group).movie(ii_m).delta_dprime...
                            ((trimmed_time_span>=time_periods_eu(ii_t_period,1))&(trimmed_time_span<=time_periods_eu(ii_t_period,2))));
                        these_delta_dprime_mix(1,ii_m)=mean(all_data.mouse(ii_mouse).p_corr(ii_pc_group).movie(ii_m).delta_dprime_mix...
                            ((trimmed_time_span>=time_periods_eu(ii_t_period,1))&(trimmed_time_span<=time_periods_eu(ii_t_period,2))));
                    end
                end
                per_mouse_delta_dprime(:,ii_mouse)=mean(these_delta_dprime,2);
                per_mouse_delta_dprime_mix(:,ii_mouse)=mean(these_delta_dprime_mix,2);
                all_delta_dprime=[all_delta_dprime these_delta_dprime];
                all_delta_dprime_mix=[all_delta_dprime_mix these_delta_dprime_mix];


                glm_ddp.data(glm_ii+1:glm_ii+length(these_delta_dprime))=these_delta_dprime;
                glm_ddp.pcorr(glm_ii+1:glm_ii+length(these_delta_dprime))=ii_pc_group*ones(1,length(these_delta_dprime));
                glm_ddp.window(glm_ii+1:glm_ii+length(these_delta_dprime))=ii_t_period*ones(1,length(these_delta_dprime));
                glm_ddp.bvsw(glm_ii+1:glm_ii+length(these_delta_dprime))=1*ones(1,length(these_delta_dprime));
                glm_ddp.mice(glm_ii+1:glm_ii+length(these_delta_dprime))=ii_mouse*ones(1,length(these_delta_dprime));
                glm_ii=glm_ii+length(these_delta_dprime);

                glm_ddp.data(glm_ii+1:glm_ii+length(these_delta_dprime_mix))=these_delta_dprime_mix;
                glm_ddp.pcorr(glm_ii+1:glm_ii+length(these_delta_dprime_mix))=ii_pc_group*ones(1,length(these_delta_dprime_mix));
                glm_ddp.window(glm_ii+1:glm_ii+length(these_delta_dprime_mix))=ii_t_period*ones(1,length(these_delta_dprime_mix));
                glm_ddp.bvsw(glm_ii+1:glm_ii+length(these_delta_dprime_mix))=0*ones(1,length(these_delta_dprime_mix));
                glm_ddp.mice(glm_ii+1:glm_ii+length(these_delta_dprime_mix))=ii_mouse*ones(1,length(these_delta_dprime_mix));
                glm_ii=glm_ii+length(these_delta_dprime_mix);

                glm_ddpbmw.data(glm_iibmw+1)=per_mouse_delta_dprime(:,ii_mouse)-per_mouse_delta_dprime_mix(:,ii_mouse);
                glm_ddpbmw.pcorr(glm_iibmw+1)=ii_pc_group;
                glm_ddpbmw.window(glm_iibmw+1)=ii_t_period;
                glm_iibmw=glm_iibmw+1;
            end

            id_ii=id_ii+1;
            input_data(id_ii).data=all_delta_dprime;
            input_data(id_ii).description=[period_labels{ii_t_period} ' between ' prof_labels{ii_pc_group}];

            id_ii=id_ii+1;
            input_data(id_ii).data=all_delta_dprime_mix;
            input_data(id_ii).description=[period_labels{ii_t_period} ' within ' prof_labels{ii_pc_group}];

            id_ii_bmw=id_ii_bmw+1;
            input_data_bmw(id_ii_bmw).data=per_mouse_delta_dprime-per_mouse_delta_dprime_mix;
            input_data_bmw(id_ii_bmw).description=[period_labels{ii_t_period}  prof_labels{ii_pc_group}];

            figure(figNo-1)

            %Bar within
            bar_offset=bar_offset+1;

            if ii_pc_group==1
                hp=bar(bar_offset,mean(all_delta_dprime_mix),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                hatchfill2(hp,'single','HatchAngle',45,'hatchcolor',[150/255 150/255 150/255]);
            else
                hp=bar(bar_offset,mean(all_delta_dprime_mix),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                hatchfill2(hp,'single','HatchAngle',45,'hatchcolor',[150/255 150/255 150/255]);
            end


            if length(all_delta_dprime_mix)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(all_delta_dprime_mix...
                    ,edges,bar_offset,rand_offset_s,'k','k',1);
            end



            %Bar between
            bar_offset=bar_offset+1;

            if ii_pc_group==1
                hp=bar(bar_offset,mean(all_delta_dprime),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                hatchfill2(hp,'single','HatchAngle',45,'hatchcolor',[0/255 0/255 0/255]);
            else
                hp=bar(bar_offset,mean(all_delta_dprime),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                hatchfill2(hp,'single','HatchAngle',45,'hatchcolor',[0/255 0/255 0/255]);
            end

            if length(all_delta_dprime)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(all_delta_dprime...
                    ,edges,bar_offset,rand_offset_s,'k','k',1);
            end

            for ii_mouse=1:length(per_mouse_delta_dprime)
                plot([bar_offset-1 bar_offset],[per_mouse_delta_dprime_mix(ii_mouse) per_mouse_delta_dprime(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
            end

            %Now do between-within
            figure(figNo)
            %Bar
            bar_offset_bmw=bar_offset_bmw+1;

            if ii_pc_group==1
                bar(bar_offset_bmw,mean(per_mouse_delta_dprime-per_mouse_delta_dprime_mix),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            else
                bar(bar_offset_bmw,mean(per_mouse_delta_dprime-per_mouse_delta_dprime_mix),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
            end

            if length(all_delta_dprime_mix)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(per_mouse_delta_dprime-per_mouse_delta_dprime_mix...
                    ,edges,bar_offset_bmw,rand_offset,'k','k',4);
            end







        end
        bar_offset=bar_offset+1;
        bar_offset_bmw=bar_offset_bmw+1;
    end

    figure(figNo-1)

    xticks([1.5 3.5 5.5 7.5 10.5 12.5 14.5 16.5])
    xticklabels({'Baseline','Pre','Odor','Reinf','Baseline','Pre','Odor','Reinf'})
    xtickangle(45)

    ylabel('delta dprime')
    title('delta dprime black=between gray=within)')

    figure(figNo)
    xticks([1 2 3 4 6 7 8 9])
    xticklabels({'Baseline','Pre','Odor','Reinf','Baseline','Pre','Odor','Reinf'})
    xtickangle(45)

    ylabel('delta dprime')

    title('delta dprime (between-within)')



    %Perform the glm for separate within and between
    fprintf(1, ['\nglm for delta dprime including mice\n'])
    fprintf(fileID, ['\nglm for delta dprime including mice\n']);

    tbl = table(glm_ddp.data',glm_ddp.pcorr',glm_ddp.window',glm_ddp.bvsw',glm_ddp.mice',...
        'VariableNames',{'d_dprime','p_correct','time_window','between_vs_within','mice'});
    mdl = fitglm(tbl,'d_dprime~p_correct+time_window+between_vs_within+mice+p_correct*time_window*between_vs_within'...
        ,'CategoricalVars',[2,3,4,5])

%     fprintf(1, ['\nNote: mice did not differ, do glm without mice as categorical variable\n'])
%     fprintf(fileID, ['\nNote: mice did not differ, do glm without mice as categorical variable\n']);
% 
%     fprintf(1, ['\nglm for delta dprime without mice as categorical variable\n'])
%     fprintf(fileID, ['\nglm for delta dprime without mice as categorical variable\n']);
% 
%     tbl = table(glm_ddp.data',glm_ddp.pcorr',glm_ddp.window',glm_ddp.bvsw',...
%         'VariableNames',{'d_dprime','p_correct','time_window','between_vs_within'});
%     mdl = fitglm(tbl,'d_dprime~p_correct+time_window+between_vs_within+p_correct*time_window*between_vs_within'...
%         ,'CategoricalVars',[2,3,4])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);



    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta dprime\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for delta dprime\n']);

    [output_data] = drgMutiRanksumorTtest(input_data, fileID,0);

    %Perform the glm for separate between-within
    fprintf(1, ['\nglm for delta dprime (between-within)\n'])
    fprintf(fileID, ['\nglm for delta dprime (between-within)\n']);

    tbl = table(glm_ddpbmw.data',glm_ddpbmw.pcorr',glm_ddpbmw.window',...
        'VariableNames',{'d_dprime_bmw','p_correct','time_window'});
    mdl = fitglm(tbl,'d_dprime_bmw~p_correct+time_window+p_correct*time_window'...
        ,'CategoricalVars',[2,3])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);



    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for delta dprime (between-within)\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for delta dprime  (between-within)\n']);

    [output_data] = drgMutiRanksumorTtest(input_data_bmw, fileID,0);


 
    pffft=1;
end
fclose(fileID)






