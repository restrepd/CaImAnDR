%drgCaImAn_analyze_batch_pre_per_dim_per_mousev2
close all force
clear all



[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dim_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdzv2 run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

fileID = fopen([choiceBatchPathName 'drgCaImAn_analyze_batch_pre_per_dim_per_mouse.txt'],'w');

no_files=handles.no_files;
moving_mean_n=10;


handles_out2=[];

for grNo=1:max(handles.group)
    handles_out2.group_no(grNo).ii_files=0;
end

ii_out=1;

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

t_from=-10;
t_to=20;

fr_per_names{1}='<45%%';
fr_per_names{2}='45-65%%';
fr_per_names{3}='65-80%%';
fr_per_names{4}='>=80%%';
% fr_per_names{5}='Rev <40%%';
% fr_per_names{6}='Rev 40-65%%';
% fr_per_names{7}='Rev 65-80%%';
% fr_per_names{8}='Rev >=80%%';

per_names{1}='<40%';
per_names{2}='40-65%';
per_names{3}='65-80%';
per_names{4}='>=80%';

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


these_groups=unique(handles.group);

if exist([choiceBatchPathName choiceFileName(1:end-2) 'v2.mat'])==0
    %Process each file separately

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



    %Find the mean ROI number
    all_noROIs=[];
    for fileNo=1:no_files
        tic
        %         pre_per_PathName=handles.PathName_pre_per{fileNo};

        if iscell(handles.PathName_out)
            pre_per_outPathName=handles.PathName_out{fileNo};
        else
            pre_per_outPathName=handles.PathName_out;
        end


        pre_per_FileName=handles.FileName_pre_per{fileNo};


        load([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])

        this_noROIs=handles_out.ii_out(ii_out).handles_out.noROIs_before_trimming;
        all_noROIs=[all_noROIs this_noROIs];

        fprintf(1, ['Number of ROIs for file No ' num2str(fileNo) ' is ' num2str(this_noROIs) '\n'])
    end

    mean_noROIs=mean(all_noROIs);
    handles_out2.mean_noROIs=mean_noROIs;
    fprintf(1, ['Mean number of ROIs ' num2str(mean_noROIs)  '\n'])

    for grNo=1:no_pcorr*length(these_groups)
        handles_out2.group_no(grNo).ii_files=0;
    end

    allPcorr=[];
    allShort=[];
    for fileNo=1:no_files
        tic
        %         pre_per_PathName=handles.PathName_pre_per{fileNo};
        if iscell(handles.PathName_out)
            pre_per_outPathName=handles.PathName_out{fileNo};
        else
            pre_per_outPathName=handles.PathName_out;
        end

        pre_per_FileName=handles.FileName_pre_per{fileNo};


        load([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])

        this_group=handles.group(fileNo);
        this_grNo=find(these_groups==this_group);

        if isfield(handles_out.ii_out(1).handles_out,'percent_correct')
            pCorr=handles_out.ii_out(1).handles_out.percent_correct;
            allPcorr=[allPcorr pCorr];
            handles_out2.pcorr_per_file(fileNo)=pCorr;
            ii_pCorr=1;
            if (pCorr>=45)&(pCorr<=65)
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

        grNo=ii_pCorr;



        handles_out2.group_no(grNo).ii_files=handles_out2.group_no(grNo).ii_files+1;
        ii_files=handles_out2.group_no(grNo).ii_files;
        if ~isfield(handles_out2.group_no(grNo),'ii_time_span')
            handles_out2.group_no(grNo).ii_time_span=[];
            handles_out2.group_no(grNo).time_span=[];
        end

        %Note that dimensionality is normalized
        handles_out2.group_no(grNo).noROIs_before_trimming(ii_files)=handles_out.ii_out(ii_out).handles_out.noROIs_before_trimming;
        %         handles_out2.group_no(grNo).this_noROIs=handles_out.ii_out(ii_out).handles_out.noROIs_before_trimming;
        handles_out2.group_no(grNo).ii_time_span(ii_files)=length(handles_out.ii_out(ii_out).handles_out.time_span);
        if length(handles_out.ii_out(1).handles_out.time_span)<700
            allShort=[allShort fileNo];
            pfft=1;
        end
        handles_out2.group_no(grNo).time_span(ii_files,1:length(handles_out.ii_out(ii_out).handles_out.time_span))=handles_out.ii_out(ii_out).handles_out.time_span;
        handles_out2.group_no(grNo).dimensionality(ii_files,1:length(handles_out.ii_out(ii_out).handles_out.dimensionality))=mean_noROIs*handles_out.ii_out(ii_out).handles_out.dimensionality/this_noROIs;
        handles_out2.group_no(grNo).dimensionalitysp(ii_files,1:length(handles_out.ii_out(ii_out).handles_out.dimensionality))=mean_noROIs*handles_out.ii_out(ii_out).handles_out.dimensionalitysp/this_noROIs;
        handles_out2.group_no(grNo).dimensionalitysm(ii_files,1:length(handles_out.ii_out(ii_out).handles_out.dimensionality))=mean_noROIs*handles_out.ii_out(ii_out).handles_out.dimensionalitysm/this_noROIs;
        handles_out2.group_no(grNo).mouseNo(ii_files)=handles_out2.mouseNo_per_file(fileNo);


        fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
    end
else
    load([choiceBatchPathName choiceFileName(1:end-2) '.mat'])
end


figureNo=0;


%Plot the timecourse for dimensionality for each group per mouse
per_mouse_dim=[];

%Note: These data were acquired at different rates. Here they are all
%resampled to a dt of 0.03 sec
dt_res=0.03;
time_span=t_from:dt_res:t_to;
fprintf(1, ['Length of time_span ' num2str(length(time_span)) '\n'])

for mouseNo=1:length(handles_out2.mouse_names)
    for grNo=[2 4]

        %         per_mouse_dim.group(grNo).mouse(mouseNo).dimensionality=[];
        %         per_mouse_dim.group(grNo).mouse(mouseNo).time_span=[];

        if ~isempty(handles_out2.group_no(grNo).ii_time_span)

            these_dim=[];
            these_dim_sp=[];
            these_dim_sm=[];


            %Extrapolate all points onto time_span
            these_mice=handles_out2.group_no(grNo).mouseNo;
            ii_included=0;
            for ii_f=1:length(handles_out2.group_no(grNo).ii_time_span)
                print_out=1;

                if (these_mice(ii_f)==mouseNo)
                    ii_included=ii_included+1;

                    this_ii_tspan=handles_out2.group_no(grNo).ii_time_span(ii_f);
                    this_time_span=handles_out2.group_no(grNo).time_span(ii_f,1:this_ii_tspan);

                    for ii_tsp=1:length(time_span)
                        if time_span(ii_tsp)<this_time_span(1)
                            these_dim(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionality(ii_f,1);
                            these_dim_sp(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionalitysp(ii_f,1);
                            these_dim_sm(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionalitysm(ii_f,1);
                        else
                            if time_span(ii_tsp)>this_time_span(end)
                                these_dim(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionality(ii_f,end);
                                these_dim_sp(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionalitysp(ii_f,end);
                                these_dim_sm(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionalitysm(ii_f,end);
                                if print_out==1
                                    fprintf(1, ['Mouse No ' num2str(mouseNo) ' Group No ' num2str(grNo) ' File No ' num2str(ii_f) ' ii_tsp ' num2str(ii_tsp) '\n'])
                                    print_out=0;
                                end
                            else
                                ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                                ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                                these_dim(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionality(ii_f,ii_0)+...
                                    (handles_out2.group_no(grNo).dimensionality(ii_f,ii_1)-handles_out2.group_no(grNo).dimensionality(ii_f,ii_0))...
                                    *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                these_dim_sp(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionalitysp(ii_f,ii_0)+...
                                    (handles_out2.group_no(grNo).dimensionalitysp(ii_f,ii_1)-handles_out2.group_no(grNo).dimensionalitysp(ii_f,ii_0))...
                                    *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                these_dim_sm(ii_included,ii_tsp)=handles_out2.group_no(grNo).dimensionalitysm(ii_f,ii_0)+...
                                    (handles_out2.group_no(grNo).dimensionalitysm(ii_f,ii_1)-handles_out2.group_no(grNo).dimensionalitysm(ii_f,ii_0))...
                                    *(time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                            end
                        end
                    end


                end
            end


            %For some reason whe have a subset of sessions with zero accuracy
            vetted_these_dim=[];
            ii_s=0;
            for ii_sessions=1:size(these_dim,1)
                if sum(these_dim(ii_sessions,:)==0)<200
                    ii_s=ii_s+1;
                    vetted_these_dim(ii_s,:)=these_dim(ii_sessions,:);
                end
            end

            these_dim=[];
            these_dim=vetted_these_dim;

            vetted_these_dim=[];
            ii_s=0;
            for ii_sessions=1:size(these_dim_sp,1)
                if sum(these_dim_sp(ii_sessions,:)==0)<200
                    ii_s=ii_s+1;
                    vetted_these_dim(ii_s,:)=these_dim_sp(ii_sessions,:);
                end
            end

            these_dim_sp=[];
            these_dim_sp=vetted_these_dim;

            vetted_these_dim=[];
            ii_s=0;
            for ii_sessions=1:size(these_dim_sm,1)
                if sum(these_dim_sm(ii_sessions,:)==0)<200
                    ii_s=ii_s+1;
                    vetted_these_dim(ii_s,:)=these_dim_sm(ii_sessions,:);
                end
            end

            these_dim_sm=[];
            these_dim_sm=vetted_these_dim;

 
            per_mouse_dim.group(grNo).mouse(mouseNo).dimensionality=these_dim;
            per_mouse_dim.group(grNo).mouse(mouseNo).mean_dim=mean(these_dim,1)';
            per_mouse_dim.group(grNo).mouse(mouseNo).mean_dim_sp=mean(these_dim_sp,1)';
            per_mouse_dim.group(grNo).mouse(mouseNo).mean_dim_sm=mean(these_dim_sm,1)';
            per_mouse_dim.group(grNo).mouse(mouseNo).time_span=time_span';



        end


    end
end



%Plot the dimensionality calculated per mouse
dim_per_group=[];

glm_dim=[];
glm_ii=0;

id_ii=0;
input_data=[];

these_mean_dim_proficient=[];
ii_proficient=0;

for grNo=[2 4]

    ii_m_included=0;
    these_mice=zeros(1,100);
    these_mean_dim_per_mouse=[];
    ii_d_included=0;
    these_dimensionality=zeros(100,length(time_span));

    for mouseNo=1:length(handles_out2.mouse_names)
        if ~isempty(per_mouse_dim.group(grNo).mouse(mouseNo).dimensionality)
            ii_m_included=ii_m_included+1;
            these_mean_dim_per_mouse(ii_m_included,:)=per_mouse_dim.group(grNo).mouse(mouseNo).mean_dim;
            delta_ii=size(per_mouse_dim.group(grNo).mouse(mouseNo).dimensionality,1);
            these_dimensionality(ii_d_included+1:ii_d_included+delta_ii,:)=per_mouse_dim.group(grNo).mouse(mouseNo).dimensionality;
            these_mice(ii_d_included+1:ii_d_included+delta_ii)=mouseNo;
            ii_d_included=ii_d_included+delta_ii;
        end
    end

    these_dimensionality(ii_d_included+1:end,:)=[];
    these_mice=these_mice(1:ii_d_included);

    %Save the data for the bar graph
    %Pre window
    these_dim_per_mouse=mean(these_mean_dim_per_mouse(:,(time_span>=-1)&(time_span<=0)),2);
    these_dim_per_session=mean(these_dimensionality(:,(time_span>=-1)&(time_span<=0)),2);

    glm_dim.data(glm_ii+1:glm_ii+length(these_dim_per_session))=these_dim_per_session;
%     if grNo<5
%         glm_dim.fwd_rev(glm_ii+1:glm_ii+length(these_dim_per_session))=0*ones(1,length(these_dim_per_session));
%     else
%         glm_dim.fwd_rev(glm_ii+1:glm_ii+length(these_dim_per_session))=1*ones(1,length(these_dim_per_session));
%     end
%     pcorr=rem(grNo-1,4)+1;
    glm_dim.pcorr(glm_ii+1:glm_ii+length(these_dim_per_session))=grNo*ones(1,length(these_dim_per_session));
    glm_dim.window(glm_ii+1:glm_ii+length(these_dim_per_session))=1*ones(1,length(these_dim_per_session));
    glm_dim.mice(glm_ii+1:glm_ii+length(these_dim_per_session))=these_mice;

    glm_ii=glm_ii+length(these_dim_per_session);


    dim_per_group.group(grNo).window(1).dims_per_session=these_dim_per_session;
    dim_per_group.group(grNo).window(1).dims_per_mouse=these_dim_per_mouse;

    id_ii=id_ii+1;
    input_data(id_ii).data=these_dim_per_session;
    input_data(id_ii).description=['Pre ' fr_per_names{grNo}];

    %Odor window
    these_dim_per_mouse=mean(these_mean_dim_per_mouse(:,(time_span>=3.1)&(time_span<=4.1)),2);
    these_dim_per_session=mean(these_dimensionality(:,(time_span>=3.1)&(time_span<=4.1)),2);

    glm_dim.data(glm_ii+1:glm_ii+length(these_dim_per_session))=these_dim_per_session;
%     if grNo<5
%         glm_dim.fwd_rev(glm_ii+1:glm_ii+length(these_dim_per_session))=0*ones(1,length(these_dim_per_session));
%     else
%         glm_dim.fwd_rev(glm_ii+1:glm_ii+length(these_dim_per_session))=1*ones(1,length(these_dim_per_session));
%     end
%     pcorr=rem(grNo-1,4)+1;
    glm_dim.pcorr(glm_ii+1:glm_ii+length(these_dim_per_session))=grNo*ones(1,length(these_dim_per_session));
    glm_dim.window(glm_ii+1:glm_ii+length(these_dim_per_session))=2*ones(1,length(these_dim_per_session));
    glm_dim.mice(glm_ii+1:glm_ii+length(these_dim_per_session))=these_mice;

    glm_ii=glm_ii+length(these_dim_per_session);

    dim_per_group.group(grNo).window(2).dims_per_session=these_dim_per_session;
    dim_per_group.group(grNo).window(2).dims_per_mouse=these_dim_per_mouse;

    id_ii=id_ii+1;
    input_data(id_ii).data=these_dim_per_session;
    input_data(id_ii).description=['Odor ' fr_per_names{grNo}];

    if size(these_mean_dim_per_mouse,1)>2


        figureNo = figureNo + 1;
        try
            close(figureNo)
        catch
        end
        hFig=figure(figureNo);
        hold on

        ax=gca;ax.LineWidth=3;
        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])


        CIpv = bootci(1000, @mean, these_dimensionality);
        meanpv=mean(these_dimensionality,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(time_span',mean(these_dimensionality,1)', CIpv', 'k');

%         for ii_session=1:size(these_mean_dim_per_mouse,1)
%             plot(time_span',smoothdata(these_mean_dim_per_mouse(ii_session,:)','gaussian',100),'Color',[100/255 100/255 100/255],'LineWidth',1)
%             if (grNo==4)||(grNo==8)
%                 ii_proficient=ii_proficient+1;
%                 these_mean_dim_per_mouse_proficient(ii_proficient,:)=these_mean_dim_per_mouse(ii_session,:);
%                 time_span_proficient=time_span;
%             end
%         end
% 
%         plot(time_span',mean(these_mean_dim_per_mouse,1)', 'k','LineWidth',1.5);

        %         ylim([0.3 1.2])
        ylim([3 8])
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

        title(['Dimensionality calculated per mouse for ' per_names{ii_pcorr}])

    end

end

%Plot a bar graph of dimensionality
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.3 .3 .5 .25])

hold on

edges=[0:0.5:12];
rand_offset=0.4;

bar_offset=0;


for grNo=[2 4]



    %Pre
    bar_offset=bar_offset+1;
    these_dim_pre=dim_per_group.group(grNo).window(1).dims_per_session;
    if grNo==2
        bar(bar_offset,mean(these_dim_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
    else
        bar(bar_offset,mean(these_dim_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    end

    if length(these_dim_pre)>2
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_dim_pre...
            ,edges,bar_offset,rand_offset,'k','k',3);
    end


    %Odor
    bar_offset=bar_offset+1;
    these_dim_odor=dim_per_group.group(grNo).window(2).dims_per_session;
    if grNo==2
        bar(bar_offset,mean(these_dim_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
    else
        bar(bar_offset,mean(these_dim_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
    end

    if length(these_dim_odor)>2
        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_dim_odor...
            ,edges,bar_offset,rand_offset,'k','k',3);
    end


     for ii_mouse=1:length(dim_per_group.group(grNo).window(2).dims_per_mouse)
                plot([bar_offset-1 bar_offset],[dim_per_group.group(grNo).window(1).dims_per_mouse(ii_mouse) dim_per_group.group(grNo).window(2).dims_per_mouse(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',3 )
     end

    bar_offset=bar_offset+1;
end

xticks([1 2 4 5])
xticklabels({'Pre','Odor','Pre','Odor'})

title(['Dimensionality'])
ylabel('Dimensionality')

xlim([0 6])


%Perform the glm 
fprintf(1, ['\nglm for dimensionality\n'])
fprintf(fileID, ['\nglm for dimensionality\n']);

tbl = table(glm_dim.data',glm_dim.pcorr',glm_dim.window',glm_dim.mice',...
    'VariableNames',{'dimensionality','percent_correct','window','mice'});
mdl = fitglm(tbl,'dimensionality~percent_correct+window+mice+percent_correct*window'...
    ,'CategoricalVars',[2,3,4])


txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);

fclose(fileID);


out_file=[choiceBatchPathName choiceFileName];
out_file=[out_file(1:end-2) 'v2.mat'];
save(out_file,'handles_out2','handles','-v7.3')
pfft=1;