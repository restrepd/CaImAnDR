%drgCaImAn_analyze_batch_pre_per_dimensionality_fsdzv2
close all
clear all



[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dim_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdzv2 run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

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

fr_per_names{1}='Fwd <40%%';
fr_per_names{2}='Fwd 40-65%%';
fr_per_names{3}='Fwd 65-80%%';
fr_per_names{4}='Fwd >=80%%';
fr_per_names{5}='Rev <40%%';
fr_per_names{6}='Rev 40-65%%';
fr_per_names{7}='Rev 65-80%%';
fr_per_names{8}='Rev >=80%%';

if exist([choiceBatchPathName choiceFileName(1:end-2) '.mat'])==0
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

      for fileNo=1:no_files
          tic
          %         pre_per_PathName=handles.PathName_pre_per{fileNo};
          if iscell(handles.PathName_out)
              pre_per_outPathName=handles.PathName_out{fileNo};
          else
              pre_per_outPathName=handles.PathName_out;
          end

          pre_per_FileName=handles.FileName_pre_per{fileNo};
        grNo=handles.group(fileNo);

        load([pre_per_outPathName pre_per_FileName(1:end-4) handles.suffix_out])

        handles_out2.group_no(grNo).ii_files=handles_out2.group_no(grNo).ii_files+1;
        ii_files=handles_out2.group_no(grNo).ii_files;

        %Note that dimensionality is normalized
        handles_out2.group_no(grNo).file(ii_files).noROIs_before_trimming=handles_out.ii_out(ii_out).handles_out.noROIs_before_trimming;
        this_noROIs=handles_out.ii_out(ii_out).handles_out.noROIs_before_trimming;
        handles_out2.group_no(grNo).file(ii_files).dimensionality=mean_noROIs*handles_out.ii_out(ii_out).handles_out.dimensionality/this_noROIs;
        handles_out2.group_no(grNo).file(ii_files).dimensionalitysp=mean_noROIs*handles_out.ii_out(ii_out).handles_out.dimensionalitysp/this_noROIs;
        handles_out2.group_no(grNo).file(ii_files).dimensionalitysm=mean_noROIs*handles_out.ii_out(ii_out).handles_out.dimensionalitysm/this_noROIs;
        handles_out2.group_no(grNo).file(ii_files).time_span=handles_out.ii_out(ii_out).handles_out.time_span;

        fprintf(1, ['Import for file No ' num2str(fileNo) ' done in ' num2str(toc) ' sec\n'])
    end
else
    load([choiceBatchPathName choiceFileName(1:end-2) '.mat'])
end


figureNo=0;


%Plot the timecourse for dimensionality
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

for grNo=1:max(handles.group)
    
    no_files=handles_out2.group_no(grNo).ii_files;
    dim_length=0;
    for ii_files=1:no_files
        dim_length=max([dim_length length(handles_out2.group_no(grNo).file(ii_files).dimensionality)]);
    end

    this_dimensionality=zeros(no_files,dim_length);
    this_time_span=zeros(1,dim_length);

    for ii_files=1:no_files
        this_dimensionality(ii_files,1:length(handles_out2.group_no(grNo).file(ii_files).dimensionality))=handles_out2.group_no(grNo).file(ii_files).dimensionality;
        this_time_span(1,1:length(handles_out2.group_no(grNo).file(ii_files).time_span))=handles_out2.group_no(grNo).file(ii_files).time_span;
    end

    subplot(max(handles.group),1,grNo)
    
    if no_files>2

        CIpv = bootci(1000, @mean, this_dimensionality);
        meanpv=mean(this_dimensionality,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(this_time_span',mean(this_dimensionality,1)', CIpv', 'b');
    else
        if no_files>0
            plot(this_time_span',mean(this_dimensionality,1)', 'b');
        end

    end
    xlabel('Time(sec)')
    title(handles.group_names{grNo})

    overall_mean=mean(mean(this_dimensionality,2));
    fprintf(1, ['\nMean dimensionality for  ' handles.group_names{grNo} ' ' num2str(overall_mean) '\n']);
end
sgtitle('Dimensionality')


%Plot the timecourse for dimensionality sp
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

for grNo=1:max(handles.group)
    
    no_files=handles_out2.group_no(grNo).ii_files;
    dim_length=0;
    for ii_files=1:no_files
        dim_length=max([dim_length length(handles_out2.group_no(grNo).file(ii_files).dimensionality)]);
    end

    this_dimensionality=zeros(no_files,dim_length);
    this_time_span=zeros(1,dim_length);
 
    for ii_files=1:no_files 
        this_dimensionality(ii_files,1:length(handles_out2.group_no(grNo).file(ii_files).dimensionality))=handles_out2.group_no(grNo).file(ii_files).dimensionalitysp;
        this_time_span(1,1:length(handles_out2.group_no(grNo).file(ii_files).time_span))=handles_out2.group_no(grNo).file(ii_files).time_span;
    end

    subplot(max(handles.group),1,grNo)
    
    if no_files>2

        CIpv = bootci(1000, @mean, this_dimensionality);
        meanpv=mean(this_dimensionality,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(this_time_span',mean(this_dimensionality,1)', CIpv', 'b');
    else
        if no_files>0
            plot(this_time_span',mean(this_dimensionality,1)', 'b');
        end

    end
    xlabel('Time(sec)')
    title(handles.group_names{grNo})

    overall_mean=mean(mean(this_dimensionality,2));
    fprintf(1, ['\nMean dimensionality for  ' handles.group_names{grNo} ' ' num2str(overall_mean) '\n']);
end
sgtitle('Dimensionality S+')



%Plot the timecourse for dimensionality sp
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

for grNo=1:max(handles.group)
    
    no_files=handles_out2.group_no(grNo).ii_files;
    dim_length=0;
    for ii_files=1:no_files
        dim_length=max([dim_length length(handles_out2.group_no(grNo).file(ii_files).dimensionality)]);
    end

    this_dimensionality=zeros(no_files,dim_length);
    this_time_span=zeros(1,dim_length);

    for ii_files=1:no_files
        this_dimensionality(ii_files,1:length(handles_out2.group_no(grNo).file(ii_files).dimensionality))=handles_out2.group_no(grNo).file(ii_files).dimensionalitysm;
        this_time_span(1,1:length(handles_out2.group_no(grNo).file(ii_files).time_span))=handles_out2.group_no(grNo).file(ii_files).time_span;
    end

    subplot(max(handles.group),1,grNo)
    
    if no_files>2

        CIpv = bootci(1000, @mean, this_dimensionality);
        meanpv=mean(this_dimensionality,1);
        CIpv(1,:)=meanpv-CIpv(1,:);
        CIpv(2,:)=CIpv(2,:)-meanpv;


        [hlpvl, hppvl] = boundedline(this_time_span',mean(this_dimensionality,1)', CIpv', 'b');
    else
        if no_files>0
            plot(this_time_span',mean(this_dimensionality,1)', 'b');
        end

    end
    xlabel('Time(sec)')
    title(handles.group_names{grNo})

    overall_mean=mean(mean(this_dimensionality,2));
    fprintf(1, ['\nMean dimensionality for  ' handles.group_names{grNo} ' ' num2str(overall_mean) '\n']);
end
sgtitle('Dimensionality S-')



%Plot the timecourse for dimensionality per file
dy=15;
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .6])

hold on

for grNo=1:max(handles.group)
    
    no_files=handles_out2.group_no(grNo).ii_files;
    dim_length=0;
    for ii_files=1:no_files
        dim_length=max([dim_length length(handles_out2.group_no(grNo).file(ii_files).dimensionality)]);
    end

    this_dimensionality=zeros(no_files,dim_length);
    this_time_span=zeros(1,dim_length);

    for ii_files=1:no_files
        this_dimensionality=[];
        this_time_span=[];
        this_dimensionality(1,1:length(handles_out2.group_no(grNo).file(ii_files).dimensionality))=handles_out2.group_no(grNo).file(ii_files).dimensionality;
        this_time_span(1,1:length(handles_out2.group_no(grNo).file(ii_files).time_span))=handles_out2.group_no(grNo).file(ii_files).time_span;
        plot(this_time_span,this_dimensionality+dy*ii_files,'k')
    end
% 
%   
%     this_time_span=this_time_span-30;
%     if no_files>2
% 
%         CIpv = bootci(1000, @mean, this_dimensionality);
%         meanpv=mean(this_dimensionality,1);
%         CIpv(1,:)=meanpv-CIpv(1,:);
%         CIpv(2,:)=CIpv(2,:)-meanpv;
% 
% 
%         [hlpvl, hppvl] = boundedline(this_time_span',mean(this_dimensionality,1)', CIpv', 'b');
%     else
%         if no_files>0
%             plot(this_time_span',mean(this_dimensionality,1)', 'b');
%         end
% 
%     end
   
end
xlabel('Time(sec)')
title('Dimensionality')

out_file=[choiceBatchPathName choiceFileName];
out_file=[out_file(1:end-2) '.mat'];
save(out_file,'handles_out2','handles','-v7.3')
pfft=1;