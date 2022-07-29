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


%Now plot for each algorithm the prediction accuracy
handles_out2=[];

for grNo=1:max(handles.group)
    handles_out2.group_no(grNo).ii_files=0;
end

ii_out=1;

if exist([choiceBatchPathName choiceFileName(1:end-2) '.mat'])==0
    %Process each file separately
    
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

%Bar graph plot for accuracy
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