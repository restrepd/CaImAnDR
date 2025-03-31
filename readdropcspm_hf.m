function rspm_out=readdropcspm_hf(PathName,FileName,show_figure)
% readdropcspm_hf
% This script file reads Splus_minus output files
% The program assumes version 1.1 output files
% First, open the file

%VERY IMPORTANT: When you bring in the files sort them by date

if nargin==0
    clear;
    close all;
    show_figure=1;
    [FileName,PathName] = uigetfile('*spm.mat','Select the go-no go file','MultiSelect','On');
end


trialNo=0;
if ischar(FileName)==1
    no_files=1;
else
    no_files=length(FileName);
end
for fileNum=1:no_files
    if ischar(FileName)==1
        fullName=[PathName,FileName];
        load(fullName)
    else
       handles=[];
       fullName=[PathName,FileName{fileNum}];
       load(fullName) 
    end
    for trNo=1:length(handles.dropcData.trialScore)
        trialNo=trialNo+1;
        handlesin.dropcData.trialScore(trialNo)=handles.dropcData.trialScore(trNo);
        handlesin.dropcData.isReinforced(trialNo)=handles.dropcData.isReinforced(trNo);
        handlesin.dropcData.odorType(trialNo)=handles.dropcData.odorType(trNo);
    end
    handlesin.dropcData.trials_per_file(fileNum)=length(handles.dropcData.trialScore);
end


%% Calculate percent correct for the sliding window
no_trials=trialNo;
sliding_window=20; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;


score=~(handlesin.dropcData.trialScore==(handlesin.dropcData.odorType-1));

 for ii=1:no_trials-sliding_window+1
%         first_time=handlesin.dropcData.trialTime(ii);
%         last_time=handlesin.dropcData.trialTime(ii+sliding_window-1);
        rspm_out.perCorr(ii+(sliding_window/2))=100*sum(score(ii:ii+sliding_window-1))/sliding_window;
        
        if ii==1
            rspm_out.perCorr(ii:(sliding_window/2))=rspm_out.perCorr(ii+(sliding_window/2));
        end
        if ii==no_trials-sliding_window+1
            rspm_out.perCorr(ii+(sliding_window/2)+1:no_trials)=rspm_out.perCorr(ii+(sliding_window/2));
        end
 end
    
% Plot percent correct vs trial 
if show_figure==1
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.25 .25 .5 .25])
    
    jj_low=find(rspm_out.perCorr<max_percent_low_beh);
    plot(jj_low,rspm_out.perCorr(jj_low),'ob')
    hold on
    jj_high=find(rspm_out.perCorr>min_precent_high_beh);
    plot(jj_high,rspm_out.perCorr(jj_high),'or')
    
    jj_mid=find((rspm_out.perCorr<=min_precent_high_beh)&(rspm_out.perCorr>=max_percent_low_beh));
    plot(jj_mid,rspm_out.perCorr(jj_mid),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
    hold on
    plot([1 no_trials],[50 50],'-k')
    
    for filN=1:fileNum
        plot([sum(handlesin.dropcData.trials_per_file(1:filN)) sum(handlesin.dropcData.trials_per_file(1:filN)) ],[0 100],'-k')
    end
    
    if ischar(FileName)==1
        title(['Percent correct vs. trial number ' FileName])
    else
        title(['Percent correct vs. trial number ' FileName{1}])
    end
    xlabel('Trial number')
    ylabel('Percent correct')
    ylim([0 100])
end

%Calculate the per file mean percent correct
ii_last=0;
for filN=1:fileNum
    rspm_out.percent_correct_per_file(filN)=100*sum(score(ii_last+1:ii_last+handlesin.dropcData.trials_per_file(filN)))/handlesin.dropcData.trials_per_file(filN);
    ii_last=ii_last+handlesin.dropcData.trials_per_file(filN);
end

pffft=1

% %% Now do intertrial intervals
% 
% rspm_out.total_time_min=handlesin.dropcData.epochTime(end)/60;
% trialNo=length(handlesin.dropcData.trialTime);
% rspm_out.iti=handlesin.dropcData.trialTime(2:trialNo)-handlesin.dropcData.trialTime(1:trialNo-1);
% rspm_out.mean_iti=mean(rspm_out.iti);
% rspm_out.mean_iti_maxPer=mean(rspm_out.iti(min_ii:max_ii));
% rspm_out.median_iti_maxPer=median(rspm_out.iti(min_ii:max_ii))
% 
% if show_figure==1
%     try
%         close 2
%     catch
%     end
%     
%     hFig2 = figure(2);
%     set(hFig2, 'units','normalized','position',[.25 .45 .5 .25])
%     
%     
%     plot(1:trialNo-1,rspm_out.iti)
%     
%     
%     title(['Intertrial interval vs. trial number' FileName])
%     xlabel('Trial number')
%     ylabel('Intertrial interval (sec)')
%     xlim([0 160])
%     ylim([0 160])
% end

pfft=1


