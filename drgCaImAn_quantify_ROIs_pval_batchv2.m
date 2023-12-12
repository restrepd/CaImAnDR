%drgCaImAn_analyze_pval_batch2.m
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

naive_pro{1}='Learning';
naive_pro{3}='Proficient';


if isfield(handles_in,'min_tr_div')
    min_tr_div=handles_in.min_tr_div;
    min_tr_resp=handles_in.min_tr_resp;
else
    min_tr_div=12;
    min_tr_resp=6;
end



%Find out which files were included (>=min_tr trials)
files_included_glm=zeros(1,length(handles_out.file));
for fNo=1:length(handles_out.file)
    if handles_out.file(fNo).output_data_odor.total_trials(1)>=min_tr_div
        files_included_glm(fNo)=1;
    end
end
files_included_glm=logical(files_included_glm);


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
for mouseNo=1:max((handles_out.mouseNo))
    for grNo=1:max(groups)
        ROI_accounting.group(grNo).mouse(mouseNo).no_sessions=0;
    end
end
for grNo=groups
    for mouseNo=unique(handles_out.mouseNo)

        switch grNo
            case 1
                files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo)&files_included_glm;
            case 2
                files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
            case 3
                files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo)&files_included_glm;
        end

        %Do ROI accounting per mouse
        for fileNo=1:length(files_included)
            if files_included(fileNo)==1
                ROI_accounting.group(grNo).mouse(mouseNo).no_sessions=ROI_accounting.group(grNo).mouse(mouseNo).no_sessions+1;
                ROI_accounting.group(grNo).mouse(mouseNo).session(ROI_accounting.group(grNo).mouse(mouseNo).no_sessions).no_ROIS=handles_out.file(fileNo).no_neurons;
            end
        end
    end


end

%Now print the output (ROI+/-SD, no sessions)
for grNo=groups
    all_mice_no_ROIS=[];
    for mouseNo=unique(handles_out.mouseNo)
        noROIs=zeros(1,ROI_accounting.group(grNo).mouse(mouseNo).no_sessions);
        if ROI_accounting.group(grNo).mouse(mouseNo).no_sessions>0
            for sessionNo=1:ROI_accounting.group(grNo).mouse(mouseNo).no_sessions
                noROIs(sessionNo)=ROI_accounting.group(grNo).mouse(mouseNo).session(sessionNo).no_ROIS;
            end
            all_mice_no_ROIS=[all_mice_no_ROIS noROIs];
            fprintf(1,['For ' naive_pro{grNo} ' mouse number ' num2str(mouseNo) ' mean number of ROIs ' num2str(mean(noROIs)) ' STD ' num2str(std(noROIs)) ' n= ' num2str(length(noROIs)) '\n'])
        else
            fprintf(1,['For ' naive_pro{grNo} ' mouse number ' num2str(mouseNo) ' no sessions\n'])

        end

    end
    fprintf(1,['For ' naive_pro{grNo} ' all mice mean number of ROIs ' num2str(mean(all_mice_no_ROIS)) ' STD ' num2str(std(all_mice_no_ROIS)) ' n= ' num2str(length(all_mice_no_ROIS)) '\n'])
end
