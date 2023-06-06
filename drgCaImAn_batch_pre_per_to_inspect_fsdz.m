function drgCaImAn_batch_pre_per_to_inspect_fsdz(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_inspect_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_inspect_fsdz run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

%Parallel batch processing for each file
all_files_present=1;
first_file=handles.first_file;
for filNum=first_file:handles.no_files
     
    
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





if all_files_present==1
    
    
    %Process each file separately
    for fileNo=first_file:length(handles.FileName_pre_per)
        tic
        first_toc=toc;

        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};

        this_handles_choices.pre_perFileName=pre_per_FileName;
        this_handles_choices.pre_per_PathName=pre_per_PathName;


        this_handles_choices.post_time=5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on
        this_handles_choices.post_shift=0; %Set to 0 if you want to train with odor on points
        this_handles_choices.pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift

        this_handles_choices.p_threshold=1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
        this_handles_choices.dt_p_threshold=4.2; %Time to be used after the odor on for the p_threshold t_test
        this_handles_choices.show_figures=1; %Show the figures
        this_handles_choices.min_trials=20; %Minimum number of trials needed to perform calculations for each percent group
        this_handles_choices.conv_dt=0.3;

        start_toc=toc;


%         drgCaImAn_inspect_traces_pre_pre(this_handles_choices);
        drgCaImAn_inspect_traces_pre_prev2(this_handles_choices);


        fprintf(1, ['Data processed for file number %d took %d minutes\n'],fileNo,(toc-start_toc)/60);


    end
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
     
end






