function drgCaImAn_batch_pre_per_to_decode_entire_session_fsdz(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

    

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_LDA_fsdz run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

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


% if exist([handles.PathName_out handles.FileName_out])==0
%     handles_out=[];
%     ii_out=0;
% else
%     load([handles.PathName_out handles.FileName_out])
%     ii_out=handles_out.last_ii_out;
%     first_file=handles_out.last_file_processed+1;
% end

figNo=0;
show_figures=1;

if all_files_present==1
    
    all_delta_odor=[];
    all_delta_odor_on_reinf_on=[];
    all_delta_reinf=[];

    %Process each file separately
    for fileNo=first_file:length(handles.FileName_pre_per)
        tic
        first_toc=toc;

        handles_out=[];
        ii_out=0;
        
        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};

        this_handles=handles;

        load([pre_per_PathName pre_per_FileName])

        all_delta_odor=[all_delta_odor delta_odor];
        all_delta_odor_on_reinf_on=[all_delta_odor_on_reinf_on delta_odor_on_reinf_on];
        all_delta_reinf=[all_delta_reinf delta_reinf];

        handles=this_handles;

    end

    fprintf(1, 'delta_odor %d\n',mean(all_delta_odor));
    fprintf(1, 'delta_odor_on_reinf_on %d\n',mean(all_delta_odor_on_reinf_on));
    fprintf(1, 'all_delta_reinf %d\n',mean(all_delta_reinf));
end






