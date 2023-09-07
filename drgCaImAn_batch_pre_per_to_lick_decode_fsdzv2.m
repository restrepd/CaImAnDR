function drgCaImAn_batch_pre_per_to_lick_decode_fsdzv2(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a

    

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_lickdec_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_lick_decode_fsdzv2 run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

% new_no_files=handles.no_files;

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

if exist(handles.PathName_out)==0
    %create directory
    mkdir(handles.PathName_out)
end

training_time_from=handles.training_time_from;
training_delta_time=handles.training_delta_time;


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
    
    
    %Process each file separately
    for fileNo=first_file:handles.no_files
        tic
        first_toc=toc;
        
        handles_out=[];
        ii_out=0;
        
        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};
        
        [percent_correct] = drgCaImAnFindPercentCorrect(pre_per_PathName, pre_per_FileName);
        
        if percent_correct>=80
            
            drgCaImAn_add_trimmed_licks_to_pre_per(pre_per_FileName,pre_per_PathName);
            
            handles_choices.pre_per_PathName=pre_per_PathName;
            handles_choices.pre_per_FileName=pre_per_FileName;
            handles_choices.processing_algorithm=handles.processing_algorithm;
            handles_choices.MLalgo_to_use=handles.MLalgo_to_use;
            handles_choices.dt_p_threshold=handles.dt_p_threshold;
            handles_choices.show_figures=handles.show_figures;
            if handles.MLalgo_to_use==4
                handles_choices.no_repeats=handles.no_repeats;
            else
                handles_choices.no_repeats=1; %Note: the only algorithm that needs repeats is ANN, other aglorithms give the same result each run
            end
            handles_choices.k_fold=handles.k_fold;
            
            handles_choices.pre_time=handles.pre_time;
            handles_choices.ii_cost=handles.ii_cost;
            
            for ii_training=1:length(training_time_from)
                
                handles_choices.p_threshold=1.1;
                handles_choices.post_time=training_delta_time;
                handles_choices.post_shift=training_time_from(ii_training);
                
                ii_out=ii_out+1;
                handles_out.ii_out(ii_out).handles_choices=handles_choices;
                handles_out.ii_out(ii_out).grNo=handles.group(fileNo);
                handles_out.ii_out(ii_out).fileNo=fileNo;
                
                start_toc=toc;
                
%                 fprintf(1, ['Started processing file number %d, condition number %d\n'],fileNo,ii_training);
                
                handles_out.ii_out(ii_out).handles_out=drgCaImAn_repeat_decode_licks_entire_session(handles_choices);
                
                fprintf(1, ['Data processed for file number %d, condition number %d\n'],fileNo,ii_training);
                
%                 fprintf(1,'Processing time for drgCaImAn_pre_per_to_LDA_fsdz_new %d hours\n',(toc-start_toc)/(60*60));
            end
            
            fprintf(1,'\nProcessing time for file No %d is %d minutes\n\n',fileNo,(toc-first_toc)/60);
            
            %Save output file
            handles_out.last_file_processed=fileNo;
            handles_out.last_ii_out=ii_out;
            handles_out.handles=handles;
            save([handles.PathName_out pre_per_FileName(1:end-4) suffix_out],'handles_out','handles_choices','-v7.3')
        else
            fprintf(1,'\n\nFile %d not procressed because percent correct behavior <80\n',fileNo);
        end
    end
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end






