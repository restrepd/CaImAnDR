function drgCaImAn_batch_pre_per_to_LDA_fsdz(choiceBatchPathName,choiceFileName)
%Note: fitcnet will not work in Matlab versions earlier than 2021a


first_file=1;

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
        fprintf(1, ['Program will be terminated because file No %d, ' pre_per_PathName pre_per_FileName ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
end


if exist([handles.PathName_out handles.FileName_out])==0
    handles_out=[];
    ii_out=0;
    first_file=1;
else
    load([handles.PathName_out handles.FileName_out])
    ii_out=length(handles_out.ii_out);
    first_file=1+(length(handles_out.ii_out)/(length(handles.p_thr_less_than)*length(handles.MLalgo)));
end

figNo=0;
show_figures=1;

if all_files_present==1
    
    
    %Process each file separately
    for fileNo=first_file:length(handles.FileName_pre_per)
        tic
        first_toc=toc;
        
        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};
         
        figNo=0;
         
        for ii_p_thr=1:length(handles.p_thr_less_than)
            for MLalgo=handles.MLalgo
                
                ii_out=ii_out+1;

                handles_out.ii_out(ii_out).p_thr_less_than=handles.p_thr_less_than(ii_p_thr);
                handles_out.ii_out(ii_out).p_thr_more_than=handles.p_thr_more_than(ii_p_thr);
                handles_out.ii_out(ii_out).MLalgo=MLalgo;
                handles_out.ii_out(ii_out).grNo=handles.group(fileNo);
                handles_out.ii_out(ii_out).pre_per_PathName=pre_per_PathName;
                handles_out.ii_out(ii_out).pre_per_FileName=pre_per_FileName;
                start_toc=toc;
                switch processing_algo
                    case 1
                        handles_out.ii_out(ii_out).handles=drgCaImAn_pre_per_to_LDA_fsdz_new(pre_per_PathName, pre_per_FileName,handles.p_thr_less_than(ii_p_thr),handles.p_thr_more_than(ii_p_thr)...
                            ,MLalgo,show_figures,handles.no_sp_sm_trials_to_use, handles.first_sp_sm_trial_no,figNo,fileNo,handles.this_cost);
                    case 2
                        handles_out.ii_out(ii_out).handles=drgCaImAn_pre_per_to_LDA_fsdz_new_pre_trained_w_post(pre_per_PathName, pre_per_FileName,handles.p_thr_less_than(ii_p_thr),handles.p_thr_more_than(ii_p_thr)...
                            ,MLalgo,show_figures,handles.no_sp_sm_trials_to_use, handles.first_sp_sm_trial_no,figNo,fileNo,handles.this_cost);
                end
                figNo=figNo+3;
                fprintf(1, ['Data processed for file number %d, p_threshold from %d to %d and MLalgo %d\n'],fileNo,handles.p_thr_more_than(ii_p_thr),handles.p_thr_less_than(ii_p_thr),MLalgo);
                
                fprintf(1,'Processing time for drgCaImAn_pre_per_to_LDA_fsdz_new %d hours\n',(toc-start_toc)/(60*60));
            end
        end
        
        fprintf(1,'\n\nProcessing time for file No %d is %d hours\n',fileNo,(toc-first_toc)/(60*60));
        
        %Save output file
        handles_out.handles=handles;
        save([handles.PathName_out handles.FileName_out],'handles_out','-v7.3')
    end
    
   
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
     
end






