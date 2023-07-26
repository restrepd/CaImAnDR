function handles=drgCaImAn_dropc_choices_20220810_bla_mmHippoCamKGrin1 

%run with drgCalmAn_batch_dropc_fsdz

%%Grin4 fsds1
%This is the experiment Amyl Acetate 10% (S+) x Acetophenone 10% (S-)
%No lick or water reward given
%Mouse was water restricted, and the mouse was trained with Amyl Acetate as 
%S+ for many days, but no water on the test day 


handles.first_file=1;


%csv files for image analysis
for ii=1:4
    handles.PathNamecsv{ii}='C:\Users\Researcher\Documents\20220810-mmHippoCamKGrin1\';
end

handles.csvFileName{1}='Slide1.sld - 1-pwer9-Frame2-spm-081022_E_10_Iter_3500_output_ncorr_ext_bla.mat';
handles.csvFileName{2}='Slide1.sld - 2-pwer9-Frame2-spm-081022_E_10_Iter_3500_output_ncorr_ext_bla.mat';
handles.csvFileName{3}='Slide1.sld - 3-pwer9-Frame2-spm-081022_E_10_Iter_3500_output_ncorr_ext_bla.mat';
handles.csvFileName{4}='Slide1.sld - 4-pwer9-Frame2-spm-081022_E_10_Iter_3500_output_ncorr_ext_bla.mat';
 
handles.no_scans_per_image(1:4)=2;


%Other files
for ii=1:4
    handles.PathName{ii}='C:\Users\Researcher\Documents\20220810-mmHippoCamKGrin1\';
end


%spm files
handles.spmFileName{1}='1-mmHippoCamKGrin1-spm20220810T141703spm.mat';
handles.spmFileName{2}='2-mmHippoCamKGrin1-spm20220810T143836spm.mat';
handles.spmFileName{3}='3-mmHippoCamKGrin1-spm20220810T152630spm.mat';
handles.spmFileName{4}='4-mmHippoCamKGrin1-spm20220810T154800spm.mat';



%rhd files
handles.rhdFileName{1}='1-CA1-mmHippoCamKGrin1-spm_220810_141644.rhd';
handles.rhdFileName{2}='2-CA1-mmHippoCamKGrin1-spm_220810_143808.rhd';
handles.rhdFileName{3}='3-CA1-mmHippoCamKGrin1-spm_220810_152614.rhd';
handles.rhdFileName{4}='4-CA1-mmHippoCamKGrin1-spm_220810_154740.rhd';



handles.no_files=length(handles.csvFileName);
