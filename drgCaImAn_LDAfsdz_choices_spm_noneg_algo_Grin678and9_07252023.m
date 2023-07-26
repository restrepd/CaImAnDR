function handles=drgCaImAn_LDAfsdz_choices_spm_noneg_algo_Grin678and9_07252023

%run with drgCalmAn_batch_dropc_fsdz
% handles.PathName_out='F:\SFTP\Ming Ma\Grin6and7_spm_raw\';
% handles.FileName_out='entire_session_deepcad6_02052022.mat';

%Path and suffix for file output
for fileNo=1:68
    handles.PathName_out{fileNo}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\';
end
for fileNo=69:131
    handles.PathName_out{fileNo}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\';
end

handles.suffix_out='_alg1.mat';
handles.MLalgo_to_use=[6]; %Compare algorithms change to 1 to 6
%1=LDA
%2=SVM
%3=Bayesian naive, I get some errors, maybe you can fix them?
%4=Artificial neural network
%5=Decision tree
%6=glm


handles.processing_algorithm=3; %Very impotrant, 1 and 2 have not been vetted
handles.post_time=[5];
handles.k_fold=[5]; %Not used for processing_algorithm=3
handles.post_shift=[0.5];
handles.pre_time=[5];

handles.p_threshold=[1.1]; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
handles.dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
handles.show_figures=0; %Show the figures
handles.ii_cost=3; %Cost for decoding models

handles.first_file=1;
handles.skip_files=[];

handles.group_names{1}='forward';
handles.group_names{2}='reverse';
% handles.group_names{3}='\\homeodor_female(kx_anesthesia)';
% handles.group_names{4}='\\homeodor_male';
% handles.group_names{5}='\\homeodor_sameodor';
% handles.group_names{6}='odor_off';

% handles.group_names{6}='\\homeodor_male_pcdh21'; %Note: I am doing pcdh1
% and male together
% handles.session_no_per_mouse(1)=1; means day1 of training
%Sex
% 0 Male
% 1 Female



handles.mouse{1}='GRIN8';
handles.odor_pair{1}='HeptanalMO';
handles.session_no_per_mouse(1)=1;
handles.group(1)=1;
handles.sex(1)=0;
handles.PathName_pre_per{1}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220323-mmHippoFGrin8\';
handles.FileName_pre_per{1}='Slide2.sld - 1-pwer20-Frame2-spm-032322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{2}='GRIN8';
handles.odor_pair{2}='HeptanalMO';
handles.session_no_per_mouse(2)=1;
handles.group(2)=1;
handles.sex(2)=0;
handles.PathName_pre_per{2}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220323-mmHippoFGrin8\';
handles.FileName_pre_per{2}='Slide2.sld - 2-pwer20-Frame2-spm-032322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{3}='GRIN8';
handles.odor_pair{3}='HeptanalMO';
handles.session_no_per_mouse(3)=1;
handles.group(3)=1;
handles.sex(3)=0;
handles.PathName_pre_per{3}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220323-mmHippoFGrin8\';
handles.FileName_pre_per{3}='Slide2.sld - 3-pwer20-Frame2-spm-032322 - 2_E_10_Iter_3500_output-1_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{4}='GRIN8';
handles.odor_pair{4}='HeptanalMO';
handles.session_no_per_mouse(4)=2;
handles.group(4)=1;
handles.sex(4)=0;
handles.PathName_pre_per{4}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220325-mmHippoFGrin8\';
handles.FileName_pre_per{4}='Slide1.sld - 1-pwer19-Frame2-spm-032522_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


% not working
% Error using classreg.learning.classif.FullClassificationModel.processCost
% Misclassification cost matrix must be 1-by-1.
handles.mouse{5}='GRIN8';
handles.odor_pair{5}='HeptanalMO';
handles.session_no_per_mouse(5)=2;
handles.group(5)=1;
handles.sex(5)=0;
handles.PathName_pre_per{5}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220325-mmHippoFGrin8\';
handles.FileName_pre_per{5}='Slide1.sld - 2-pwer19-Frame2-spm-032522_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{6}='GRIN8';
handles.odor_pair{6}='HeptanalMO';
handles.session_no_per_mouse(6)=3;
handles.group(6)=1;
handles.sex(6)=0;
handles.PathName_pre_per{6}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220328-mmHippoFGrin8\';
handles.FileName_pre_per{6}='Slide1.sld - 1-pwer16-Frame2-spm-032822_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


%%not working
% Error using classreg.learning.classif.FullClassificationModel.processCost
% Misclassification cost matrix must be 1-by-1.

handles.mouse{7}='GRIN8';
handles.odor_pair{7}='HeptanalMO';
handles.session_no_per_mouse(7)=3;
handles.group(7)=1;
handles.sex(7)=0;
handles.PathName_pre_per{7}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220328-mmHippoFGrin8\';
handles.FileName_pre_per{7}='Slide1.sld - 2-pwer16-Frame2-spm-032822 - 2_E_10_Iter_3500_output-1_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{8}='GRIN8';
handles.odor_pair{8}='HeptanalMO';
handles.session_no_per_mouse(8)=3;
handles.group(8)=1;
handles.sex(8)=0;
handles.PathName_pre_per{8}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220328-mmHippoFGrin8\';
handles.FileName_pre_per{8}='Slide1.sld - 3-pwer16-Frame2-spm-032822 - 1_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{9}='GRIN8';
handles.odor_pair{9}='HeptanalMO';
handles.session_no_per_mouse(9)=4;
handles.group(9)=1;
handles.sex(9)=0;
handles.PathName_pre_per{9}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220330-mmHippoFGrin8\';
handles.FileName_pre_per{9}='Slide1.sld - 1-pwer19-Frame2-spm-033022_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{10}='GRIN8';
handles.odor_pair{10}='HeptanalMO';
handles.session_no_per_mouse(10)=4;
handles.group(10)=1;
handles.sex(10)=0;
handles.PathName_pre_per{10}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220330-mmHippoFGrin8\';
handles.FileName_pre_per{10}='Slide1.sld - 2-pwer19-Frame2-spm-033022_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{11}='GRIN8';
handles.odor_pair{11}='HeptanalMO';
handles.session_no_per_mouse(11)=4;
handles.group(11)=1;
handles.sex(11)=0;
handles.PathName_pre_per{11}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220330-mmHippoFGrin8\';
handles.FileName_pre_per{11}='Slide1.sld - 4-pwer19-Frame2-spm-033022_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{12}='GRIN8';
handles.odor_pair{12}='HeptanalMO';
handles.session_no_per_mouse(12)=4;
handles.group(12)=1;
handles.sex(12)=0;
handles.PathName_pre_per{12}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220330-mmHippoFGrin8\';
handles.FileName_pre_per{12}='Slide1.sld - 5-pwer19-Frame2-spm-033022_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{13}='GRIN8';
handles.odor_pair{13}='HeptanalMO';
handles.session_no_per_mouse(13)=4;
handles.group(13)=1;
handles.sex(13)=0;
handles.PathName_pre_per{13}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220330-mmHippoFGrin8\';
handles.FileName_pre_per{13}='Slide1.sld - 6-pwer19-Frame2-spm-033022_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{14}='GRIN8';
handles.odor_pair{14}='HeptanalMO';
handles.session_no_per_mouse(14)=5;
handles.group(14)=1;
handles.sex(14)=0;
handles.PathName_pre_per{14}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220401-mmHippoFGrin8\';
handles.FileName_pre_per{14}='Slide1.sld - 1-pwer20-Frame2-spm-040122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{15}='GRIN8';
handles.odor_pair{15}='HeptanalMO';
handles.session_no_per_mouse(15)=5;
handles.group(15)=1;
handles.sex(15)=0;
handles.PathName_pre_per{15}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220401-mmHippoFGrin8\';
handles.FileName_pre_per{15}='Slide1.sld - 2-pwer20-Frame2-spm-040122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{16}='GRIN8';
handles.odor_pair{16}='HeptanalMO';
handles.session_no_per_mouse(16)=5;
handles.group(16)=1;
handles.sex(16)=0;
handles.PathName_pre_per{16}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220401-mmHippoFGrin8\';
handles.FileName_pre_per{16}='Slide1.sld - 3-pwer20-Frame2-spm-040122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{17}='GRIN8';
handles.odor_pair{17}='HeptanalMO';
handles.session_no_per_mouse(17)=5;
handles.group(17)=1;
handles.sex(17)=0;
handles.PathName_pre_per{17}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220401-mmHippoFGrin8\';
handles.FileName_pre_per{17}='Slide1.sld - 4-pwer20-Frame2-spm-040122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{18}='GRIN8';
handles.odor_pair{18}='HeptanalMO';
handles.session_no_per_mouse(18)=6;
handles.group(18)=1;
handles.sex(18)=0;
handles.PathName_pre_per{18}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220404-mmHippoFGrin8\';
handles.FileName_pre_per{18}='Slide1.sld - 1-pwer20-Frame2-spm-040422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{19}='GRIN8';
handles.odor_pair{19}='MOHeptanal';
handles.session_no_per_mouse(19)=6;
handles.group(19)=2;
handles.sex(19)=0;
handles.PathName_pre_per{19}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220404-mmHippoFGrin8\';
handles.FileName_pre_per{19}='Slide1.sld - 2-pwer20-Frame2-reverspm-040422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{20}='GRIN8';
handles.odor_pair{20}='MOHeptanal';
handles.session_no_per_mouse(20)=6;
handles.group(20)=2;
handles.sex(20)=0;
handles.PathName_pre_per{20}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220404-mmHippoFGrin8\';
handles.FileName_pre_per{20}='Slide1.sld - 3-pwer20-Frame2-reverspm-040422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


% This resulted in an error in drgCaImAn_SVZ_entire_session
% Error using ClassificationNaiveBayes\findNoDataCombos (line 347)
% A normal distribution cannot be fit for the combination of class 1 and
% predictor these_training_measurements. The data has zero variance.
handles.mouse{21}='GRIN8';
handles.odor_pair{21}='MOHeptanal';
handles.session_no_per_mouse(21)=6;
handles.group(21)=2;
handles.sex(21)=0;
handles.PathName_pre_per{21}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220404-mmHippoFGrin8\';
handles.FileName_pre_per{21}='Slide1.sld - 4-pwer20-Frame2-reverspm-040422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{22}='GRIN8';
handles.odor_pair{22}='MOHeptanal';
handles.session_no_per_mouse(22)=7;
handles.group(22)=2;
handles.sex(22)=0;
handles.PathName_pre_per{22}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220406-mmHippoFGrin8\';
handles.FileName_pre_per{22}='Slide1.sld - 1-pwer20-Frame2-reverspm-040622_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{23}='GRIN8';
handles.odor_pair{23}='MOHeptanal';
handles.session_no_per_mouse(23)=7;
handles.group(23)=2;
handles.sex(23)=0;
handles.PathName_pre_per{23}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220406-mmHippoFGrin8\';
handles.FileName_pre_per{23}='Slide1.sld - 2-pwer20-Frame2-reverspm-040622_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{24}='GRIN8';
handles.odor_pair{24}='MOHeptanal';
handles.session_no_per_mouse(24)=7;
handles.group(24)=2;
handles.sex(24)=0;
handles.PathName_pre_per{24}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220406-mmHippoFGrin8\';
handles.FileName_pre_per{24}='Slide1.sld - 3-pwer20-Frame2-reverspm-040622_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{25}='GRIN8';
handles.odor_pair{25}='MOHeptanal';
handles.session_no_per_mouse(25)=7;
handles.group(25)=2;
handles.sex(25)=0;
handles.PathName_pre_per{25}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220406-mmHippoFGrin8\';
handles.FileName_pre_per{25}='Slide1.sld - 4-pwer20-Frame2-reverspm-040622_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{26}='GRIN8';
handles.odor_pair{26}='MOHeptanal';
handles.session_no_per_mouse(26)=7;
handles.group(26)=2;
handles.sex(26)=0;
handles.PathName_pre_per{26}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220406-mmHippoFGrin8\';
handles.FileName_pre_per{26}='Slide1.sld - 5-pwer20-Frame2-reverspm-040622_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{27}='GRIN8';
handles.odor_pair{27}='MOHeptanal';
handles.session_no_per_mouse(27)=8;
handles.group(27)=2;
handles.sex(27)=0;
handles.PathName_pre_per{27}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220408-mmHippoFGrin8\';
handles.FileName_pre_per{27}='Slide1.sld - 1-pwer20-Frame2-reverspm-040822 - 1_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{28}='GRIN8';
handles.odor_pair{28}='MOHeptanal';
handles.session_no_per_mouse(28)=8;
handles.group(28)=2;
handles.sex(28)=0;
handles.PathName_pre_per{28}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220408-mmHippoFGrin8\';
handles.FileName_pre_per{28}='Slide1.sld - 2-pwer20-Frame2-reverspm-040822_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{29}='GRIN8';
handles.odor_pair{29}='MOHeptanal';
handles.session_no_per_mouse(29)=8;
handles.group(29)=2;
handles.sex(29)=0;
handles.PathName_pre_per{29}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220408-mmHippoFGrin8\';
handles.FileName_pre_per{29}='Slide1.sld - 3-pwer20-Frame2-reverspm-040822_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{30}='GRIN8';
handles.odor_pair{30}='MOHeptanal';
handles.session_no_per_mouse(30)=8;
handles.group(30)=2;
handles.sex(30)=0;
handles.PathName_pre_per{30}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220408-mmHippoFGrin8\';
handles.FileName_pre_per{30}='Slide1.sld - 4-pwer20-Frame2-reverspm-040822_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{31}='GRIN8';
handles.odor_pair{31}='MOHeptanal';
handles.session_no_per_mouse(31)=8;
handles.group(31)=2;
handles.sex(31)=0;
handles.PathName_pre_per{31}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220408-mmHippoFGrin8\';
handles.FileName_pre_per{31}='Slide1.sld - 5-pwer20-Frame2-reverspm-040822_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{32}='GRIN8';
handles.odor_pair{32}='MOHeptanal';
handles.session_no_per_mouse(32)=9;
handles.group(32)=2;
handles.sex(32)=0;
handles.PathName_pre_per{32}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220411-mmHippoFGrin8\';
handles.FileName_pre_per{32}='Slide1.sld - 1-pwer21-Frame2-reverspm-041122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{33}='GRIN8';
handles.odor_pair{33}='MOHeptanal';
handles.session_no_per_mouse(33)=9;
handles.group(33)=2;
handles.sex(33)=0;
handles.PathName_pre_per{33}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220411-mmHippoFGrin8\';
handles.FileName_pre_per{33}='Slide1.sld - 2-pwer21-Frame2-reverspm-041122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{34}='GRIN8';
handles.odor_pair{34}='MOHeptanal';
handles.session_no_per_mouse(34)=9;
handles.group(34)=2;
handles.sex(34)=0;
handles.PathName_pre_per{34}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220411-mmHippoFGrin8\';
handles.FileName_pre_per{34}='Slide1.sld - 3-pwer21-Frame2-reverspm-041122 - 1_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{35}='GRIN8';
handles.odor_pair{35}='MOHeptanal';
handles.session_no_per_mouse(35)=9;
handles.group(35)=2;
handles.sex(35)=0;
handles.PathName_pre_per{35}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220411-mmHippoFGrin8\';
handles.FileName_pre_per{35}='Slide1.sld - 4-pwer21-Frame2-reverspm-041122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{36}='GRIN8';
handles.odor_pair{36}='MOHeptanal';
handles.session_no_per_mouse(36)=9;
handles.group(36)=2;
handles.sex(36)=0;
handles.PathName_pre_per{36}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220411-mmHippoFGrin8\';
handles.FileName_pre_per{36}='Slide1.sld - 5-pwer21-Frame2-reverspm-041122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{37}='GRIN8';
handles.odor_pair{37}='MOHeptanal';
handles.session_no_per_mouse(37)=10;
handles.group(37)=2;
handles.sex(37)=0;
handles.PathName_pre_per{37}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220413-mmHippoFGrin8\';
handles.FileName_pre_per{37}='Slide1.sld - 1-pwer21-Frame2-reverspm-041322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{38}='GRIN8';
handles.odor_pair{38}='MOHeptanal';
handles.session_no_per_mouse(38)=10;
handles.group(38)=2;
handles.sex(38)=0;
handles.PathName_pre_per{38}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220413-mmHippoFGrin8\';
handles.FileName_pre_per{38}='Slide1.sld - 2-pwer21-Frame2-reverspm-041322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{39}='GRIN8';
handles.odor_pair{39}='MOHeptanal';
handles.session_no_per_mouse(39)=10;
handles.group(39)=2;
handles.sex(39)=0;
handles.PathName_pre_per{39}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220413-mmHippoFGrin8\';
handles.FileName_pre_per{39}='Slide1.sld - 3-pwer21-Frame2-reverspm-041322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{40}='GRIN8';
handles.odor_pair{40}='MOHeptanal';
handles.session_no_per_mouse(40)=10;
handles.group(40)=2;
handles.sex(40)=0;
handles.PathName_pre_per{40}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220413-mmHippoFGrin8\';
handles.FileName_pre_per{40}='Slide1.sld - 4-pwer21-Frame2-reverspm-041322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{41}='GRIN8';
handles.odor_pair{41}='MOHeptanal';
handles.session_no_per_mouse(41)=10;
handles.group(41)=2;
handles.sex(41)=0;
handles.PathName_pre_per{41}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220413-mmHippoFGrin8\';
handles.FileName_pre_per{41}='Slide1.sld - 5-pwer21-Frame2-reverspm-041322_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';






handles.mouse{42}='GRIN9';
handles.odor_pair{42}='HeptanalMO';
handles.session_no_per_mouse(42)=1;
handles.group(42)=1;
handles.sex(42)=1;
handles.PathName_pre_per{42}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{42}='Slide1.sld - 1-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{43}='GRIN9';
handles.odor_pair{43}='HeptanalMO';
handles.session_no_per_mouse(43)=1;
handles.group(43)=1;
handles.sex(43)=1;
handles.PathName_pre_per{43}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{43}='Slide1.sld - 2-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{44}='GRIN9';
handles.odor_pair{44}='HeptanalMO';
handles.session_no_per_mouse(44)=1;
handles.group(44)=1;
handles.sex(44)=1;
handles.PathName_pre_per{44}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{44}='Slide1.sld - 3-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{45}='GRIN9';
handles.odor_pair{45}='HeptanalMO';
handles.session_no_per_mouse(45)=1;
handles.group(45)=1;
handles.sex(45)=1;
handles.PathName_pre_per{45}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{45}='Slide1.sld - 4-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{46}='GRIN9';
handles.odor_pair{46}='HeptanalMO';
handles.session_no_per_mouse(46)=1;
handles.group(46)=1;
handles.sex(46)=1;
handles.PathName_pre_per{46}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{46}='Slide1.sld - 5-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{47}='GRIN9';
handles.odor_pair{47}='HeptanalMO';
handles.session_no_per_mouse(47)=1;
handles.group(47)=1;
handles.sex(47)=1;
handles.PathName_pre_per{47}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{47}='Slide1.sld - 6-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{48}='GRIN9';
handles.odor_pair{48}='HeptanalMO';
handles.session_no_per_mouse(48)=1;
handles.group(48)=1;
handles.sex(48)=1;
handles.PathName_pre_per{48}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220427-mmHippoFGrin9\';
handles.FileName_pre_per{48}='Slide1.sld - 7-pwer21-Frame2-spm-042722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{49}='GRIN9';
handles.odor_pair{49}='HeptanalMO';
handles.session_no_per_mouse(49)=2;
handles.group(49)=1;
handles.sex(49)=1;
handles.PathName_pre_per{49}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220429-mmHippoFGrin9\';
handles.FileName_pre_per{49}='Slide1.sld - 1-pwer21-Frame2-spm-042922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{50}='GRIN9';
handles.odor_pair{50}='HeptanalMO';
handles.session_no_per_mouse(50)=2;
handles.group(50)=1;
handles.sex(50)=1;
handles.PathName_pre_per{50}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220429-mmHippoFGrin9\';
handles.FileName_pre_per{50}='Slide1.sld - 2-pwer21-Frame2-spm-042922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{51}='GRIN9';
handles.odor_pair{51}='MOHeptanal';
handles.session_no_per_mouse(51)=2;
handles.group(51)=2;
handles.sex(51)=1;
handles.PathName_pre_per{51}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220429-mmHippoFGrin9\';
handles.FileName_pre_per{51}='Slide1.sld - 3-pwer21-Frame2-reversespm-042922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{52}='GRIN9';
handles.odor_pair{52}='MOHeptanal';
handles.session_no_per_mouse(52)=2;
handles.group(52)=2;
handles.sex(52)=1;
handles.PathName_pre_per{52}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220429-mmHippoFGrin9\';
handles.FileName_pre_per{52}='Slide1.sld - 4-pwer21-Frame2-reversespm-042922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{53}='GRIN9';
handles.odor_pair{53}='MOHeptanal';
handles.session_no_per_mouse(53)=2;
handles.group(53)=2;
handles.sex(53)=1;
handles.PathName_pre_per{53}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220429-mmHippoFGrin9\';
handles.FileName_pre_per{53}='Slide1.sld - 5-pwer21-Frame2-reversespm-042922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{54}='GRIN9';
handles.odor_pair{54}='MOHeptanal';
handles.session_no_per_mouse(54)=2;
handles.group(54)=2;
handles.sex(54)=1;
handles.PathName_pre_per{54}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220429-mmHippoFGrin9\';
handles.FileName_pre_per{54}='Slide1.sld - 6-pwer21-Frame2-reversespm-042922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


% This resulted in an error in drgCaImAn_SVZ_entire_session
% Error using ClassificationNaiveBayes\findNoDataCombos (line 347)
% A normal distribution cannot be fit for the combination of class 1 and
% predictor these_training_measurements. The data has zero variance.
handles.mouse{55}='GRIN9';
handles.odor_pair{55}='MOHeptanal';
handles.session_no_per_mouse(55)=3;
handles.group(55)=2;
handles.sex(55)=1;
handles.PathName_pre_per{55}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220502-mmHippoFGrin9\';
handles.FileName_pre_per{55}='Slide1.sld - 1-pwer21-Frame2-reversespm-050222_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{56}='GRIN9';
handles.odor_pair{56}='MOHeptanal';
handles.session_no_per_mouse(56)=3;
handles.group(56)=2;
handles.sex(56)=1;
handles.PathName_pre_per{56}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220502-mmHippoFGrin9\';
handles.FileName_pre_per{56}='Slide1.sld - 2-pwer21-Frame2-reversespm-050222_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{57}='GRIN9';
handles.odor_pair{57}='MOHeptanal';
handles.session_no_per_mouse(57)=3;
handles.group(57)=2;
handles.sex(57)=1;
handles.PathName_pre_per{57}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220502-mmHippoFGrin9\';
handles.FileName_pre_per{57}='Slide1.sld - 3-pwer21-Frame2-reversespm-050222_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{58}='GRIN9';
handles.odor_pair{58}='MOHeptanal';
handles.session_no_per_mouse(58)=3;
handles.group(58)=2;
handles.sex(58)=1;
handles.PathName_pre_per{58}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220502-mmHippoFGrin9\';
handles.FileName_pre_per{58}='Slide1.sld - 4-pwer21-Frame2-reversespm-050222_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{59}='GRIN9';
handles.odor_pair{59}='MOHeptanal';
handles.session_no_per_mouse(59)=3;
handles.group(59)=2;
handles.sex(59)=1;
handles.PathName_pre_per{59}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220502-mmHippoFGrin9\';
handles.FileName_pre_per{59}='Slide1.sld - 5-pwer21-Frame2-reversespm-050222_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{60}='GRIN9';
handles.odor_pair{60}='MOHeptanal';
handles.session_no_per_mouse(60)=4;
handles.group(60)=2;
handles.sex(60)=1;
handles.PathName_pre_per{60}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220504-mmHippoFGrin9\';
handles.FileName_pre_per{60}='Slide1.sld - 1-pwer21-Frame2-reverspm-050422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{61}='GRIN9';
handles.odor_pair{61}='MOHeptanal';
handles.session_no_per_mouse(61)=4;
handles.group(61)=2;
handles.sex(61)=1;
handles.PathName_pre_per{61}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220504-mmHippoFGrin9\';
handles.FileName_pre_per{61}='Slide1.sld - 2-pwer21-Frame2-reverspm-050422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{62}='GRIN9';
handles.odor_pair{62}='MOHeptanal';
handles.session_no_per_mouse(62)=4;
handles.group(62)=2;
handles.sex(62)=1;
handles.PathName_pre_per{62}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220504-mmHippoFGrin9\';
handles.FileName_pre_per{62}='Slide1.sld - 3-pwer21-Frame2-reverspm-050422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{63}='GRIN9';
handles.odor_pair{63}='MOHeptanal';
handles.session_no_per_mouse(63)=4;
handles.group(63)=2;
handles.sex(63)=1;
handles.PathName_pre_per{63}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220504-mmHippoFGrin9\';
handles.FileName_pre_per{63}='Slide1.sld - 4-pwer21-Frame2-reverspm-050422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{64}='GRIN9';
handles.odor_pair{64}='MOHeptanal';
handles.session_no_per_mouse(64)=4;
handles.group(64)=2;
handles.sex(64)=1;
handles.PathName_pre_per{64}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220504-mmHippoFGrin9\';
handles.FileName_pre_per{64}='Slide1.sld - 5-pwer21-Frame2-reverspm-050422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{65}='GRIN9';
handles.odor_pair{65}='MOHeptanal';
handles.session_no_per_mouse(65)=4;
handles.group(65)=2;
handles.sex(65)=1;
handles.PathName_pre_per{65}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220504-mmHippoFGrin9\';
handles.FileName_pre_per{65}='Slide1.sld - 6-pwer21-Frame2-reverspm-050422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{66}='GRIN9';
handles.odor_pair{66}='MOHeptanal';
handles.session_no_per_mouse(66)=5;
handles.group(66)=2;
handles.sex(66)=1;
handles.PathName_pre_per{66}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220505-mmHippoFGrin9\';
handles.FileName_pre_per{66}='Slide1.sld - 1-pwer21-frame2-reverspm-050522_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{67}='GRIN9';
handles.odor_pair{67}='MOHeptanal';
handles.session_no_per_mouse(67)=5;
handles.group(67)=2;
handles.sex(67)=1;
handles.PathName_pre_per{67}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220505-mmHippoFGrin9\';
handles.FileName_pre_per{67}='Slide1.sld - 2-pwer21-frame2-reverspm-050522_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{68}='GRIN9';
handles.odor_pair{68}='MOHeptanal';
handles.session_no_per_mouse(68)=5;
handles.group(68)=2;
handles.sex(68)=1;
handles.PathName_pre_per{68}='F:\SFTP\Ming Ma\Grin8and9_spm_nonneg\20220505-mmHippoFGrin9\';
handles.FileName_pre_per{68}='Slide1.sld - 3-pwer21-frame2-reverspm-050522_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

%6 abd 7
%Note: The _dec2 files are stored in the F:\SFTP\Ming
%Ma\Grin6and7_spm_nonneg\ directory

handles.mouse{69}='GRIN6';
handles.odor_pair{69}='HeptanalMO';
handles.session_no_per_mouse(69)=1;
handles.group(69)=1;
handles.sex(69)=1;
handles.PathName_pre_per{69}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220124-mmHippoFGrin6\';
handles.FileName_pre_per{69}='1-20pwer-Frame2-spm_XY0_Z0_T0000_C0_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{70}='GRIN6';
handles.odor_pair{70}='HeptanalMO';
handles.session_no_per_mouse(70)=1;
handles.group(70)=1;
handles.sex(70)=1;
handles.PathName_pre_per{70}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220124-mmHippoFGrin6\';
handles.FileName_pre_per{70}='2-20pwer-Frame3-spm_XY0_Z0_T00000_C0_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{71}='GRIN6';
handles.odor_pair{71}='HeptanalMO';
handles.session_no_per_mouse(71)=1;
handles.group(71)=1;
handles.sex(71)=1;
handles.PathName_pre_per{71}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220124-mmHippoFGrin6\';
handles.FileName_pre_per{71}='3-20pwer-Frame3-spm_XY0_Z0_T00000_C0_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{72}='GRIN6';
handles.odor_pair{72}='HeptanalMO';
handles.session_no_per_mouse(72)=2;
handles.group(72)=1;
handles.sex(72)=1;
handles.PathName_pre_per{72}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin6\';
handles.FileName_pre_per{72}='Slide1.sld - 1-pwer18-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{73}='GRIN6';
handles.odor_pair{73}='HeptanalMO';
handles.session_no_per_mouse(73)=2;
handles.group(73)=1;
handles.sex(73)=1;
handles.PathName_pre_per{73}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin6\';
handles.FileName_pre_per{73}='Slide1.sld - 2-pwer18-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{74}='GRIN6';
handles.odor_pair{74}='HeptanalMO';
handles.session_no_per_mouse(74)=2;
handles.group(74)=1;
handles.sex(74)=1;
handles.PathName_pre_per{74}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin6\';
handles.FileName_pre_per{74}='Slide1.sld - 3-pwer18-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{75}='GRIN6';
handles.odor_pair{75}='HeptanalMO';
handles.session_no_per_mouse(7)=2;
handles.group(75)=1;
handles.sex(75)=1;
handles.PathName_pre_per{75}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin6\';
handles.FileName_pre_per{75}='Slide1.sld - 4-pwer18-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{76}='GRIN6';
handles.odor_pair{76}='HeptanalMO';
handles.session_no_per_mouse(76)=3;
handles.group(76)=1;
handles.sex(76)=1;
% handles.PathName_pre_per{76}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin6\';
handles.PathName_pre_per{76}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin6_rec\';
% handles.FileName_pre_per{76}='Slide1.sld - 1-18pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';
handles.FileName_pre_per{76}='Slide1.sld -20220128-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{77}='GRIN6';
handles.odor_pair{77}='HeptanalMO';
handles.session_no_per_mouse(77)=3;
handles.group(77)=1;
handles.sex(77)=1;
handles.PathName_pre_per{77}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin6\';
handles.FileName_pre_per{77}='Slide1.sld - 2-18pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{78}='GRIN6';
handles.odor_pair{78}='HeptanalMO';
handles.session_no_per_mouse(78)=3;
handles.group(78)=1;
handles.sex(78)=1;
handles.PathName_pre_per{78}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin6\';
handles.FileName_pre_per{78}='Slide1.sld - 3-18pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{79}='GRIN6';
handles.odor_pair{79}='HeptanalMO';
handles.session_no_per_mouse(79)=3;
handles.group(79)=1;
handles.sex(79)=1;
handles.PathName_pre_per{79}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin6\';
handles.FileName_pre_per{79}='Slide1.sld - 4-18pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{80}='GRIN6';
handles.odor_pair{80}='HeptanalMO';
handles.session_no_per_mouse(80)=4;
handles.group(80)=1;
handles.sex(80)=1;
% handles.PathName_pre_per{80}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin6\';
% handles.FileName_pre_per{80}='Slide1.sld - 1-18pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';
handles.PathName_pre_per{80}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin6_rec\';
handles.FileName_pre_per{80}='Slide1.sld -20220131-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{81}='GRIN6';
handles.odor_pair{81}='MOHeptanal';
handles.session_no_per_mouse(81)=4;
handles.group(81)=2;
handles.sex(81)=1;
handles.PathName_pre_per{81}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin6\';
handles.FileName_pre_per{81}='Slide1.sld - 2-18pwer-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{82}='GRIN6';
handles.odor_pair{82}='MOHeptanal';
handles.session_no_per_mouse(82)=4;
handles.group(82)=2;
handles.sex(82)=1;
handles.PathName_pre_per{82}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin6\';
handles.FileName_pre_per{82}='Slide1.sld - 3-18pwer-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{83}='GRIN6';
handles.odor_pair{83}='MOHeptanal';
handles.session_no_per_mouse(83)=4;
handles.group(83)=2;
handles.sex(83)=1;
handles.PathName_pre_per{83}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin6\';
handles.FileName_pre_per{83}='Slide1.sld - 4-18pwer-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{84}='GRIN6';
handles.odor_pair{84}='MOHeptanal';
handles.session_no_per_mouse(84)=5;
handles.group(84)=2;
handles.sex(84)=1;
handles.PathName_pre_per{84}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin6\';
handles.FileName_pre_per{84}='Slide1.sld - 1-pwer18-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{85}='GRIN6';
handles.odor_pair{85}='MOHeptanal';
handles.session_no_per_mouse(85)=5;
handles.group(85)=2;
handles.sex(85)=1;
handles.PathName_pre_per{85}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin6\';
handles.FileName_pre_per{85}='Slide1.sld - 2-pwer18-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{86}='GRIN6';
handles.odor_pair{86}='MOHeptanal';
handles.session_no_per_mouse(86)=5;
handles.group(86)=2;
handles.sex(86)=1;
handles.PathName_pre_per{86}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin6\';
handles.FileName_pre_per{86}='Slide1.sld - 3-pwer18-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{87}='GRIN6';
handles.odor_pair{87}='MOHeptanal';
handles.session_no_per_mouse(87)=5;
handles.group(87)=2;
handles.sex(87)=1;
handles.PathName_pre_per{87}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin6\';
handles.FileName_pre_per{87}='Slide1.sld - 4-pwer18-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{88}='GRIN6';
handles.odor_pair{88}='MOHeptanal';
handles.session_no_per_mouse(88)=5;
handles.group(88)=2;
handles.sex(88)=1;
handles.PathName_pre_per{88}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin6\';
handles.FileName_pre_per{88}='Slide1.sld - 5-pwer18-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


% This resulted in an error in drgCaImAn_SVZ_entire_session
% Error using ClassificationNaiveBayes\findNoDataCombos (line 347)
% A normal distribution cannot be fit for the combination of class 1 and
% predictor these_training_measurements. The data has zero variance.
% handles.mouse{21}='GRIN6';
% handles.odor_pair{21}='MOHeptanal';
% handles.session_no_per_mouse(21)=6;
% handles.group(21)=2;
% handles.sex(21)=1;
% handles.PathName_pre_per{21}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin6\';
% handles.FileName_pre_per{21}='Slide1.sld - 1-pwer18-Frame2-revesespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{89}='GRIN6';
handles.odor_pair{89}='MOHeptanal';
handles.session_no_per_mouse(89)=6;
handles.group(89)=2;
handles.sex(89)=1;
handles.PathName_pre_per{89}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin6\';
handles.FileName_pre_per{89}='Slide1.sld - 2-pwer18-Frame2-revesespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{90}='GRIN6';
handles.odor_pair{90}='MOHeptanal';
handles.session_no_per_mouse(90)=6;
handles.group(90)=2;
handles.sex(90)=1;
handles.PathName_pre_per{90}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin6\';
handles.FileName_pre_per{90}='Slide1.sld - 3-pwer18-Frame2-revesespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{91}='GRIN6';
handles.odor_pair{91}='MOHeptanal';
handles.session_no_per_mouse(91)=6;
handles.group(91)=2;
handles.sex(91)=1;
handles.PathName_pre_per{91}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin6\';
handles.FileName_pre_per{91}='Slide1.sld - 4-pwer18-Frame2-revesespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{92}='GRIN6';
handles.odor_pair{92}='MOHeptanal';
handles.session_no_per_mouse(92)=6;
handles.group(92)=2;
handles.sex(92)=1;
handles.PathName_pre_per{92}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin6\';
handles.FileName_pre_per{92}='Slide1.sld - 5-pwer18-Frame2-revesespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{93}='GRIN6';
handles.odor_pair{93}='MOHeptanal';
handles.session_no_per_mouse(93)=7;
handles.group(93)=2;
handles.sex(93)=1;
handles.PathName_pre_per{93}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220207-mmHippoFGrin6\';
handles.FileName_pre_per{93}='Slide1.sld - 1-pwer18-Frame2-reversespm-020722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{94}='GRIN6';
handles.odor_pair{94}='MOHeptanal';
handles.session_no_per_mouse(94)=7;
handles.group(94)=2;
handles.sex(94)=1;
handles.PathName_pre_per{94}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220207-mmHippoFGrin6\';
handles.FileName_pre_per{94}='Slide1.sld - 2-pwer18-Frame2-reversespm-020722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{95}='GRIN6';
handles.odor_pair{95}='MOHeptanal';
handles.session_no_per_mouse(95)=7;
handles.group(95)=2;
handles.sex(95)=1;
handles.PathName_pre_per{95}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220207-mmHippoFGrin6\';
handles.FileName_pre_per{95}='Slide1.sld - 3-pwer18-Frame2-reversespm-020722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{96}='GRIN6';
handles.odor_pair{96}='MOHeptanal';
handles.session_no_per_mouse(96)=8;
handles.group(96)=2;
handles.sex(96)=1;
handles.PathName_pre_per{96}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin6\';
handles.FileName_pre_per{96}='Slide1.sld - 1-pwer19-Frame2-reversespm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{97}='GRIN6';
handles.odor_pair{97}='MOHeptanal';
handles.session_no_per_mouse(97)=8;
handles.group(97)=2;
handles.sex(97)=1;
handles.PathName_pre_per{97}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin6\';
handles.FileName_pre_per{97}='Slide1.sld - 2-pwer19-Frame2-reversespm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{98}='GRIN6';
handles.odor_pair{98}='MOHeptanal';
handles.session_no_per_mouse(98)=8;
handles.group(98)=2;
handles.sex(98)=1;
handles.PathName_pre_per{98}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin6\';
handles.FileName_pre_per{98}='Slide1.sld - 3-pwer19-Frame2-reversespm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{99}='GRIN6';
handles.odor_pair{99}='MOHeptanal';
handles.session_no_per_mouse(99)=8;
handles.group(99)=2;
handles.sex(99)=1;
handles.PathName_pre_per{99}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin6\';
handles.FileName_pre_per{99}='Slide1.sld - 4-pwer19-Frame2-reversespm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';



handles.mouse{100}='GRIN6';
handles.odor_pair{100}='MOHeptanal';
handles.session_no_per_mouse(100)=9;
handles.group(100)=2;
handles.sex(100)=1;
handles.PathName_pre_per{100}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220211-mmHippoFGrin6\';
handles.FileName_pre_per{100}='Slide1.sld - 1-pwer18-Frame2-reversespm-021122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{101}='GRIN6';
handles.odor_pair{101}='MOHeptanal';
handles.session_no_per_mouse(101)=9;
handles.group(101)=2;
handles.sex(101)=1;
handles.PathName_pre_per{101}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220211-mmHippoFGrin6\';
handles.FileName_pre_per{101}='Slide1.sld - 2-pwer18-Frame2-reversespm-021122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{102}='GRIN6';
handles.odor_pair{102}='MOHeptanal';
handles.session_no_per_mouse(102)=9;
handles.group(102)=2;
handles.sex(102)=1;
handles.PathName_pre_per{102}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220211-mmHippoFGrin6\';
handles.FileName_pre_per{102}='Slide1.sld - 3-pwer18-Frame2-reversespm-021122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{103}='GRIN6';
handles.odor_pair{103}='MOHeptanal';
handles.session_no_per_mouse(103)=9;
handles.group(103)=2;
handles.sex(103)=1;
handles.PathName_pre_per{103}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220211-mmHippoFGrin6\';
handles.FileName_pre_per{103}='Slide1.sld - 4-pwer18-Frame2-reversespm-021122_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{104}='GRIN7';
handles.odor_pair{104}='HeptanalMO';
handles.session_no_per_mouse(104)=1;
handles.group(104)=1;
handles.sex(104)=1;
handles.PathName_pre_per{104}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220124-mmHippoFGrin7\';
handles.FileName_pre_per{104}='Slide1.sld - 1-22pwer-Frame3-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{105}='GRIN7';
handles.odor_pair{105}='HeptanalMO';
handles.session_no_per_mouse(105)=1;
handles.group(105)=1;
handles.sex(105)=1;
handles.PathName_pre_per{105}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220124-mmHippoFGrin7\';
handles.FileName_pre_per{105}='Slide1.sld - 2-22pwer-Frame3-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{106}='GRIN7';
handles.odor_pair{106}='HeptanalMO';
handles.session_no_per_mouse(106)=1;
handles.group(106)=1;
handles.sex(106)=1;
handles.PathName_pre_per{106}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220124-mmHippoFGrin7\';
handles.FileName_pre_per{106}='Slide1.sld - 3-22pwer-Frame3-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{107}='GRIN7';
handles.odor_pair{107}='HeptanalMO';
handles.session_no_per_mouse(107)=2;
handles.group(107)=1;
handles.sex(107)=1;
handles.PathName_pre_per{107}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin7\';
handles.FileName_pre_per{107}='Slide2.sld - 1-pwer21-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{108}='GRIN7';
handles.odor_pair{108}='HeptanalMO';
handles.session_no_per_mouse(108)=2;
handles.group(108)=1;
handles.sex(108)=1;
handles.PathName_pre_per{108}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin7\';
handles.FileName_pre_per{108}='Slide2.sld - 2-pwer21-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{109}='GRIN7';
handles.odor_pair{109}='HeptanalMO';
handles.session_no_per_mouse(109)=2;
handles.group(109)=1;
handles.sex(109)=1;
handles.PathName_pre_per{109}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220126-mmHippoFGrin7\';
handles.FileName_pre_per{109}='Slide2.sld - 3-pwer21-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{110}='GRIN7';
handles.odor_pair{110}='HeptanalMO';
handles.session_no_per_mouse(110)=3;
handles.group(110)=1;
handles.sex(110)=1;
handles.PathName_pre_per{110}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin7\';
handles.FileName_pre_per{110}='Slide2.sld - 1-22pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{111}='GRIN7';
handles.odor_pair{111}='HeptanalMO';
handles.session_no_per_mouse(111)=3;
handles.group(111)=1;
handles.sex(111)=1;
handles.PathName_pre_per{111}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin7\';
handles.FileName_pre_per{111}='Slide2.sld - 2-22pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{112}='GRIN7';
handles.odor_pair{112}='HeptanalMO';
handles.session_no_per_mouse(112)=3;
handles.group(112)=1;
handles.sex(112)=1;
handles.PathName_pre_per{112}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin7\';
handles.FileName_pre_per{112}='Slide2.sld - 3-22pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{113}='GRIN7';
handles.odor_pair{113}='HeptanalMO';
handles.session_no_per_mouse(113)=3;
handles.group(113)=1;
handles.sex(113)=1;
handles.PathName_pre_per{113}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220128-mmHippoFGrin7\';
handles.FileName_pre_per{113}='Slide2.sld - 4-22pwer-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{114}='GRIN7';
handles.odor_pair{114}='HeptanalMO';
handles.session_no_per_mouse(114)=4;
handles.group(114)=1;
handles.sex(114)=1;
handles.PathName_pre_per{114}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin7\';
handles.FileName_pre_per{114}='slide1.sld - 1-pwer22-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{115}='GRIN7';
handles.odor_pair{115}='MOHeptanal';
handles.session_no_per_mouse(115)=4;
handles.group(115)=2;
handles.sex(115)=1;
handles.PathName_pre_per{115}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin7\';
handles.FileName_pre_per{115}='Slide1.sld - 2-pwer22-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{116}='GRIN7';
handles.odor_pair{116}='MOHeptanal';
handles.session_no_per_mouse(116)=4;
handles.group(116)=2;
handles.sex(116)=1;
handles.PathName_pre_per{116}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220131-mmHippoFGrin7\';
handles.FileName_pre_per{116}='Slide1.sld - 3-pwer22-Frame2-reversespm _E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{117}='GRIN7';
handles.odor_pair{117}='MOHeptanal';
handles.session_no_per_mouse(117)=5;
handles.group(117)=2;
handles.sex(117)=1;
handles.PathName_pre_per{117}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin7\';
handles.FileName_pre_per{117}='Slide2.sld - 1-pwer21-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{118}='GRIN7';
handles.odor_pair{118}='MOHeptanal';
handles.session_no_per_mouse(118)=5;
handles.group(118)=2;
handles.sex(118)=1;
handles.PathName_pre_per{118}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin7\';
handles.FileName_pre_per{118}='Slide2.sld - 2-pwer21-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{119}='GRIN7';
handles.odor_pair{119}='MOHeptanal';
handles.session_no_per_mouse(119)=5;
handles.group(119)=2;
handles.sex(119)=1;
handles.PathName_pre_per{119}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin7\';
handles.FileName_pre_per{119}='Slide2.sld - 3-pwer21-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{120}='GRIN7';
handles.odor_pair{120}='MOHeptanal';
handles.session_no_per_mouse(120)=5;
handles.group(120)=2;
handles.sex(120)=1;
handles.PathName_pre_per{120}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin7\';
handles.FileName_pre_per{120}='Slide2.sld - 4-pwer21-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{121}='GRIN7';
handles.odor_pair{121}='MOHeptanal';
handles.session_no_per_mouse(121)=5;
handles.group(121)=2;
handles.sex(121)=1;
handles.PathName_pre_per{121}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220203-mmHippoFGrin7\';
handles.FileName_pre_per{121}='Slide2.sld - 5-pwer21-Frame2-reversespm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


% This resulted in an error in drgCaImAn_SVZ_entire_session
% Error using ClassificationNaiveBayes\findNoDataCombos (line 347)
% A normal distribution cannot be fit for the combination of class 1 and
% predictor these_training_measurements. The data has zero variance.
handles.mouse{122}='GRIN7';
handles.odor_pair{122}='MOHeptanal';
handles.session_no_per_mouse(122)=6;
handles.group(122)=2;
handles.sex(122)=1;
handles.PathName_pre_per{122}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin7\';
handles.FileName_pre_per{122}='Slide2.sld - 1-pwer21-Frame2-reversespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{123}='GRIN7';
handles.odor_pair{123}='MOHeptanal';
handles.session_no_per_mouse(123)=6;
handles.group(123)=2;
handles.sex(123)=1;
handles.PathName_pre_per{123}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin7\';
handles.FileName_pre_per{123}='Slide2.sld - 2-pwer21-Frame2-reversespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{124}='GRIN7';
handles.odor_pair{124}='MOHeptanal';
handles.session_no_per_mouse(124)=6;
handles.group(124)=2;
handles.sex(124)=1;
handles.PathName_pre_per{124}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin7\';
handles.FileName_pre_per{124}='Slide2.sld - 3-pwer21-Frame2-reversespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{125}='GRIN7';
handles.odor_pair{125}='MOHeptanal';
handles.session_no_per_mouse(125)=6;
handles.group(125)=2;
handles.sex(125)=1;
handles.PathName_pre_per{125}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220204-mmHippoFGrin7\';
handles.FileName_pre_per{125}='Slide2.sld - 4-pwer21-Frame2-reversespm-020422_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{126}='GRIN7';
handles.odor_pair{126}='MOHeptanal';
handles.session_no_per_mouse(126)=7;
handles.group(126)=2;
handles.sex(126)=1;
handles.PathName_pre_per{126}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220207-mmHippoFGrin7\';
handles.FileName_pre_per{126}='Slide2.sld - 1-pwer21-Frame2-reversespm-020722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{127}='GRIN7';
handles.odor_pair{127}='MOHeptanal';
handles.session_no_per_mouse(127)=7;
handles.group(127)=2;
handles.sex(127)=1;
handles.PathName_pre_per{127}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220207-mmHippoFGrin7\';
handles.FileName_pre_per{127}='Slide2.sld - 2-pwer21-Frame2-reversespm-020722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{128}='GRIN7';
handles.odor_pair{128}='MOHeptanal';
handles.session_no_per_mouse(128)=7;
handles.group(128)=2;
handles.sex(128)=1;
handles.PathName_pre_per{128}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220207-mmHippoFGrin7\';
handles.FileName_pre_per{128}='Slide2.sld - 3-pwer21-Frame2-reversespm-020722_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.mouse{129}='GRIN7';
handles.odor_pair{129}='MOHeptanal';
handles.session_no_per_mouse(129)=8;
handles.group(129)=2;
handles.sex(129)=1;
handles.PathName_pre_per{129}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin7\';
handles.FileName_pre_per{129}='Slide2.sld - 1-pwer21-Frame2-reverspm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{130}='GRIN7';
handles.odor_pair{130}='MOHeptanal';
handles.session_no_per_mouse(130)=8;
handles.group(130)=2;
handles.sex(130)=1;
handles.PathName_pre_per{130}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin7\';
handles.FileName_pre_per{130}='Slide2.sld - 2-pwer21-Frame2-reverspm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';

handles.mouse{131}='GRIN7';
handles.odor_pair{131}='MOHeptanal';
handles.session_no_per_mouse(131)=8;
handles.group(131)=2;
handles.sex(131)=1;
handles.PathName_pre_per{131}='F:\SFTP\Ming Ma\Grin6and7_spm_nonneg\20220209-mmHippoFGrin7\';
handles.FileName_pre_per{131}='Slide2.sld - 3-pwer21-Frame2-reverspm-020922_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per.mat';


handles.no_files=length(handles.FileName_pre_per);
