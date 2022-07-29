function handles=drgCaImAn_LDAfsdz_choices_deepcad_raw_entire_02102022

%run with drgCalmAn_batch_dropc_fsdz
% handles.PathName_out='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
% handles.FileName_out='entire_session_deepcad6_02052022.mat';

%Path and suffix for file output
handles.PathName_out='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.suffix_out='_dec2.mat';

%Set the cost

handles.processing_algorithm=3;

handles.MLalgo_to_use=[1 4 6]; %Compare algorithms
handles.post_time=[5];
handles.k_fold=[5]; %Not used for processing_algorithm=3
handles.post_shift=[0];
handles.pre_time=[5];

handles.p_threshold=[0.1 0.2 0.3 1.1]; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
handles.dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
handles.show_figures=0; %Show the figures
handles.ii_cost=3; %Cost for decoding models

handles.first_file=1;

handles.group_names{1}='AAAP';
handles.group_names{2}='homeodor_female';
handles.group_names{3}='homeodor_female(kx_anesthesia)';
handles.group_names{4}='homeodor_male';
handles.group_names{5}='homeodor_sameodor';
handles.group_names{6}='odor_off';

% handles.group_names{6}='homeodor_male_pcdh21'; %Note: I am doing pcdh1
% and male together

%Sex
% 0 Male
% 1 Female

handles.mouse{1}='GRIN4';
handles.odor_pair{1}='homeodor_malepcdh2';
handles.group(1)=4;
handles.sex(1)=1;
handles.PathName_pre_per{1}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{1}='20210621_Grin4_fsds_home_pcdh2_XY1624311429_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';


handles.mouse{2}='GRIN1';
handles.odor_pair{2}='homeodor_pcdh2';
handles.group(2)=4;
handles.sex(2)=1;
handles.PathName_pre_per{2}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{2}='Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{3}='GRIN7';
handles.odor_pair{3}='homeodor_female';
handles.group(3)=2;
handles.sex(3)=1;
handles.PathName_pre_per{3}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{3}='20211214_Grin7_home_otherFemaleGRIN6_XY1639506542_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{4}='GRIN4';
handles.odor_pair{4}='homeodor_male';
handles.group(4)=4;
handles.sex(4)=1;
handles.PathName_pre_per{4}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{4}='20210409_Grin4_home_male_XY1617998540_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{5}='GRIN3';
handles.odor_pair{5}='homeodor_male';
handles.group(5)=4;
handles.sex(5)=1;
handles.PathName_pre_per{5}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{5}='20210330_Grin3_fsds_Home_otherMale_XY1617138791_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{6}='GRIN3';
handles.odor_pair{6}='homeodor_male';
handles.group(6)=4;
handles.sex(6)=1;
handles.PathName_pre_per{6}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{6}='20210408_Grin3_Home_othermale_XY1617911391_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{7}='GRIN1';
handles.odor_pair{7}='homeodor_male';
handles.group(7)=4;
handles.sex(7)=1;
handles.PathName_pre_per{7}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{7}='20210415_Grin1_home_othermale_XY1618518579_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{8}='GRIN1';
handles.odor_pair{8}='homeodor_male';
handles.group(8)=4;
handles.sex(8)=1;
handles.PathName_pre_per{8}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{8}='20210402_Grin1_fsds_home_othermale_XY1617394264_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{9}='GRIN6';
handles.odor_pair{9}='homeodor_male';
handles.group(9)=4;
handles.sex(9)=1;
handles.PathName_pre_per{9}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{9}='20211213_Grin6_Home_OtherMalePCDH2_XY1639425701_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';


handles.mouse{10}='GRIN7';
handles.odor_pair{10}='homeodor_male';
handles.group(10)=4;
handles.sex(10)=1;
handles.PathName_pre_per{10}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{10}='20211213_Grin7_Home_OtherMalePCDH2_XY1639434384_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{11}='GRIN6';
handles.odor_pair{11}='homeodor_male';
handles.group(11)=4;
handles.sex(11)=1;
handles.PathName_pre_per{11}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{11}='20211214_Grin6_home_othermalePCDH2_XY1639514792_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{12}='GRIN7';
handles.odor_pair{12}='homeodor_male';
handles.group(12)=4;
handles.sex(12)=1;
handles.PathName_pre_per{12}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{12}='20211214_Grin7_home_otherFemaleGRIN6_XY1639506542_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';


% handles.mouse{13}='GRIN1';
% handles.odor_pair{13}='homeodor_female';
% handles.group(13)=2;
% handles.sex(13)=1;
% handles.PathName_pre_per{13}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
% handles.FileName_pre_per{13}='20210803_FCM2_Home_OtherGrin1_2_XY1628020657_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{13}='GRIN6';
handles.odor_pair{13}='AAAP';
handles.group(13)=1;
handles.sex(13)=1;
handles.PathName_pre_per{13}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{13}='20211208_Grin6_ISO_ACT10p_XY1639001826_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{14}='GRIN4';
handles.odor_pair{14}='homeodor_sameodor';
handles.group(14)=5;
handles.sex(14)=1;
handles.PathName_pre_per{14}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{14}='20210422_Grin4_home_home_control_XY1619120309_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{15}='GRIN3';
handles.odor_pair{15}='homeodor_sameodor';
handles.group(15)=5;
handles.sex(15)=1;
handles.PathName_pre_per{15}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{15}='20210420_Grin3_home_home_control_XY1618950608_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{16}='GRIN1';
handles.odor_pair{16}='homeodor-sameodor';
handles.group(16)=5;
handles.sex(16)=1;
handles.PathName_pre_per{16}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{16}='20210415_Grin1_home_home_control_XY1618521884_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{17}='GRIN7';
handles.odor_pair{17}='AAAP';
handles.group(17)=1;
handles.sex(17)=1;
handles.PathName_pre_per{17}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{17}='20211208_Grin7_ISO_ACET10p_XY1638992533_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

% 
% handles.mouse{17}='GRIN1';
% handles.odor_pair{17}='homeodor_female(kx_anesthesia)';
% handles.group(17)=3;
% handles.sex(17)=1;
% handles.PathName_pre_per{17}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
% handles.FileName_pre_per{17}='20210317_Grin1_fsds_Home_OtherGrin3_kxanestesia_XY1616018201_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{18}='GRIN6';
handles.odor_pair{18}='homeodor_female';
handles.group(18)=2;
handles.sex(18)=1;
handles.PathName_pre_per{18}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{18}='20211213_Grin6_Home_OtherFemaleGRIN7_XY1639429269_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{19}='GRIN4';
handles.odor_pair{19}='homeodor_female';
handles.group(19)=2;
handles.sex(19)=1;
handles.PathName_pre_per{19}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{19}='20210308_Grin4_fsds_homeodor_otherGRIN3_XY1615240058_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{20}='GRIN1';
handles.odor_pair{20}='homeodor_female';
handles.group(20)=2;
handles.sex(20)=1;
handles.PathName_pre_per{20}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{20}='20210302_Grin1_homeodor_otherGRIN3_XY1614719876_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';


% This resulted in an error in drgCaImAn_SVZ_entire_session
% Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
% A normal distribution cannot be fit for the combination of class 1 and
% predictor these_training_measurements. The data has zero variance.
handles.mouse{21}='GRIN4';
handles.odor_pair{21}='homeodor_female(kx_anesthesia)';
handles.group(21)=3;
handles.sex(21)=1;
handles.PathName_pre_per{21}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{21}='20210409_Grin4_home_otherGRIN3_KX_XY1618002992_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{22}='GRIN3';
handles.odor_pair{22}='AAAP';
handles.group(22)=1;
handles.sex(22)=1;
handles.PathName_pre_per{22}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{22}='20210323_Grin3_fsds_Amyl_Acethop_XY1616532966_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{23}='GRIN1';
handles.odor_pair{23}='AAAP';
handles.group(23)=1;
handles.sex(23)=1;
handles.PathName_pre_per{23}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{23}='20210317_Grin1_fsds_AmylSp_AcetopheSm_XY0_Z0_T00000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{24}='GRIN4';
handles.odor_pair{24}='AAAP';
handles.group(24)=1;
handles.sex(24)=1;
handles.PathName_pre_per{24}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{24}='20210308_Grin4_fsds_amyl_acetophenone_XY1615243443_Z0_T0000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{25}='GRIN4';
handles.odor_pair{25}='AAAP';
handles.group(25)=1;
handles.sex(25)=1;
handles.PathName_pre_per{25}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{25}='20210228_Grin4_fsds1_amyl_acetophen_XY1614548423_Z0_T0000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{26}='GRIN7';
handles.odor_pair{26}='AAAP';
handles.group(26)=1;
handles.sex(26)=1;
handles.PathName_pre_per{26}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{26}='20211206_Grin7_ISO_ACT10p_XY1638830043_Z0_T0000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{27}='GRIN2';
handles.odor_pair{27}='odor_off';
handles.group(27)=6;
handles.sex(27)=1;
handles.PathName_pre_per{27}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{27}='20201119_Grin2_Ca1_spm3_XY1605816795_Z0_T0000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';

handles.mouse{28}='GRIN2';
handles.odor_pair{28}='odor_off';
handles.group(28)=6;
handles.sex(28)=1;
handles.PathName_pre_per{28}='E:\SFTP\Two_odor_DeepCAD_E_10_Iter_3500\';
handles.FileName_pre_per{28}='20201119_Grin2_Ca1_spm1_XY1605812694_Z0_T0000_C0f_E_10_Iter_3500_output_ncorr_ext_batchv2_pre_per.mat';


handles.no_files=length(handles.FileName_pre_per);
