function handles=drgCaImAn_dropc_choices_Grin

%run with drgCalmAn_batch_dropc_fsdz

%%Grin1 fsds1
%This is the experiment of own home odor (as S+) x other mouse home odor (as S-)
%no lick or water reward given, no water restriction 

handles.first_file=1;


%csv or .mat extract files for image analysis
% handles.PathNamecsv{1}='/Volumes/Diego Mac Drive/SFTP/CalmAn_20210302/';
% handles.csvFileName{1}='Results_Grin3_homeodor_ZX1614725527_bestROIs.csv';

handles.PathNamecsv{1}='E:\SFTP\ncorr_ext\20210621\';
handles.csvFileName{1}='20210621_Grin4_fsds_home_pcdh2_XY1624311429_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{2}='E:\SFTP\ncorr_ext\20210617\';
handles.csvFileName{2}='Grin1_fsds_home_otherPcdh2_XY1623961175_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{3}='E:\SFTP\ncorr_ext\20210803\';
handles.csvFileName{3}='20210803_FCM2_Home_OtherGrin1_2_XY1628020657_Z0_T00000_C0_ncorr_ext.mat'
handles.PathNamecsv{4}='E:\SFTP\ncorr_ext\20210422\';
handles.csvFileName{4}='20210422_Grin4_home_home_control_XY1619120309_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{5}='E:\SFTP\ncorr_ext\20210420\';
handles.csvFileName{5}='20210420_Grin3_home_home_control_XY1618950608_Z0_T00000_C0_ncorr_ext.mat'
handles.PathNamecsv{6}='E:\SFTP\ncorr_ext\20210415\';
handles.csvFileName{6}='20210415_Grin1_home_home_control_XY1618521884_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{7}='E:\SFTP\ncorr_ext\20210409\';
handles.csvFileName{7}='20210409_Grin4_home_male_XY1617998540_Z0_T00000_C0_ncorr_ext.mat'
handles.PathNamecsv{8}='E:\SFTP\ncorr_ext\20210330\';
handles.csvFileName{8}='20210330_Grin3_fsds_Home_otherMale_XY1617138791_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{9}='E:\SFTP\ncorr_ext\20210408\';
handles.csvFileName{9}='20210408_Grin3_Home_othermale_XY1617911391_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{10}='E:\SFTP\ncorr_ext\20210415\';
handles.csvFileName{10}='20210415_Grin1_home_othermale_XY1618518579_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{11}='E:\SFTP\ncorr_ext\20210402\';
handles.csvFileName{11}='20210402_Grin1_fsds_home_othermale_XY1617394264_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{12}='E:\SFTP\ncorr_ext\20210409\';
handles.csvFileName{12}='20210409_Grin4_home_otherGRIN3_KX_XY1618002992_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{13}='E:\SFTP\ncorr_ext\20210408\';
handles.csvFileName{13}='20210408_Grin3_Home_otheGRIN4_KX_XY1617916254_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{14}='E:\SFTP\ncorr_ext\20210317\';
handles.csvFileName{14}='20210317_Grin1_fsds_Home_OtherGrin3_kxanestesia_XY1616018201_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{15}='E:\SFTP\ncorr_ext\20210308\';
handles.csvFileName{15}='20210308_Grin4_fsds_homeodor_otherGRIN3_XY1615240058_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{16}='E:\SFTP\ncorr_ext\20210302\';
handles.csvFileName{16}='20210302_Grin1_homeodor_otherGRIN3_XY1614719876_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{17}='E:\SFTP\ncorr_ext\20210302\';
handles.csvFileName{17}='20210302_Grin3_homeodor_otherGRIN1_2_XY1614725527_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{18}='E:\SFTP\ncorr_ext\20210323\';
handles.csvFileName{18}='20210323_Grin3_fsds_Amyl_Acethop_XY1616532966_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{19}='E:\SFTP\ncorr_ext\20210317\';
handles.csvFileName{19}='20210317_Grin1_fsds_AmylSp_AcetopheSm_XY0_Z0_T00000_C0_ncorr_ext.mat';
handles.PathNamecsv{20}='E:\SFTP\ncorr_ext\20210308\';
handles.csvFileName{20}='20210308_Grin4_fsds_amyl_acetophenone_XY1615243443_Z0_T0000_C0_ncorr_ext.mat';

%This is the time per image in seconds
handles.dt(1)=0.31494;
handles.dt(2)=0.31587;
handles.dt(3)=0.31362;
handles.dt(4)=0.31347;
handles.dt(5)=0.31450;
handles.dt(6)=0.31360;
handles.dt(7)=0.31369;
handles.dt(8)=0.31411;
handles.dt(9)=0.31567;
handles.dt(10)=0.31357;
handles.dt(11)=0.31398;
handles.dt(12)=0.31344;
handles.dt(13)=0.31499;
handles.dt(14)=0.31429;
handles.dt(15)=0.31336;
handles.dt(16)=0.31314;
handles.dt(17)=0.31333;
handles.dt(18)=mean([0.31333 0.31364]);
handles.dt(19)=0.31364;
handles.dt(20)=0.31405;

%Other files
handles.PathName{1}='E:\SFTP\ncorr_ext\20210621\';
handles.PathName{2}='E:\SFTP\ncorr_ext\20210617\';
handles.PathName{3}='E:\SFTP\ncorr_ext\20210803\';
handles.PathName{4}='E:\SFTP\ncorr_ext\20210422\';
handles.PathName{5}='E:\SFTP\ncorr_ext\20210420\';
handles.PathName{6}='E:\SFTP\ncorr_ext\20210415\';
handles.PathName{7}='E:\SFTP\ncorr_ext\20210409\';
handles.PathName{8}='E:\SFTP\ncorr_ext\20210330\';
handles.PathName{9}='E:\SFTP\ncorr_ext\20210408\';
handles.PathName{10}='E:\SFTP\ncorr_ext\20210415\';
handles.PathName{11}='E:\SFTP\ncorr_ext\20210402\';
handles.PathName{12}='E:\SFTP\ncorr_ext\20210409\';
handles.PathName{13}='E:\SFTP\ncorr_ext\20210408\';
handles.PathName{14}='E:\SFTP\ncorr_ext\20210317\';
handles.PathName{15}='E:\SFTP\ncorr_ext\20210308\';
handles.PathName{16}='E:\SFTP\ncorr_ext\20210302\';
handles.PathName{17}='E:\SFTP\ncorr_ext\20210302\';
handles.PathName{18}='E:\SFTP\ncorr_ext\20210323\';
handles.PathName{19}='E:\SFTP\ncorr_ext\20210317\';
handles.PathName{20}='E:\SFTP\ncorr_ext\20210308\';

%spm filesdr
handles.spmFileName{1}='20210621_Grin4_fsds_Home_pcdh2_210621_154452_210621_154818spm.mat';
handles.spmFileName{2}='20210617_Grin1_fsds_home_otherPcdh2_20210617T143117spm.mat';
handles.spmFileName{3}='20210803_FCM2_home_otherGrin1_210803_140922_2spm.mat';
handles.spmFileName{4}='20210422_Grin4_home_home_conrol_210422_134926.mat';
handles.spmFileName{5}='20210420_Grin3_home_home_control_210420_144055.mat';
handles.spmFileName{6}='20210415_Grin1_home_home_control_210415_153534.mat';
handles.spmFileName{7}='20210409_Grin4_fsds_home_othermale_210409_141255.mat';
handles.spmFileName{8}='20210330_Grin3_fsds_homeodor_otherMales_210330_152355.mat';
handles.spmFileName{9}='20210408_Grin3_home_othermale_210408_140028.mat';
handles.spmFileName{10}='20210415_Grin1_home_othermale_210415_144024.mat';
handles.spmFileName{11}='20210402_Grin1_FSDS_Home_othermale_210402_142128.mat';
handles.spmFileName{12}='20210409_Grin4_fsds_home_otherGRIN3_KX_210409_152715.mat';
handles.spmFileName{13}='20210408_Grin3_home_otherGRIN4_KX_210408_152137.mat';
handles.spmFileName{14}='20210317_Grin1_fsds_Home_OtherGrin3_kxanestesia_210317_160722.mat';
handles.spmFileName{15}='20210308_Grin4_fsds_homeodor_210308_145810.mat';
handles.spmFileName{16}='20210302_Grin1_Homeodor_otherodor.mat';
handles.spmFileName{17}='20210302_Grin3_Homeodor_otherodor_210302_160234_2.mat';
handles.spmFileName{18}='20210323_Grin3_FSDS_AmylAcetate_Acetophenone_210323_150624.mat';
handles.spmFileName{19}='20210317_Grin1_fsds_AmylAcetate_Acetophenone_210317_144542.mat';
handles.spmFileName{20}='20210308_Grin4_fsds_amyl_acetophenone_210308_155437.mat';

%rhd files
handles.rhdFileName{1}='20210621_Grin4_fsds_Home_pcdh2_210621_154452_210621_154818.rhd';
handles.rhdFileName{2}='20210617_Grin1_fsds_home_otherPcdh2_20210617T143117spm.rhd';
handles.rhdFileName{3}='20210803_FCM2_home_otherGrin1_210803_140922_2.rhd';
handles.rhdFileName{4}='20210422_Grin4_home_home_conrol_210422_134926.rhd';
handles.rhdFileName{5}='20210420_Grin3_home_home_control_210420_144055.rhd';
handles.rhdFileName{6}='20210415_Grin1_home_home_control_210415_153534.rhd';
handles.rhdFileName{7}='20210409_Grin4_fsds_home_othermale_210409_141255.rhd';
handles.rhdFileName{8}='20210330_Grin3_fsds_homeodor_otherMales_210330_152355.rhd';
handles.rhdFileName{9}='20210408_Grin3_home_othermale_210408_140028.rhd';
handles.rhdFileName{10}='20210415_Grin1_home_othermale_210415_144024.rhd';
handles.rhdFileName{11}='20210402_Grin1_FSDS_Home_othermale_210402_142128.rhd';
handles.rhdFileName{12}='20210409_Grin4_fsds_home_otherGRIN3_KX_210409_152715.rhd';
handles.rhdFileName{13}='20210408_Grin3_home_otherGRIN4_KX_210408_152137.rhd';
handles.rhdFileName{14}='20210317_Grin1_fsds_Home_OtherGrin3_kxanestesia_210317_160722.rhd';
handles.rhdFileName{15}='20210308_Grin4_fsds_homeodor_210308_145810.rhd';
handles.rhdFileName{16}='20210302_Grin1_Homeodor_otherodor_210302_142808.rhd';
handles.rhdFileName{17}='20210302_Grin3_Homeodor_otherodor_210302_160234_2.rhd';
handles.rhdFileName{18}='20210323_Grin3_FSDS_AmylAcetate_Acetophenone_210323_150624.rhd';
handles.rhdFileName{19}='20210317_Grin1_fsds_AmylAcetate_Acetophenone_210317_144542.rhd';
handles.rhdFileName{20}='20210308_Grin4_fsds_amyl_acetophenone_210308_155437.rhd';

handles.no_files=length(handles.csvFileName);
