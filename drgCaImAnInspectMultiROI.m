clear all
close all
 
%Which ii_out do you want to inspect
ii_out=1;

figNo=0;
dt_span=15;
MLalgo=6;
time_windows=[2 4.1];

%Load file
[pre_perFileName,pre_perPathName] = uigetfile({'*pre_per_rdec.mat'},'Select the .m file with all the choices for analysis');

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);



handles_out2=handles_out.ii_out(ii_out).handles_out;
dt=handles_out2.dt;
no_time_points=size(handles_out2.ROI(1).MLalgo(MLalgo).per_trial_sp_timecourse,2);
no_ROI_draws=length(handles_out2.ROI);
time_span=[0:dt:dt*no_time_points]-dt_span+dt;

%This is odor window
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on
%         accuracy_tr=handles_out2.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
accuracy_per_ROI=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI=[accuracy_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end

histogram(accuracy_per_ROI)



accuracy_per_ROI_sh=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_sh=[accuracy_per_ROI_sh mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end

histogram(accuracy_per_ROI_sh)


title(['Acccuracy histogram'])
xlabel('Accuracy')

fprintf(1,'Mean accuracy all runs= %d\n',mean(accuracy_per_ROI))
fprintf(1,'Shuffled trial mean accuracy all runs= %d\n',mean(accuracy_per_ROI_sh))

pfft=1;



