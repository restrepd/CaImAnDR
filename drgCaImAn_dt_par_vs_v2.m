% drgCaImAn_dt_par_vs_v2.m
% Compares timecourses for accuracy between per image TTL-synchronized and
%single TTL synchronized
close all
clear all

load('F:\SFTP\Ming Ma\CamKGrin1\drgCaImAn_LDAfsdz_choices_nn_v2_mmHippoCamKGrin1_08_24_22.mat')
accuracy_out_v2=accuracy_out;


load('F:\SFTP\Ming Ma\CamKGrin1\drgCaImAn_LDAfsdz_choices_nn_par_mmHippoCamKGrin1_08_18_22.mat')
accuracy_out_par=accuracy_out;

pfft=1;

figureNo=0;
figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

%get the v2
these_per_corr_v2=[];
grNo=4;
for ii=1:size(accuracy_out_v2.groupNo(grNo).these_per_corr,1)
    these_per_corr_v2(ii,:)=accuracy_out_v2.groupNo(grNo).these_per_corr(ii,:);
end

ii_end=size(accuracy_out_v2.groupNo(grNo).these_per_corr,1);

grNo=8;
for ii=1:size(accuracy_out_v2.groupNo(grNo).these_per_corr,1)
    these_per_corr_v2(ii+ii_end,:)=accuracy_out_v2.groupNo(grNo).these_per_corr(ii,:);
end

time_span_v2=accuracy_out_v2.groupNo(grNo).time_span;


CIpv = bootci(1000, @mean, these_per_corr_v2);
meanpv=mean(these_per_corr_v2,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;


[hlpvl, hppvl] = boundedline(time_span_v2',mean(these_per_corr_v2,1)', CIpv', 'r');

%get the par
these_per_corr_par=[];
grNo=4;
for ii=1:size(accuracy_out_par.groupNo(grNo).these_per_corr,1)
    these_per_corr_par(ii,:)=accuracy_out_par.groupNo(grNo).these_per_corr(ii,:);
end

ii_end=size(accuracy_out_par.groupNo(grNo).these_per_corr,1);

grNo=8;
for ii=1:size(accuracy_out_par.groupNo(grNo).these_per_corr,1)
    these_per_corr_par(ii+ii_end,:)=accuracy_out_par.groupNo(grNo).these_per_corr(ii,:);
end
 

CIpv = bootci(1000, @mean, these_per_corr_par);
meanpv=mean(these_per_corr_par,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;

time_span_par=accuracy_out_par.groupNo(grNo).time_span;

[hlpvl, hppvl] = boundedline(time_span_par',mean(these_per_corr_par,1)', CIpv', 'b');



plot(time_span_v2',mean(these_per_corr_v2,1)', 'r','LineWidth',1.5);
plot(time_span_par',mean(these_per_corr_par,1)', 'b','LineWidth',1.5);

ylim([0.3 1.2])
this_ylim=ylim;
plot([0 0],this_ylim,'-k')
xlim([-10 20])

xlabel('Time(sec)')


title(['Decoding accuracy for red: singleTTL, blue: per image TTL'])

%Determine the time shift by minimizing the absolute difference from -1 to
%2 sec
mean_per_corr_v2=mean(these_per_corr_v2,1);


mean_per_corr_par=mean(these_per_corr_par,1);




%Now find the optimal shift
abs_dpcorr=[];
t_shifts=[];
ii=0;
for t_shift=0:0.05:2
    this_mean_per_corr_v2=mean_per_corr_v2((time_span_v2>=1-t_shift)&(time_span_v2<=2-t_shift));
    this_mean_per_corr_par=mean_per_corr_par((time_span_par>=1)&(time_span_par<=2));
    ii=ii+1;
    t_shifts(ii)=t_shift;
    abs_dpcorr(ii)=sum(abs(this_mean_per_corr_v2-this_mean_per_corr_par));
end

[minabs,minii]=min(abs_dpcorr);
optimal_t_shift=t_shifts(minii);

fprintf(1, ['\nOptimal time shift from fit = ' num2str(t_shifts(minii)) '\n\n']);


figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

%get the v2



CIpv = bootci(1000, @mean, these_per_corr_v2);
meanpv=mean(these_per_corr_v2,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;


[hlpvl, hppvl] = boundedline(time_span_v2'+optimal_t_shift,mean(these_per_corr_v2,1)', CIpv', 'r');



CIpv = bootci(1000, @mean, these_per_corr_par);
meanpv=mean(these_per_corr_par,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;

time_span_par=accuracy_out_par.groupNo(grNo).time_span;

[hlpvl, hppvl] = boundedline(time_span_par',mean(these_per_corr_par,1)', CIpv', 'b');



plot(time_span_v2'+optimal_t_shift,mean(these_per_corr_v2,1)', 'r','LineWidth',1.5);
plot(time_span_par',mean(these_per_corr_par,1)', 'b','LineWidth',1.5);

ylim([0.3 1.2])
this_ylim=ylim;
plot([0 0],this_ylim,'-k')
xlim([-10 20])

xlabel('Time(sec)')
 

title(['Fit-shifted decoding accuracy for red: singleTTL, blue: per image TTL'])

 
%Now get the lag between the first image TTL and the first FV in the
%metadata
choiceBatchPathName='F:\SFTP\Ming Ma\CamKGrin1\';
choiceFileName='drgCaImAn_LDAfsdz_choices_nn_par_mmHippoCamKGrin1_08_18_22.m';
addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;



%Parallel batch processing for each file
t_shift_rhd=[];
for fileNo=1:handles_choice.no_files
    pre_per_PathName=handles_choice.PathName_pre_per{fileNo};
    pre_per_FileName=handles_choice.FileName_pre_per{fileNo};
    load([pre_per_PathName  pre_per_FileName])
    t_shift_rhd(fileNo)=first_digital_in_time_rhd-first_imge_ttl_time_rhd;
    pffft=1;
end


fprintf(1, ['\nOptimal time shift from rhd = ' num2str(mean(t_shift_rhd)) '\n\n']);

optimal_t_shift=mean(t_shift_rhd);

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

%get the v2



CIpv = bootci(1000, @mean, these_per_corr_v2);
meanpv=mean(these_per_corr_v2,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;


[hlpvl, hppvl] = boundedline(time_span_v2'+optimal_t_shift,mean(these_per_corr_v2,1)', CIpv', 'r');



CIpv = bootci(1000, @mean, these_per_corr_par);
meanpv=mean(these_per_corr_par,1);
CIpv(1,:)=meanpv-CIpv(1,:);
CIpv(2,:)=CIpv(2,:)-meanpv;

time_span_par=accuracy_out_par.groupNo(grNo).time_span;

[hlpvl, hppvl] = boundedline(time_span_par',mean(these_per_corr_par,1)', CIpv', 'b');



plot(time_span_v2'+optimal_t_shift,mean(these_per_corr_v2,1)', 'r','LineWidth',1.5);
plot(time_span_par',mean(these_per_corr_par,1)', 'b','LineWidth',1.5);

ylim([0.3 1.2])
this_ylim=ylim;
plot([0 0],this_ylim,'-k')
xlim([-10 20])

xlabel('Time(sec)')
 

title(['rhd-shifted decoding accuracy for red: singleTTL, blue: per image TTL'])
pffft=1;


