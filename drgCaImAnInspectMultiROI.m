%drgCaImAnInspectMultiROI
%This function inspects the output of multi ROI decoding
clear all
close all
 
%Which ii_out and ROIs do you want to inspect
ii_out=1;
ROI1=50;
ROI2=62;

%Definition of other variables
figNo=0;
dt_span=15;
MLalgo=6;
time_windows=[3.1 4.1];

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

%Load file
[pre_perFileName,pre_perPathName] = uigetfile({'*pre_per_rdec.mat'},'Select the .m file with all the choices for analysis');

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAnInspectMultiROI run for ' pre_perFileName '\n\n']);



handles_out2=handles_out.ii_out(ii_out).handles_out;
dt=handles_out2.dt;
no_time_points=length(handles_out2.time_span);
no_ROI_draws=handles_out2.no_ROI_draws;
time_span=handles_out2.time_span;
dFF_per_trial_sm=handles_out2.dFF_per_trial_sm;
dFF_per_trial_sp=handles_out2.dFF_per_trial_sp;


%This is odor window
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

edges=[0:0.1:1];

hold on
%         accuracy_tr=handles_out2.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
accuracy_per_ROI=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI=[accuracy_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end

histogram(accuracy_per_ROI,edges,'Normalization','Probability')



accuracy_per_ROI_sh=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_sh=[accuracy_per_ROI_sh mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end

histogram(accuracy_per_ROI_sh,edges,'Normalization','Probability')


title(['Acccuracy histogram'])
xlabel('Accuracy')

fprintf(1,'Mean accuracy all runs= %d\n',mean(accuracy_per_ROI))
fprintf(1,'Shuffled trial mean accuracy all runs= %d\n',mean(accuracy_per_ROI_sh))

%Plot pseudocolor accuracy vs dFF fraction positive
%This is odor window
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on
colors=colormap(jet);

ii_t_start=find(time_span>=time_windows(1),1,'first');
ii_t_end=find(time_span<=time_windows(2),1,'last');

for iiROI=1:no_ROI_draws

    for ii_t=ii_t_start:ii_t_end
        dFF_sm_nonzero=0;
        n_dFF_sm=0;
        for trNo=1:size(dFF_per_trial_sm,1)
            n_dFF_sm=n_dFF_sm+1;
            if dFF_per_trial_sm(trNo,iiROI,ii_t)>0.05
                dFF_sm_nonzero=dFF_sm_nonzero+1;
            end
        end
        dFF_sm_nonzero_f=dFF_sm_nonzero/n_dFF_sm;

        dFF_sp_nonzero=0;
        n_dFF_sp=0;
        for trNo=1:size(dFF_per_trial_sp,1)
            n_dFF_sp=n_dFF_sp+1;
            if dFF_per_trial_sp(trNo,iiROI,ii_t)>0.05
                dFF_sp_nonzero=dFF_sp_nonzero+1;
            end
        end
        dFF_sp_nonzero_f=dFF_sp_nonzero/n_dFF_sp;
        if ((accuracy_per_ROI(iiROI)<0.3)||(accuracy_per_ROI(iiROI)>0.7))
            ii_color=1+ceil(accuracy_per_ROI(iiROI)*255);
            plot(dFF_sm_nonzero_f,dFF_sp_nonzero_f,'o','MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                'MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                'MarkerSize',5)
        end
    end
end

title('Pseudocolor accuracy')
xlabel('Fraction dFF non zero S-')
ylabel('Fraction dFF non zero S+')

figNo=figNo+1;

try
    close(figNo)
catch
end


%Plot the per-trial timecourse
hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.83 .1 .05 .3])

prain=[0:1/99:1];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap jet
shading interp
ax=gca;
set(ax,'XTickLabel','')



%Plot prediction

%     for iiMLalgo=MLalgo_to_use
%This is odor window

%Do Sp first
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on

edges=[0:0.1:1];

%S+
prediction_sp_per_ROI=[];
for iiROI=1:no_ROI_draws
    prediction_sp_per_ROI=[prediction_sp_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end

histogram(prediction_sp_per_ROI, edges)

%S-
prediction_sm_per_ROI=[];
for iiROI=1:no_ROI_draws
    prediction_sm_per_ROI=[prediction_sm_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end

histogram(prediction_sm_per_ROI, edges)


title(['Histogram of label predictions'])
xlabel('Prediction')
legend('S+','S-')


%Plot prediction S+ vs S-
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on


for iiROI=1:no_ROI_draws
    ii_color=1+ceil(accuracy_per_ROI(iiROI)*255);
    plot(prediction_sm_per_ROI(iiROI), prediction_sp_per_ROI(iiROI),'o',...
        'MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
        'MarkerEdgeColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
        'MarkerSize',5)
end


title(['Prediction S+ vs. S-, color accuracy'])
ylabel('Prediction S+')
xlabel('Prediction S-')

%Plot the accuracy for ROI1
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on

this_correct_predict_sh=handles_out2.ROI(ROI1).MLalgo(MLalgo).this_correct_predict_sh;
this_correct_predict=handles_out2.ROI(ROI1).MLalgo(MLalgo).this_correct_predict;

CIsm = bootci(1000, @mean, this_correct_predict_sh);
meansm=mean(this_correct_predict_sh,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');


CIsm = bootci(1000, @mean, this_correct_predict);
meansm=mean(this_correct_predict,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);

plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);

text(30,0.75,'Shuffled','Color','k')
text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])


ylim([-0.1 1.1])
xlim([-7 15])
this_ylim=ylim;

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor on markers
plot([0 0],this_ylim,'-k')
odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

title(['Accuracy for ii_ROI '  num2str(ROI1)])
xlabel('Time(sec)')
ylabel('Accuracy')

%Now plot the accuracy for ROI2
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on

this_correct_predict_sh=handles_out2.ROI(ROI2).MLalgo(MLalgo).this_correct_predict_sh;
this_correct_predict=handles_out2.ROI(ROI2).MLalgo(MLalgo).this_correct_predict;

CIsm = bootci(1000, @mean, this_correct_predict_sh);
meansm=mean(this_correct_predict_sh,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');


CIsm = bootci(1000, @mean, this_correct_predict);
meansm=mean(this_correct_predict,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);

plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);

text(30,0.75,'Shuffled','Color','k')
text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])


ylim([-0.1 1.1])
xlim([-7 15])
this_ylim=ylim;

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor on markers
plot([0 0],this_ylim,'-k')
odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

title(['Accuracy for ii_ROI '  num2str(ROI2)])
xlabel('Time(sec)')
ylabel('Accuracy')

%Calculate accuracy in a wider post-stimulus window and latency

%This is post-stimulus window
time_windows=[-1.5 10];
pre_time_window=[-3 -1.5];

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

edges=[0:0.1:1];

hold on
%         accuracy_tr=handles_out2.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
accuracy_per_ROIw=[];
accuracy_per_ROIpre=[];
accuracy_per_ROIstd=[];
latency_per_ROI=[];
cropped_time_span_ii=find((time_span>=time_windows(1))&(time_span<=time_windows(2)));
cropped_time_span=time_span(cropped_time_span_ii);
for iiROI=1:no_ROI_draws
    accuracy_per_ROIw=[accuracy_per_ROIw prctile(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),1),95)];
    accuracy_per_ROIpre=[accuracy_per_ROIpre mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),1))];
    accuracy_per_ROIstd=[accuracy_per_ROIstd std(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),1))];
    this_acc_timecourse=mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),1);
    this_latency_ii=find(this_acc_timecourse>accuracy_per_ROIpre(end)+2*accuracy_per_ROIstd(end),1,'first');
    if ~isempty(this_latency_ii)
        latency_per_ROI=[latency_per_ROI cropped_time_span(this_latency_ii)];
    else
        latency_per_ROI=[latency_per_ROI NaN];
    end
end

histogram(accuracy_per_ROIw,edges,'Normalization','Probability')

histogram(accuracy_per_ROI_sh,edges,'Normalization','Probability')


title(['Acccuracy histogram, wide window'])
xlabel('Accuracy')

fprintf(1,'Mean accuracy wide window all runs= %d\n',mean(accuracy_per_ROIw))
fprintf(1,'Shuffled trial mean accuracy all runs= %d\n',mean(accuracy_per_ROI_sh))


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

edges=[-1.5:1:10.5];
histogram(latency_per_ROI(accuracy_per_ROIw>0.65),edges)

title(['Histogram for latency for ROIs with accuracy per ROIw > 0.65'])
xlabel('Time(sec)')
ylabel('Counts')

%Plot the accuracy for ROIs with accuracy_per_ROIw>0.65
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .6])

hold on

ii_shift=0;

for iiROI=1:no_ROI_draws
    if accuracy_per_ROIw(iiROI)>0.65
        ii_shift=ii_shift+1;
        this_correct_predict_sh=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh;
        this_correct_predict=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict;

        plot(time_span',mean(this_correct_predict_sh,1)'+ii_shift,'-k','DisplayName','Shuffled')
        plot(time_span',mean(this_correct_predict,1)'+ii_shift, '-','Color',[0 114/255 178/255]);
    end
end

text(30,0.75,'Shuffled','Color','k')
text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])


ylim([-0.1 ii_shift+1.1])
xlim([-7 15])
this_ylim=ylim;

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor on markers
plot([0 0],this_ylim,'-k')
odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

title(['Accuracy for ROIs with accuracy per ROIw > 0.65'])
xlabel('Time(sec)')
ylabel('Accuracy')

pfft=1; 



