function handles_outd=drgCaImAnInspectDecodingSelectROI(handles_choices)
%drgCaImAnInspectDecodingSingleROI
%This function inspects the decoding results of multi ROI decoding
%When run with show_figures=1 it plots results for a couple of ROIs


if exist('handles_choices')==0
    clear all
    close all

    %Which ii_out and ROIs do you want to inspect
    ii_out=1;

    %time window for accuracy
    time_windows=[3.1 4.1];
    time_windows_pre=[-1 0];

    %These are the windows used for latency measurement
    time_window_lat=[-1.5 10];
    pre_time_window=[-7 -1.5];

    %Thresholds defining accuracy zones
    %     acc_thr=[0.3 0.65]; %Done initially
    acc_thr=[0.35 0.65]; %Done initially

    dt_span=15;
    MLalgo=6;

   
    %This dt_lat is used to find the maximum accuracy attained within dt_lat
    %after latency start
    dt_lat=1;
    sustained_dt=0.2;
    lat_fact=1;
    mad_pre_accuracy=0.05;

    %When run with show_figures=1 it plots results for a couple of ROIs
    show_figures=1;

    %The variables below need to be defined if show_figures=1
    ROI1=6; %50 has a nice accuracy curve starting at odor on
    %1 has an accuracy close to zero and has low dFF zero fraction for both S+ and S-
    %22 has accuracy below 0.5 with high dFF zero fraction for both S+ and S-
    %62 has a nice increase 2 sec after odor on
    %6 starts going up with the FV click

    time_to_inspect1=4; %Time to inspect the dFF

     %Load file
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per_rdec.mat'},'Select the .m file with all the choices for analysis');



else
    ii_out=handles_choices.ii_out;
    show_figures=handles_choices.show_figures;
    pre_perFileName=handles_choices.pre_perFileName;
    pre_perPathName=handles_choices.pre_perPathName;
    time_windows=handles_choices.time_windows;
    time_window_lat=handles_choices.time_window_lat;
    pre_time_window=handles_choices.pre_time_window;
    acc_thr=handles_choices.acc_thr;
    dt_span=handles_choices.dt_span;
    MLalgo=handles_choices.MLalgo;
    time_windows_pre=handles_choices.time_windows_pre;

    ROI1=1;
    time_to_inspect1=4; %Time to inspect the dFF

      %This dt_lat is used to find the maximum accuracy attained within dt_lat
    %after latency start
    dt_lat=1;
    sustained_dt=0.2;
    lat_fact=1;
    mad_pre_accuracy=0.05;
end



%Definition of other variables
figNo=0;



delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;



load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAnInspectMultiROI run for ' pre_perFileName '\n\n']);



handles_out2=handles_out.ii_out(ii_out).handles_out;
dt=handles_out2.dt;
no_time_points=length(handles_out2.time_span);
no_ROI_draws=handles_out2.no_ROI_draws;
time_span=handles_out2.time_span;
dFF_per_trial_sm=handles_out2.dFF_per_trial_sm;
dFF_per_trial_sp=handles_out2.dFF_per_trial_sp;

if ROI1>no_ROI_draws
    ROI1=1;
end

this_ii_time_span=find(time_span>=time_to_inspect1,1,'first');

%Calculate accuracy in the time window
accuracy_per_ROI=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI=[accuracy_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end
handles_outd.accuracy_per_ROI=accuracy_per_ROI;

accuracy_per_ROI_sh=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_sh=[accuracy_per_ROI_sh mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end
handles_outd.accuracy_per_ROI_sh=accuracy_per_ROI_sh;

%Calculate accuracy in the pre time window
accuracy_per_ROI_pre=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_pre=[accuracy_per_ROI_pre mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows_pre(1))&(time_span<=time_windows_pre(2))),2))];
end
handles_outd.accuracy_per_ROI_pre=accuracy_per_ROI_pre;

accuracy_per_ROI_sh_pre=[];
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_sh_pre=[accuracy_per_ROI_sh_pre mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows_pre(1))&(time_span<=time_windows_pre(2))),2))];
end
handles_outd.accuracy_per_ROI_sh_pre=accuracy_per_ROI_sh;

%Quantify outliers in the odor window

ii_t_start=find(time_span>=time_windows(1),1,'first');
ii_t_end=find(time_span<=time_windows(2),1,'last');


dFF_sp_outlier_frac_per_ROI=[];
dFF_sm_outlier_frac_per_ROI=[];

for iiROI=1:no_ROI_draws
    dFF_sp_outlier_frac_this_ROI=[];
    dFF_sm_outlier_frac_this_ROI=[];
    for ii_t=ii_t_start:ii_t_end
        
        this_dFF_per_trial_sm=zeros(1,size(dFF_per_trial_sm,1));
        this_dFF_per_trial_sm(1,:)=dFF_per_trial_sm(:,iiROI,ii_t);
        dFF_sm_outlier_frac_this_ROI=[dFF_sm_outlier_frac_this_ROI sum(isoutlier(this_dFF_per_trial_sm,"median",ThresholdFactor=3))/length(this_dFF_per_trial_sm)];
        
        this_dFF_per_trial_sp=zeros(1,size(dFF_per_trial_sp,1));
        this_dFF_per_trial_sp(1,:)=dFF_per_trial_sp(:,iiROI,ii_t);
        dFF_sp_outlier_frac_this_ROI=[dFF_sp_outlier_frac_this_ROI sum(isoutlier(this_dFF_per_trial_sp,"median",ThresholdFactor=3))/length(this_dFF_per_trial_sp)];
         
        pfft=1;
    end
    dFF_sp_outlier_frac_per_ROI(iiROI)=mean(dFF_sp_outlier_frac_this_ROI);
    dFF_sm_outlier_frac_per_ROI(iiROI)=mean(dFF_sm_outlier_frac_this_ROI);
end


if show_figures==1
    %Accuracy histogram  odor window
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    edges=[0:0.05:1];

    hold on

    histogram(accuracy_per_ROI,edges,'Normalization','Probability')
    histogram(accuracy_per_ROI_sh,edges,'Normalization','Probability')

    title(['Acccuracy histogram, odor window'])
    xlabel('Accuracy')

    fprintf(1,'Mean accuracy odor all runs= %d\n',mean(accuracy_per_ROI))
    fprintf(1,'Shuffled trial mean accuracy odor all runs= %d\n',mean(accuracy_per_ROI_sh))

     %Accuracy histogram  pre window
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    edges=[0:0.05:1];

    hold on

    histogram(accuracy_per_ROI_pre,edges,'Normalization','Probability')
    histogram(accuracy_per_ROI_sh_pre,edges,'Normalization','Probability')

    title(['Acccuracy histogram, pre window'])
    xlabel('Accuracy')

    fprintf(1,'Mean accuracy pre all runs= %d\n',mean(accuracy_per_ROI_pre))
    fprintf(1,'Shuffled trial mean accuracy pre all runs= %d\n',mean(accuracy_per_ROI_sh_pre))

    %Plot pseudocolor accuracy vs dFF outlier fraction 
    %This is odor window
    rand_shift=0.5*rand(1,no_ROI_draws);
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on
    colors=colormap(jet);

    for iiROI=1:no_ROI_draws
        ii_color=1+ceil(accuracy_per_ROI(iiROI)*255);
        plot(dFF_sm_outlier_frac_per_ROI(iiROI)+rand_shift(iiROI),dFF_sp_outlier_frac_per_ROI(iiROI)+rand_shift(iiROI),'o','MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
            'MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
            'MarkerSize',5)
    end

    title('Pseudocolor accuracy')
    xlabel('Fraction of outliers for dFF S-')
    ylabel('Fraction of outliers for dFF S+')

    %Plot pseudocolor accuracy vs dFF fraction positive
    %This is odor window
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on
    colors=colormap(jet);

    ii_t_start=find(time_span>=time_windows(1),1,'first');
    ii_t_end=find(time_span<=time_windows(2),1,'last');


    dFF_sp_nonzero_f_ROI1=[];
    dFF_sm_nonzero_f_ROI1=[];

    for iiROI=1:no_ROI_draws

        for ii_t=ii_t_start:ii_t_end
            dFF_sm_nonzero=0;
            n_dFF_sm=0;
            for trNo=1:size(dFF_per_trial_sm,1)
                n_dFF_sm=n_dFF_sm+1;
                if dFF_per_trial_sm(trNo,iiROI,ii_t)>0.01
                    dFF_sm_nonzero=dFF_sm_nonzero+1;
                end
            end
            dFF_sm_nonzero_f=dFF_sm_nonzero/n_dFF_sm;

            dFF_sp_nonzero=0;
            n_dFF_sp=0;
            for trNo=1:size(dFF_per_trial_sp,1)
                n_dFF_sp=n_dFF_sp+1;
                if dFF_per_trial_sp(trNo,iiROI,ii_t)>0.01
                    dFF_sp_nonzero=dFF_sp_nonzero+1;
                end
            end
            dFF_sp_nonzero_f=dFF_sp_nonzero/n_dFF_sp;
            if ((accuracy_per_ROI(iiROI)<acc_thr(1)||(accuracy_per_ROI(iiROI)>acc_thr(2))))
                ii_color=1+ceil(accuracy_per_ROI(iiROI)*255);
                plot(dFF_sm_nonzero_f,dFF_sp_nonzero_f,'o','MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                    'MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                    'MarkerSize',5)
            end
            if (iiROI==ROI1)&(ii_t==this_ii_time_span)
                dFF_sp_nonzero_f_ROI1=dFF_sp_nonzero_f;
                dFF_sm_nonzero_f_ROI1=dFF_sm_nonzero_f;
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

end

%Calculate prediction per ROI
%S+
prediction_sp_per_ROI=[];
for iiROI=1:no_ROI_draws
    prediction_sp_per_ROI=[prediction_sp_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end
handles_outd.prediction_sp_per_ROI=prediction_sp_per_ROI;

%S-
prediction_sm_per_ROI=[];
for iiROI=1:no_ROI_draws
    prediction_sm_per_ROI=[prediction_sm_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2))];
end
handles_outd.prediction_sm_per_ROI=prediction_sm_per_ROI;

%Plot prediction
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    edges=[0:0.05:1];

   %S+
    histogram(prediction_sp_per_ROI, edges)

    %S-


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
    ax=gca;ax.LineWidth=3;

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
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    this_correct_predict_sh=handles_out2.ROI(ROI1).MLalgo(MLalgo).this_correct_predict_sh;
    this_correct_predict=handles_out2.ROI(ROI1).MLalgo(MLalgo).this_correct_predict;

    CIsm = bootci(1000, @mean, this_correct_predict_sh);
    meansm=mean(this_correct_predict_sh,1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'r');


    CIsm = bootci(1000, @mean, this_correct_predict);
    meansm=mean(this_correct_predict,1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'k');

    plot(time_span',mean(this_correct_predict_sh,1)','-r','LineWidth',1.5)
    plot(time_span',mean(this_correct_predict,1)', '-k','LineWidth',1.5);

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

    %Plot the dFF for ROI1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    this_time_span=time_span;

    this_dFF_per_trial_sm=zeros(size(dFF_per_trial_sm,1),size(dFF_per_trial_sm,3));
    this_dFF_per_trial_sm(:,:)=dFF_per_trial_sm(:,ROI1,:);

    this_dFF_per_trial_sp=zeros(size(dFF_per_trial_sp,1),size(dFF_per_trial_sp,3));
    this_dFF_per_trial_sp(:,:)=dFF_per_trial_sp(:,ROI1,:);


    CIsm = bootci(1000, @mean, this_dFF_per_trial_sm);
    meansm=mean(this_dFF_per_trial_sm,1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(this_time_span',mean(this_dFF_per_trial_sm,1)', CIsm',  'cmap',[238/255 111/255 179/255]);


    CIsm = bootci(1000, @mean, this_dFF_per_trial_sp);
    meansm=mean(this_dFF_per_trial_sp,1);
    CIsm(1,:)=meansm-CIsm(1,:);
    CIsm(2,:)=CIsm(2,:)-meansm;

    [hlsm, hpsm] = boundedline(this_time_span',mean(this_dFF_per_trial_sp,1)', CIsm', 'cmap',[80/255 194/255 255/255]);

    plot(this_time_span',mean(this_dFF_per_trial_sm,1)','Color',[238/255 111/255 179/255],'LineWidth',1.5)
    plot(this_time_span',mean(this_dFF_per_trial_sp,1)', 'Color',[80/255 194/255 255/255],'LineWidth',1.5);

    xlim([-7 15])
    
    this_ylim=ylim;
    ylim([-0.2*(this_ylim(2)-this_ylim(1)) this_ylim(2)])
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

    title(['dFF for ii_ROI '  num2str(ROI1)])
    xlabel('Time(sec)')
    ylabel('dFF')

    %Now understand why the result gives a speciific accuracy
    no_sp_trials=size(handles_out2.dFF_per_trial_sp,1);
    dFF_for_sp_trials=zeros(1,no_sp_trials);
    dFF_for_sp_trials=handles_out2.dFF_per_trial_sp(:,ROI1,this_ii_time_span);
    no_sm_trials=size(handles_out2.dFF_per_trial_sm,1);
    dFF_for_sm_trials=zeros(1,no_sm_trials);
    dFF_for_sm_trials=handles_out2.dFF_per_trial_sm(:,ROI1,this_ii_time_span);

    max_dFF=max([dFF_for_sp_trials; dFF_for_sm_trials]);
    min_dFF=min([dFF_for_sp_trials; dFF_for_sm_trials]);
    edges=[min_dFF-0.05*(max_dFF-min_dFF):0.05*1.2*(max_dFF-min_dFF):max_dFF]+0.05*(max_dFF-min_dFF);
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    hold on

    histogram(dFF_for_sp_trials,edges,'Normalization','Probability','FaceColor',[80/255 194/255 255/255],'EdgeColor',[0 0 0])
    histogram(dFF_for_sm_trials,edges,'Normalization','Probability','FaceColor',[238/255 111/255 179/255],'EdgeColor',[0 0 0])

    title(['Histogram for dFF for ii_ROI '  num2str(ROI1)])
    xlabel('dFF')
    ylabel('Probability')

    pffft=1;

%     %Now plot the accuracy for ROI2
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
%     ax=gca;ax.LineWidth=3;
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
% 
%     hold on
% 
%     this_correct_predict_sh=handles_out2.ROI(ROI2).MLalgo(MLalgo).this_correct_predict_sh;
%     this_correct_predict=handles_out2.ROI(ROI2).MLalgo(MLalgo).this_correct_predict;
% 
%     CIsm = bootci(1000, @mean, this_correct_predict_sh);
%     meansm=mean(this_correct_predict_sh,1);
%     CIsm(1,:)=meansm-CIsm(1,:);
%     CIsm(2,:)=CIsm(2,:)-meansm;
% 
%     [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'r');
% 
% 
%     CIsm = bootci(1000, @mean, this_correct_predict);
%     meansm=mean(this_correct_predict,1);
%     CIsm(1,:)=meansm-CIsm(1,:);
%     CIsm(2,:)=CIsm(2,:)-meansm;
% 
%     [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'k');
% 
%     plot(time_span',mean(this_correct_predict_sh,1)','-r','LineWidth',1.5)
%     plot(time_span',mean(this_correct_predict,1)', '-k','LineWidth',1.5);
% 
% %     text(30,0.75,'Shuffled','Color','r')
% %     text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])
% 
% 
%     ylim([-0.1 1.1])
%     xlim([-7 15])
%     this_ylim=ylim;
% 
%     %FV
%     rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
%     plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])
% 
%     %Odor on markers
%     plot([0 0],this_ylim,'-k')
%     odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')
% 
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')
% 
%     title(['Accuracy for ii_ROI '  num2str(ROI2)])
%     xlabel('Time(sec)')
%     ylabel('Accuracy')

    %Calculate accuracy in a wider post-stimulus window and latency
end

%Calculate latency
accuracy_per_ROIw=[];
latency_per_ROI=[];
cropped_time_span_ii=find((time_span>=time_window_lat(1))&(time_span<=time_window_lat(2)));
cropped_time_span=time_span(cropped_time_span_ii);

%Calculate the mean absolute deviation (as opposed to the STD)
all_pre_accuracy_values=[];
for iiROI=1:no_ROI_draws
    all_pre_accuracy_values=[all_pre_accuracy_values mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=pre_time_window(1))&(time_span<=pre_time_window(2))),1)];
end

% mad_pre_accuracy=mad(all_pre_accuracy_values);

% sustained_dt=0.2;
sustained_ii=ceil(sustained_dt/(time_span(2)-time_span(1)));
% lat_fact=1;

for iiROI=1:no_ROI_draws
    accuracy_per_ROIw=[accuracy_per_ROIw prctile(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_window_lat(1))&(time_span<=time_window_lat(2))),1),95)];

    this_acc_timecourse=mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_window_lat(1))&(time_span<=time_window_lat(2))),1);

    not_found=1;
    ii_shift=1;
    while not_found==1
        this_latency_delta_ii=find(this_acc_timecourse(ii_shift:end)>0.5+lat_fact*mad_pre_accuracy,1,'first');
        if (~isempty(this_latency_delta_ii))&(ii_shift+this_latency_delta_ii+sustained_ii<=length(this_acc_timecourse))
            if sum(this_acc_timecourse(ii_shift+this_latency_delta_ii-1:ii_shift+this_latency_delta_ii+sustained_ii)<0.5+lat_fact*mad_pre_accuracy)==0
                not_found=0;
                latency_per_ROI(iiROI)=cropped_time_span(ii_shift+this_latency_delta_ii-1);
            else
                ii_shift=ii_shift+this_latency_delta_ii;
                if ii_shift>length(this_acc_timecourse)
                    not_found=0;
                    latency_per_ROI(iiROI)=NaN;
                end
            end
        else
            not_found=0;
            latency_per_ROI(iiROI)=NaN;
        end
    end

end

handles_outd.latency_per_ROI=latency_per_ROI;
handles_outd.accuracy_per_ROIw=accuracy_per_ROIw;
handles_outd.mad_pre_accuracy=mad_pre_accuracy;

if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    edges=[0:0.05:1];

    hold on

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
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

    edges=[-1.5:0.5:10.5];
    histogram(latency_per_ROI(accuracy_per_ROIw>acc_thr(2)),edges)

    title(['Histogram for latency for ROIs with accuracy per ROIw > ' num2str(acc_thr(2))])
    xlabel('Time(sec)')
    ylabel('Counts')

    %Plot the accuracy for ROIs with accuracy_per_ROIw>acc_thr(2)
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .6])

    hold on

    ii_shift=0;

    for iiROI=1:no_ROI_draws
        if accuracy_per_ROIw(iiROI)>acc_thr(2)
            ii_shift=ii_shift+1;
            this_correct_predict_sh=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh;
            this_correct_predict=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict;

            plot(time_span',mean(this_correct_predict_sh,1)'+ii_shift,'-r','LineWidth',1.5)
            plot(time_span',mean(this_correct_predict,1)'+ii_shift, '-k','LineWidth',1.5);
            plot([latency_per_ROI(iiROI) latency_per_ROI(iiROI)],[ii_shift ii_shift+1],'-b')
            pffft=1;
        end
    end

    text(30,0.75,'Shuffled','Color','k')
    text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])


    ylim([0.8 ii_shift+1.1])
    xlim([-7 15])
    this_ylim=ylim;

    %FV
%     rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 0],[1 1],'Color',[0.7 0.7 0.7],'LineWidth',5);
   
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])
 
    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 mean(delta_odor)],[1 1],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[1 1],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

    title(['Accuracy for ROIs with accuracy per ROIw > ' num2str(acc_thr(2))])
    xlabel('Time(sec)')
    ylabel('Accuracy')
    

%     %Plot the second derivative of accuracy for ROIs with accuracy_per_ROIw>acc_thr(2)
%     %This did not help
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
% 
%     hFig = figure(figNo);
%     ax=gca;ax.LineWidth=3;
% 
%     set(hFig, 'units','normalized','position',[.05 .1 .3 .6])
% 
%     hold on
% 
%     ii_shift=0;
% 
%     for iiROI=1:no_ROI_draws
%         if accuracy_per_ROIw(iiROI)>acc_thr(2)
%             
% %             this_correct_predict_sh=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh;
%             this_correct_predict=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict;
%             DDthis_correct_predict = movingslope(movingslope(mean(this_correct_predict),7),7);
% 
%             plot(time_span',DDthis_correct_predict'+ii_shift,'-k','LineWidth',1.5)
% 
%             ii_shift=ii_shift+1.2*(max(DDthis_correct_predict)-min(DDthis_correct_predict));
% 
%         end
%     end
% 
%     text(30,0.75,'Shuffled','Color','k')
%     text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])
% 
% 
%     ylim([-0.1 ii_shift+1.1])
%     xlim([-7 15])
%     this_ylim=ylim;
% 
%     %FV
% %     rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
%     plot([-1.5 0],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'Color',[0.7 0.7 0.7],'LineWidth',5);
%    
%     plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])
% 
%     %Odor on markers
%     plot([0 0],this_ylim,'-k')
%     odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
%     plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')
% 
%     %Reinforcement markers
%     plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
%     reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
%     plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')
% 
%     title(['DD accuracy for ROIs with accuracy per ROIw > ' num2str(acc_thr(2))])
%     xlabel('Time(sec)')
%     ylabel('DD Accuracy')
%  
%     pffft=1; 


    

    %Plot the d prime for ROIs with accuracy_per_ROIw>acc_thr(2)
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .3 .6])
    ax=gca;ax.LineWidth=3;
    
    hold on

    ii_shift=0;

    for iiROI=1:no_ROI_draws
        if accuracy_per_ROIw(iiROI)>acc_thr(2)
            ii_shift=ii_shift+3;
            this_dFF_sm_per_trial=zeros(size(handles_out2.dFF_per_trial_sm,1),size(handles_out2.dFF_per_trial_sm,3));
            this_dFF_sm_per_trial(:,:)=handles_out2.dFF_per_trial_sm(:,iiROI,:);
            this_dFF_sp_per_trial=zeros(size(handles_out2.dFF_per_trial_sp,1),size(handles_out2.dFF_per_trial_sp,3));
            this_dFF_sp_per_trial(:,:)=handles_out2.dFF_per_trial_sp(:,iiROI,:);
            this_d_prime=zeros(1,size(handles_out2.dFF_per_trial_sp,3));
            for ii_t=1:size(handles_out2.dFF_per_trial_sp,3)
                these_dFF_sp=zeros(1,size(handles_out2.dFF_per_trial_sp,1));
                these_dFF_sp(1,:)=this_dFF_sp_per_trial(:,ii_t);
                these_dFF_sm=zeros(1,size(handles_out2.dFF_per_trial_sm,1));
                these_dFF_sm(1,:)=this_dFF_sm_per_trial(:,ii_t);

                this_d_prime(1,ii_t)=(mean(these_dFF_sp)-mean(these_dFF_sm))/sqrt((std(these_dFF_sp)^2 + std(these_dFF_sm)^2 )/2);

                if (sum(these_dFF_sm==0)==length(these_dFF_sm))&(sum(these_dFF_sp==0)==length(these_dFF_sp))
                    this_d_prime(1,ii_t)=0;
                end
                %             if isnan(this_d_prime(1,ii_t))
                %                 pffft=1;
                %             end
            end
            plot(this_time_span',mean(this_d_prime,1)'+ii_shift,'-k','LineWidth',1.5)
        end
    end


    xlim([-7 15])
    this_ylim=ylim;

    %FV
%     rectangle(Position=[-1.5,this_ylim(1)+0.05*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
    plot([-1.5 0],[this_ylim(1)+0.03*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.03*(this_ylim(2)-this_ylim(1))],'Color',[0.7 0.7 0.7],'LineWidth',5);
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.03*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.03*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.03*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.03*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

    ylim([0 37])
    title(['d prime for ROIs with accuracy per ROIw > ' num2str(acc_thr(2))])
    xlabel('Time(sec)')
    ylabel('d prime')

end
pfft=1;



