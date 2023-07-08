%drgCaImAnSoundDecoding does decoding for sound recording
close all
clear all

rng('shuffle')
training_window=[-0.45 -0.35];
prediction_window=[-4 8];

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoiceSound*.m'},'Select the .m file with all the choices for sound decoding');


fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_decode_entire_session_multi_ROI_fsdz for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;
ii_splus=0;
ii_sminus=0;
ii_all=0;
dB=[];
figNo=0;

show_max=1;
FV_click_ii=[];
FV_click_dB=[];

%Note that I (Diiego) decided to shift the first click to zero
load('/Users/restrepd/Documents/Projects/SFTP/SoundTestFinal/FV_click_ii.mat')

for fileNo=1:length(handles.FileName)
    load([handles.PathName{fileNo} handles.FileName{fileNo}])

    %Downsample the time and frequency
    if fileNo==1
        t_down=[0:0.05:22];
        dt=t_down(2)-t_down(1);
        f_down=y(1):(y(end)-y(1))/99:y(end);
        df=f_down(2)-f_down(1);
        sum_this_dBSp=zeros(length(f_down),length(t_down));
        sum_this_dBSm=zeros(length(f_down),length(t_down));
    end

    this_dB=zeros(length(f_down),length(t_down));

    for t_ii=1:length(t_down)
        for f_ii=1:length(f_down)
            these_z=zeros(sum((y>=f_down(f_ii)-(df/2))&(y<f_down(f_ii)+(df/2))),sum((x>=t_down(t_ii)-(dt/2))&(x<t_down(t_ii)+(dt/2))));
            these_z(:,:)=z((y>=f_down(f_ii)-(df/2))&(y<f_down(f_ii)+(df/2)),(x>=t_down(t_ii)-(dt/2))&(x<t_down(t_ii)+(dt/2)));
            these_dB=zeros(sum((y>=f_down(f_ii)-(df/2))&(y<f_down(f_ii)+(df/2))),sum((x>=t_down(t_ii)-(dt/2))&(x<t_down(t_ii)+(dt/2))));
            these_dB(:,:)=10*log10(real(these_z).^2+imag(these_z).^2);
            this_dB(f_ii,t_ii)=mean(these_dB(:));
        end
    end

    if handles.odor_stimulus(fileNo)==1
        ii_splus=ii_splus+1;
        dB.splus(ii_splus).this_dB=this_dB;
        dB.splus(ii_splus).fileNo=fileNo;
        sum_this_dBSp=sum_this_dBSp+this_dB;
    else
        ii_sminus=ii_sminus+1;
        dB.sminus(ii_sminus).this_dB=this_dB;
        dB.sminus(ii_sminus).fileNo=fileNo;
        sum_this_dBSm=sum_this_dBSm+this_dB;
    end

    ii_all=ii_all+1;
    ii_from=find(t_down>=-4+7,1,'first');
    data=mean(this_dB(f_down>=10000,ii_from:end));

%     [maxFV,ii_max]=max(data);
%     if ~isempty(find(data>(mean(data)+0.5*(maxFV-mean(data))),1,'first'))
%         FV_click_ii(ii_all)=find(data>(mean(data)+0.5*(maxFV-mean(data))),1,'first');
%     else
%         FV_click_ii(ii_all)=ii_max;
%     end
%   
% 
%     FV_click_ii(ii_all)=FV_click_ii(ii_all)+ii_from+1;




    if show_max==1
        %Show the collapsed figure
        figNo=2;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);
        ax=gca;ax.LineWidth=3;

        set(hFig, 'units','normalized','position',[.05 .55 .7 .3])

        hold on
        plot(data,'-b')
        plot([FV_click_ii(ii_all)-ii_from FV_click_ii(ii_all)-ii_from],[min(data) max(data)],'-k')
    end

    %Show the pseudocolor figure
    figNo=1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);
    ax=gca;ax.LineWidth=3;

    set(hFig, 'units','normalized','position',[.05 .1 .7 .3])

    hold on

    time=repmat(t_down,length(f_down),1)-7;
    frequency=repmat(f_down', 1, length(t_down));


    drg_pcolor(time,frequency,this_dB)


    colormap fire
    shading interp
    %     caxis([mindFF maxdFF]);
    xlabel('Time (sec)')
    ylabel('Frequency (kHz)');
    if handles.odor_stimulus(fileNo)==1
        title(['Decibels for file No ' num2str(fileNo) ' S+'])
        fprintf(1,['Processed file No ' num2str(fileNo) ' S+\n'])
    else
        title(['Decibels for file No ' num2str(fileNo) ' S-'])
        fprintf(1,['Processed file No ' num2str(fileNo) ' S-\n'])
    end
    if show_max==1
        plot([t_down(FV_click_ii(ii_all))-7 t_down(FV_click_ii(ii_all))-7],[3000 10000],'-k','LineWidth',2)
        FV_click_dB=[FV_click_dB this_dB(:,FV_click_ii(ii_all))];
    end
    xlim([-4 8])
    ylim([f_down(1) f_down(end)])
    pffft=1; 
end

max_dB=max([max(sum_this_dBSp/ii_splus) max(sum_this_dBSm/ii_sminus)]);
min_dB=min([min(sum_this_dBSp/ii_splus) min(sum_this_dBSm/ii_sminus)]);


%Plot average spectrogram for S+

%First shift this_dB to place the click at zero
sum_this_dBSp=zeros(length(f_down),length(t_down));
sum_this_dBSm=zeros(length(f_down),length(t_down));
median_FV_click_ii=median(FV_click_ii);

click_dBSp=zeros(length(f_down),ii_splus);
for ii=1:ii_splus
    this_dB=[];
    this_dB=dB.splus(ii).this_dB;
    click_dBSp(:,ii)=this_dB(:,FV_click_ii(dB.splus(ii).fileNo));
    K=FV_click_ii(dB.splus(ii).fileNo)-median_FV_click_ii+9;
    this_dB_shifted = circshift(this_dB', K )';
    sum_this_dBSp=sum_this_dBSp+this_dB_shifted;
end

click_dBSm=zeros(length(f_down),ii_sminus);
for ii=1:ii_sminus
    this_dB=[];
    this_dB=dB.sminus(ii).this_dB;
    click_dBSm(:,ii)=this_dB(:,FV_click_ii(dB.sminus(ii).fileNo));
    K=FV_click_ii(dB.sminus(ii).fileNo)-median_FV_click_ii+9;
    this_dB_shifted = circshift(this_dB', K )';
    sum_this_dBSm=sum_this_dBSm+this_dB_shifted;
end


min_dB=-35;
max_dB=-20;
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.05 .1 .7 .3])

time=repmat(t_down,length(f_down),1)-7;
frequency=repmat(f_down', 1, length(t_down));


drg_pcolor(time,frequency,sum_this_dBSp/ii_splus)


colormap fire
shading interp
caxis([min_dB max_dB]);
xlim([-4 8])
xlabel('Time (sec)')
ylabel('Frequency (kHz)');

title(['Decibels for average S+'])


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.05 .1 .7 .3])

time=repmat(t_down,length(f_down),1)-7;
frequency=repmat(f_down', 1, length(t_down));


drg_pcolor(time,frequency,sum_this_dBSm/ii_sminus)


colormap fire
shading interp
caxis([min_dB max_dB]);
xlim([-4 8])
xlabel('Time (sec)')
ylabel('Frequency (kHz)');

title(['Decibels for average S-'])

%Plot the rainbow
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.83 .1 .05 .3])

prain=[0:(max_dB-min_dB)/99:max_dB-min_dB];
drg_pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
colormap fire
shading interp
ax=gca;
set(ax,'XTickLabel','')
ylabel('dB')

%Now decode S+ and S=
%leave one trial out
%Store the training data in a table.
no_trials=ii_splus+ii_sminus;

ii_start=find(t_down-7>=training_window(1),1,'first');
ii_end=find(t_down-7<=training_window(2),1,'last');
ii_predict_start=find(t_down-7>=prediction_window(1),1,'first');
ii_predict_end=find(t_down-7<=prediction_window(2),1,'last');
predictions=[];

for kk=1:no_trials

    jj=0;
    these_training_measurements=[];
    these_training_decisions=[];

    %Leave one out
    for kk_other=1:no_trials
        if kk_other~=kk
            %I will do S+ first and S- last
            if kk_other<=ii_splus
                is_splus=1;
                ii_sp=kk_other;
            else
                is_splus=0;
                ii_sm=kk_other-ii_splus;
            end

            ii=FV_click_ii(kk_other);
            if is_splus==1
%                 for ii=ii_start:ii_end
                    jj=jj+1;
                    these_training_measurements=[these_training_measurements dB.splus(ii_sp).this_dB(:,ii)];
                    these_training_decisions(jj,1)=1;
%                 end
            else
%                 for ii=ii_start:ii_end
                    jj=jj+1;
                    these_training_measurements=[these_training_measurements dB.sminus(ii_sm).this_dB(:,ii)];
                    these_training_decisions(jj,1)=0;
%                 end
            end
        end
    end
    tblTrn=[];
    tblTrn = array2table(these_training_measurements);

    %Store the decisions in Y
    Y=these_training_decisions;

    Mdl = fitglm(these_training_measurements',Y','Distribution','binomial');

    %Now predict outcome for the entire timecourse for trial kk
    %I will do S+ first and S- last
    if kk<=ii_splus
        is_splus=1;
        ii_sp=kk;
    else
        is_splus=0;
        ii_sm=kk-ii_splus;
    end

    jj=0;
    label_traces=[];
    accuracy=[];
    ii_t=FV_click_ii(kk);
    if is_splus==1
%         for ii_t=ii_predict_start:ii_predict_end
            jj=jj+1;
            this_time_point=zeros(1,length(f_down));
            this_time_point(1,:)=dB.splus(ii_sp).this_dB(:,ii_t);
            [label,score] = predict(Mdl,this_time_point);
            if label>0.5
                label_traces(jj)=1;
                accuracy(jj)=1;
            else
                label_traces(jj)=0;
                accuracy(jj)=0;
            end
%         end
    else
%         for ii_t=ii_predict_start:ii_predict_end
            jj=jj+1;
            this_time_point=zeros(1,length(f_down));
            this_time_point(1,:)=dB.sminus(ii_sm).this_dB(:,ii_t);
            [label,score] = predict(Mdl,this_time_point);
            if label>0.5
                label_traces(jj)=1;
                accuracy(jj)=0;
            else
                label_traces(jj)=0;
                accuracy(jj)=1;
            end
%         end
    end

    predictions.trial(kk).accuracy=accuracy;
    predictions.trial(kk).prediction=label_traces;
    predictions.trial(kk).is_splus=is_splus;

end

%Now do shuffled decoding
for sh_ii=1:5
    for kk=1:no_trials

        jj=0;
        these_training_measurements=[];
        these_training_decisions=[];

        %Leave one out
        for kk_other=1:no_trials
            if kk_other~=kk
                %I will do S+ first and S- last
                if kk_other<=ii_splus
                    is_splus=1;
                    ii_sp=kk_other;
                else
                    is_splus=0;
                    ii_sm=kk_other-ii_splus;
                end


                ii=FV_click_ii(kk_other);

                if is_splus==1
%                     for ii=ii_start:ii_end
                        jj=jj+1;
                        these_training_measurements=[these_training_measurements dB.splus(ii_sp).this_dB(:,ii)];
                        these_training_decisions(jj,1)=1;
%                     end
                else
%                     for ii=ii_start:ii_end
                        jj=jj+1;
                        these_training_measurements=[these_training_measurements dB.sminus(ii_sm).this_dB(:,ii)];
                        these_training_decisions(jj,1)=0;
%                     end
                end
            end
        end
        tblTrn=[];
        tblTrn = array2table(these_training_measurements);

        %Store the decisions in Y
        these_training_decisions=these_training_decisions(randperm(length(these_training_decisions)));
        Y=these_training_decisions;

        Mdl = fitglm(these_training_measurements',Y','Distribution','binomial');

        %Now predict outcome for the entire timecourse for trial kk
        %I will do S+ first and S- last
        if kk<=ii_splus
            is_splus=1;
            ii_sp=kk;
        else
            is_splus=0;
            ii_sm=kk-ii_splus;
        end

        jj=0;
        label_traces=[];
        accuracy=[];
        ii_t=FV_click_ii(kk);
        if is_splus==1
%             for ii_t=ii_predict_start:ii_predict_end
                jj=jj+1;
                this_time_point=zeros(1,length(f_down));
                this_time_point(1,:)=dB.splus(ii_sp).this_dB(:,ii_t);
                [label,score] = predict(Mdl,this_time_point);
                if label>0.5
                    label_traces(jj)=1;
                    accuracy(jj)=1;
                else
                    label_traces(jj)=0;
                    accuracy(jj)=0;
                end
%             end
        else
%             for ii_t=ii_predict_start:ii_predict_end
                jj=jj+1;
                this_time_point=zeros(1,length(f_down));
                this_time_point(1,:)=dB.sminus(ii_sm).this_dB(:,ii_t);
                [label,score] = predict(Mdl,this_time_point);
                if label>0.5
                    label_traces(jj)=1;
                    accuracy(jj)=0;
                else
                    label_traces(jj)=0;
                    accuracy(jj)=1;
                end
%             end
        end

        predictions.trial(kk).sh(sh_ii).accuracy=accuracy;
        predictions.trial(kk).sh(sh_ii).prediction=label_traces;
        predictions.trial(kk).sh(sh_ii).is_splus=is_splus;

    end
end

time_span=t_down(ii_predict_start:ii_predict_end)-7;
this_correct_predict=[];
this_correct_predict_sh=[];
for kk=1:no_trials
    this_correct_predict=[this_correct_predict; predictions.trial(kk).accuracy];
    for sh_ii=1:5
        this_correct_predict_sh=[this_correct_predict_sh; predictions.trial(kk).sh(sh_ii).accuracy];
    end
end

%Plot the accuracy
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on


bar_offset=1;

CIsm = bootci(1000, @mean, this_correct_predict_sh);
meansm=mean(this_correct_predict_sh,1);

bar(bar_offset,meansm,'r')
plot([bar_offset bar_offset],CIsm,'-k','LineWidth',3)



CIsm = bootci(1000, @mean, this_correct_predict);
meansm=mean(this_correct_predict,1);

bar_offset=2;

bar(bar_offset,meansm,'FaceColor',[0.7 .7 .7])
plot([bar_offset bar_offset],CIsm,'-k','LineWidth',3)

xticks([1 2])
xticklabels({'Shuffled','Original'})

title(['Accuracy'])
ylabel('Accuracy')

p = ranksum(this_correct_predict_sh,this_correct_predict);
fprintf(1,['p-value for ranksum between accuracy and shifted accuracy ' num2str(p) '\n'])

%Plot dB at the click
click_dBSp=click_dBSp-min_dB;
click_dBSm=click_dBSm-min_dB;
click_dBSp=click_dBSp';
click_dBSm=click_dBSm';

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

hold on



CIsm = bootci(1000, @mean, click_dBSm);
meansm=mean(click_dBSm,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(f_down,mean(click_dBSm,1)', CIsm',  'cmap',[238/255 111/255 179/255]);


CIsm = bootci(1000, @mean, click_dBSp);
meansm=mean(click_dBSp,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

[hlsm, hpsm] = boundedline(f_down,mean(click_dBSp,1)', CIsm', 'cmap',[80/255 194/255 255/255]);

plot(f_down,mean(click_dBSm,1)','Color',[238/255 111/255 179/255],'LineWidth',1.5)
plot(f_down,mean(click_dBSp,1)', 'Color',[80/255 194/255 255/255],'LineWidth',1.5);

xlim([min(f_down) max(f_down)])

title(['dB at final valve click'])
xlabel('Frequency (Hz)')
ylabel('dB')
pffft=1;

