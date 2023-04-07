%drgCaImAn_analyze_PID
close all
clear all

figNo=0;

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dropc_choices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles_choice.choiceFileName=choiceFileName;
handles_choice.choiceBatchPathName=choiceBatchPathName;

%Load the data and place them on a common time base
t_start=-7;
t_end=15;
dt_display=0.01; %10 msec
time_disp=t_start:dt_display:t_end;
led_traces=zeros(200,length(time_disp));
ii_trace_LED=0;
PID_traces=zeros(200,length(time_disp));
ii_trace_PID=0;
all_delta_odor=[];
all_delta_odor_on_reinf_on=[];
all_delta_reinf=[];


for fileNo=handles_choice.first_file:handles_choice.no_files
    
    fnameca=handles_choice.csvFileName{fileNo};
    pre_per_name=[handles_choice.PathNamecsv{fileNo} fnameca(1:end-4) '_batchv2_pre_per.mat'];
    load(pre_per_name)
    all_delta_odor=[all_delta_odor delta_odor];
    all_delta_odor_on_reinf_on=[all_delta_odor_on_reinf_on delta_odor_on_reinf_on];
    all_delta_reinf=[all_delta_reinf delta_reinf];

    for ii_trace=1:size(splus_traces,1)
        %Make sure there is an LED signal and that this is not just the
        %signal going up
        if sum(splus_traces(ii_trace,:)>200)>20
            ii_trace_LED=ii_trace_LED+1;
            for ii_t=1:length(time_disp)
                jj_before=find(time_to_eventSp<time_disp(ii_t),1,'last');
                jj_after=find(time_to_eventSp>time_disp(ii_t),1,'first');
                led_traces(ii_trace_LED,ii_t)=splus_traces(ii_trace,jj_before)+(time_disp(ii_t)-time_to_eventSp(jj_before))*...
                    (splus_traces(ii_trace,jj_after)-splus_traces(ii_trace,jj_before))/(time_to_eventSp(jj_after)-time_to_eventSp(jj_before));
            end
        end
    end

     for ii_trace=1:size(splus_PID_traces,1)
        %Make sure there is an LED signal and that this is not just the
        %signal going up
        if sum(splus_PID_traces(ii_trace,:)>0.5)>20
            ii_trace_PID=ii_trace_PID+1;
            for ii_t=1:length(time_disp)
                jj_before=find(time_to_event_PIDsp<time_disp(ii_t),1,'last');
                jj_after=find(time_to_event_PIDsp>time_disp(ii_t),1,'first');
                PID_traces(ii_trace_PID,ii_t)=mean(splus_PID_traces(ii_trace,...
                    (time_to_event_PIDsp>time_disp(ii_t)-dt_display)&(time_to_event_PIDsp<time_disp(ii_t)+dt_display)));
            end
        end
     end

end

led_traces=led_traces(1:ii_trace_LED,:);
PID_traces=PID_traces(1:ii_trace_PID,:);


pffft=1;


%Plot the LED timecourse
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

ax=gca;

set(hFig, 'units','normalized','position',[.2 .2 .4 .4])
hold on

try
    CIeu = bootci(1000, @mean, led_traces);
    meaneu=mean(led_traces,1);
    CIeu(1,:)=meaneu-CIeu(1,:);
    CIeu(2,:)=CIeu(2,:)-meaneu;


    [hlsp, hpsp] = boundedline(time_disp',meaneu', CIeu', 'cmap',[0/255 0/255 0/255]);
catch
    plot(time_disp',meaneu', '-', 'Color',[0/255 0/255 0/255]);
end


this_ylim=ylim;

%Place epoch markers

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor
rectangle(Position=[0,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(all_delta_odor),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0 0 0], EdgeColor=[0 0 0])
plot([0 0],[this_ylim],'-k')
plot([mean(all_delta_odor) mean(all_delta_odor)],[this_ylim],'-k')

%Reinforcement
rectangle(Position=[mean(all_delta_odor_on_reinf_on),this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(all_delta_reinf),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[1 0 0], EdgeColor=[1 0 0])
plot([mean(all_delta_odor_on_reinf_on)+mean(all_delta_reinf) mean(all_delta_odor_on_reinf_on)+mean(all_delta_reinf)],[this_ylim],'-r')
plot([mean(all_delta_odor_on_reinf_on) mean(all_delta_odor_on_reinf_on)],[this_ylim],'-r')


% legend('Within S+','Within S-', 'Between')

plot([time_disp(1) time_disp(end)],[log10(0.05) log10(0.05)],'LineWidth',2)

xlim([-7 15])
ax.LineWidth=3;
title('LED time course ')
xlabel('Time(sec)')
ylabel('AU')

%Find the onset of the LED signal
LED_onset_times=[];
for ii_t=1:ii_trace_LED
    this_trace=zeros(1,length(time_disp));
    this_trace(1,:)=led_traces(ii_t,:);
    this_trace=this_trace-mean(this_trace(1:200));
    this_SD=std(this_trace(1:200));
    onset_ii=find(this_trace>mean(this_trace(1:200))+10*this_SD,1,'first');
    LED_onset_times(ii_t)=time_disp(onset_ii);
end

 fprintf(1, '\nMean onset time for the LED signal = %d\n',mean(LED_onset_times));


%Plot the PID timecourse
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

ax=gca;

set(hFig, 'units','normalized','position',[.2 .2 .4 .4])
hold on

try
    CIeu = bootci(1000, @mean, PID_traces);
    meaneu=mean(PID_traces,1);
    CIeu(1,:)=meaneu-CIeu(1,:);
    CIeu(2,:)=CIeu(2,:)-meaneu;


    [hlsp, hpsp] = boundedline(time_disp',meaneu', CIeu', 'cmap',[0/255 0/255 0/255]);
catch
    plot(time_disp',meaneu', '-', 'Color',[0/255 0/255 0/255]);
end


this_ylim=ylim;

%Place epoch markers

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor
rectangle(Position=[0,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(all_delta_odor),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0 0 0], EdgeColor=[0 0 0])
plot([0 0],[this_ylim],'-k')
plot([mean(all_delta_odor) mean(all_delta_odor)],[this_ylim],'-k')

%Reinforcement
rectangle(Position=[mean(all_delta_odor_on_reinf_on),this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(all_delta_reinf),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[1 0 0], EdgeColor=[1 0 0])
plot([mean(all_delta_odor_on_reinf_on)+mean(all_delta_reinf) mean(all_delta_odor_on_reinf_on)+mean(all_delta_reinf)],[this_ylim],'-r')
plot([mean(all_delta_odor_on_reinf_on) mean(all_delta_odor_on_reinf_on)],[this_ylim],'-r')


% legend('Within S+','Within S-', 'Between')

plot([time_disp(1) time_disp(end)],[log10(0.05) log10(0.05)],'LineWidth',2)

xlim([-7 15])
ax.LineWidth=3;
title('PID time course ')
xlabel('Time(sec)')
ylabel('AU')


%Find the onset of the LED signal
PID_onset_times=[];
for ii_t=1:ii_trace_LED
    this_trace=zeros(1,length(time_disp));
    this_trace(1,:)=PID_traces(ii_t,:);
    this_trace=this_trace-mean(this_trace(1:200));
    this_SD=0.002; %Note: SD is zero!
    onset_ii=find(this_trace>mean(this_trace(1:200))+10*this_SD,1,'first');
    PID_onset_times(ii_t)=time_disp(onset_ii);
end

 fprintf(1, '\nMean onset time for the PID signal = %d\n',mean(PID_onset_times));

%Plot each PID timecourse separately
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

ax=gca;

set(hFig, 'units','normalized','position',[.2 .2 .4 .4])
hold on



for ii_trial=1:size(PID_traces,1)
plot(time_disp',PID_traces(ii_trial,:)', '-', 'Color',[150/255 150/255 150/255],'LineWidth',1);
end

plot(time_disp',mean(PID_traces,1)', '-', 'Color',[0/255 0/255 0/255],'LineWidth',3);

this_ylim=ylim;

%Place epoch markers

%FV
rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

%Odor
rectangle(Position=[0,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(all_delta_odor),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0 0 0], EdgeColor=[0 0 0])
plot([0 0],[this_ylim],'-k')
plot([mean(all_delta_odor) mean(all_delta_odor)],[this_ylim],'-k')

%Reinforcement
rectangle(Position=[mean(all_delta_odor_on_reinf_on),this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),mean(all_delta_reinf),0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[1 0 0], EdgeColor=[1 0 0])
plot([mean(all_delta_odor_on_reinf_on)+mean(all_delta_reinf) mean(all_delta_odor_on_reinf_on)+mean(all_delta_reinf)],[this_ylim],'-r')
plot([mean(all_delta_odor_on_reinf_on) mean(all_delta_odor_on_reinf_on)],[this_ylim],'-r')


% legend('Within S+','Within S-', 'Between')

plot([time_disp(1) time_disp(end)],[log10(0.05) log10(0.05)],'LineWidth',2)

xlim([-7 15])
ax.LineWidth=3;
title('PID time course ')
xlabel('Time(sec)')
ylabel('AU')


pffft=1;
    