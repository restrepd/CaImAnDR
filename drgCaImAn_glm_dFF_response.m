function [sig_resp,p_val,total_trials]=drgCaImAn_glm_dFF_response(handles_in)
%Performs a glm for dFF with time and spm as the two independent variables
 
%First time for divergence

time_span=handles_in.time_span;
t_start=handles_in.t_start;
t_end=handles_in.t_end;
pre_t_start=handles_in.pre_t_start;
pre_t_end=handles_in.pre_t_end;
dFF=handles_in.dFF;
min_tr=handles_in.min_tr_resp;


show_figures=0;

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [-1.5 0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];


glm_resp_ii=0;
glm_resp=[];

for t=pre_t_start:0.5:pre_t_end

    %Pre
    these_dFF=zeros(1,size(dFF,2));
    these_dFF(1,:)=mean(dFF((time_span>=t)&(time_span<t+1),:));
    glm_resp.data(glm_resp_ii+1:glm_resp_ii+length(these_dFF))=these_dFF;
    glm_resp.pre_odor(glm_resp_ii+1:glm_resp_ii+length(these_dFF))=0*ones(1,length(these_dFF));
    glm_resp.time(glm_resp_ii+1:glm_resp_ii+length(these_dFF))=t*ones(1,length(these_dFF));
    glm_resp_ii=glm_resp_ii+length(these_dFF);

end

for t=t_start:0.5:t_end
    %Odor
    these_dFF=zeros(1,size(dFF,2));
    these_dFF(1,:)=mean(dFF((time_span>=t)&(time_span<t+1),:));
    glm_resp.data(glm_resp_ii+1:glm_resp_ii+length(these_dFF))=these_dFF;
    glm_resp.pre_odor(glm_resp_ii+1:glm_resp_ii+length(these_dFF))=1*ones(1,length(these_dFF));
    glm_resp.time(glm_resp_ii+1:glm_resp_ii+length(these_dFF))=t*ones(1,length(these_dFF));
    glm_resp_ii=glm_resp_ii+length(these_dFF);

end

tbl = table(glm_resp.data',glm_resp.pre_odor',glm_resp.time',...
    'VariableNames',{'dFF','pre_odor','time'});
mdl = fitglm(tbl,'dFF~pre_odor+time+pre_odor*time'...
    ,'CategoricalVars',[2]);




if show_figures==1
    try
        close(1)
    catch
    end
    hFig=figure(1);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.05 .2 .2 .2])

    %dFF
    CIpvsm = bootci(1000, {@mean, dFF'},'alpha',0.05);
    meanpvsm=mean(dFF',1);
    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;



%     %S+
%     CIpvsp = bootci(1000, {@mean, dFF'},'alpha',0.05);
%     meanpvsp=mean(dFF',1);
%     CIpvsp(1,:)=meanpvsp-CIpvsp(1,:);
%     CIpvsp(2,:)=CIpvsp(2,:)-meanpvsp;

    hold on

    [hlpvl, hppvl] = boundedline(time_span,mean(dFF'), CIpvsm','cmap',[158/255 31/255 99/255]);
%     [hlpvl, hppvl] = boundedline(time_span, mean(dFF'), CIpvsp','cmap',[0 114/255 178/255]);

% 
%     plot(time_span',mean(dFF')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
%     plot(time_span',mean(dFF')','Color',[0 114/255 178/255],'LineWidth',1.5);
% 
%     plot(time_span',mean(dFF')' +CIpvsm(2,:)','-r');

    this_ylim=ylim;
    for ii_te=1:length(timeEvents)
        plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
    end


    xlim([-7 15])

    xlabel('Time(sec)')
    ylabel('dFF')


end

%Now use the rule in Taxidis: At least 3 or 10% of the trials have spikes
%with the additional rule that the total number of trials should be >20
total_trials=size(dFF,2);

%Estimate number of trials with spikes
spike_trials=0;
for ii=1:size(dFF,2)
    if sum(dFF((time_span>=pre_t_start)&(time_span<=t_end),ii))>0
        spike_trials=spike_trials+1;
    end
end


include_neuron=1;
if (total_trials<min_tr)
    include_neuron=0;
end
 
if (spike_trials<0.25*total_trials)
    include_neuron=0;
end

if include_neuron==1
    p_val=min([mdl.Coefficients.pValue(2) mdl.Coefficients.pValue(4)]);
    if mdl.Coefficients.pValue(2)<=0.05
        sig_resp=1;
    else
        sig_resp=0;
    end
else
    p_val=1;
    sig_resp=0;
end

if (p_val<0.05)||(mdl.Coefficients.pValue(4)<=0.05)
    pffft=1;
end

