function [sig_div,p_val,total_trials]=drgCaImAn_glm_dFF_div(handles_in)
%Performs a glm for dFF with time and spm as the two independent variables
 
%First time for divergence

time_span=handles_in.time_span;
t_start=handles_in.t_start;
t_end=handles_in.t_end;
dFFsminus=handles_in.dFFsminus;
dFFsplus=handles_in.dFFsplus;
min_tr=handles_in.min_tr_div;


show_figures=0;

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [-1.5 0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];


glm_div_ii=0;
glm_div=[];

for t=t_start:0.5:t_end

    %dFFsminus
    these_dFFsminus=zeros(1,size(dFFsminus,2));
    these_dFFsminus(1,:)=mean(dFFsminus((time_span>=t)&(time_span<t+1),:));
    glm_div.data(glm_div_ii+1:glm_div_ii+length(these_dFFsminus))=these_dFFsminus;
    glm_div.spm(glm_div_ii+1:glm_div_ii+length(these_dFFsminus))=0*ones(1,length(these_dFFsminus));
    glm_div.time(glm_div_ii+1:glm_div_ii+length(these_dFFsminus))=t*ones(1,length(these_dFFsminus));
    glm_div_ii=glm_div_ii+length(these_dFFsminus);

    %dFFsplus
    these_dFFsplus=zeros(1,size(dFFsplus,2));
    these_dFFsplus(1,:)=mean(dFFsplus((time_span>=t)&(time_span<t+1),:));
    glm_div.data(glm_div_ii+1:glm_div_ii+length(these_dFFsplus))=these_dFFsplus;
    glm_div.spm(glm_div_ii+1:glm_div_ii+length(these_dFFsplus))=1*ones(1,length(these_dFFsplus));
    glm_div.time(glm_div_ii+1:glm_div_ii+length(these_dFFsplus))=t*ones(1,length(these_dFFsplus));
    glm_div_ii=glm_div_ii+length(these_dFFsplus);

end

tbl = table(glm_div.data',glm_div.spm',glm_div.time',...
    'VariableNames',{'dFF','spm','time'});
mdl = fitglm(tbl,'dFF~spm+time+spm*time'...
    ,'CategoricalVars',[2]);




if show_figures==1
    try
        close(1)
    catch
    end
    hFig=figure(1);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.02 .1 .3 .3])

    try
        %S-
        CIpvsm = bootci(1000, {@mean, dFFsminus'},'alpha',0.05);
        meanpvsm=mean(dFFsminus',1);
        CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
        CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;



        %S+
        CIpvsp = bootci(1000, {@mean, dFFsplus'},'alpha',0.05);
        meanpvsp=mean(dFFsplus',1);
        CIpvsp(1,:)=meanpvsp-CIpvsp(1,:);
        CIpvsp(2,:)=CIpvsp(2,:)-meanpvsp;

        hold on

        [hlpvl, hppvl] = boundedline(time_span,mean(dFFsminus'), CIpvsm','cmap',[158/255 31/255 99/255]);
        [hlpvl, hppvl] = boundedline(time_span, mean(dFFsplus'), CIpvsp','cmap',[0 114/255 178/255]);

    catch
    end

    plot(time_span',mean(dFFsminus')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
    plot(time_span',mean(dFFsplus')','Color',[0 114/255 178/255],'LineWidth',1.5);

%     plot(time_span',mean(dFFsminus')' +CIpvsm(2,:)','-r');

        this_ylim=ylim;
    fraction_marker=0.05;
    %FV
    fvhl=plot([-1.5 0],[this_ylim(1)+fraction_marker*(this_ylim(2)-this_ylim(1)) this_ylim(1)+fraction_marker*(this_ylim(2)-this_ylim(1))],'-'...
        ,'Color',[0.9 0.9 0.9],'LineWidth',5);
    plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

    %Odor on markers
    plot([0 0],this_ylim,'-k')
    odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+fraction_marker*(this_ylim(2)-this_ylim(1)) this_ylim(1)+fraction_marker*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+fraction_marker*(this_ylim(2)-this_ylim(1)) this_ylim(1)+fraction_marker*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')



    xlim([-7 15])

    xlabel('Time(sec)')
    ylabel('dFF')


end

sig_div=0;

%Now use the rule in Taxidis: At least 3 or 10% of the trials have spikes
%with the additional rule that the total number of trials should be >20
total_trials=size(dFFsminus,2)+size(dFFsplus,2);

%Estimate number of trials with spikes
spike_trials=0;
for ii=1:size(dFFsminus,2)
    if sum(dFFsminus((time_span>=t_start)&(time_span<=t_end),ii))>0
        spike_trials=spike_trials+1;
    end
end
for ii=1:size(dFFsplus,2)
    if sum(dFFsplus((time_span>=t_start)&(time_span<=t_end),ii))>0
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
    p_val=mdl.Coefficients.pValue(2);
    if mdl.Coefficients.pValue(2)<=0.05
        sig_div=1;
    end
else
    p_val=1;
    sig_div=0;
end

