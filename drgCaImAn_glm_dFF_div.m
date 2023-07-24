function [sig_div,p_val]=drgCaImAn_glm_dFF_div(handles_in)
%Performs a glm for dFF with time and spm as the two independent variables

%First time for divergence

time_span=handles_in.time_span;
t_start=handles_in.t_start;
t_end=handles_in.t_end;
dFFsminus=handles_in.dFFsminus;
dFFsplus=handles_in.dFFsplus;


show_figures=0;

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];


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
    set(hFig, 'units','normalized','position',[.2 .2 .4 .4])

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


    plot(time_span',mean(dFFsminus')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
    plot(time_span',mean(dFFsplus')','Color',[0 114/255 178/255],'LineWidth',1.5);

    plot(time_span',mean(dFFsminus')' +CIpvsm(2,:)','-r');

    this_ylim=ylim;
    for ii_te=1:length(timeEvents)
        plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
    end


    xlim([-7 15])

    xlabel('Time(sec)')
    ylabel('dFF')


end

sig_div=0;
p_val=mdl.Coefficients.pValue(2);
if mdl.Coefficients.pValue(2)<=0.05
    sig_div=1;
end

