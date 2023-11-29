function [sig_div,div_t,total_trials]=drgCaImAn_bootstrap_dFF_div(handles_in)

%First time for divergence
all_done=0;
time_span=handles_in.time_span;
t_start=handles_in.t_start;
t_end=handles_in.t_end;
pre_start=handles_in.pre_start;
pre_end=handles_in.pre_end;
dt_required=handles_in.dt_required;
dFFsminus=handles_in.dFFsminus;
dFFsplus=handles_in.dFFsplus;
delta_dFF_for_sig=0.3; %0.3 works well
safety_factor=1; %1.5 with delta_dFF_for_sig=0 too lenient

total_trials=size(dFFsminus,2)+size(dFFsplus,2);

show_figures=0;

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

ii=find(time_span>=t_start,1,'first');

% dt_required=0.5;
delta_ii_required=ceil(dt_required/(time_span(2)-time_span(1)));
ii_included=0;

%S-
CIpvsm = bootci(1000, {@mean, dFFsminus'},'alpha',0.05);
meanpvsm=mean(dFFsminus',1);
CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

pre_dFFsminus=dFFsminus((time_span>=pre_start)&(time_span<=pre_end),:);
SDsm=std(pre_dFFsminus(:));

%S+
CIpvsp = bootci(1000, {@mean, dFFsplus'},'alpha',0.05);
meanpvsp=mean(dFFsplus',1);
CIpvsp(1,:)=meanpvsp-CIpvsp(1,:);
CIpvsp(2,:)=CIpvsp(2,:)-meanpvsp;

pre_dFFsplus=dFFsplus((time_span>=pre_start)&(time_span<=pre_end),:);
SDsp=std(pre_dFFsplus(:));

% delta_dFF_for_sig=sqrt((SDsm^2 +SDsp^2)/2);
 
ii_div=NaN;
div_t=NaN;
sig_div=0;

while all_done==0

    delta_ii_sp_larger_than_sm=find((meanpvsp(ii:end)-CIpvsp(1,ii:end))-delta_dFF_for_sig>(meanpvsm(ii:end)+CIpvsm(2,ii:end)),1,'first');
    delta_ii_sm_larger_than_sp=find((meanpvsm(ii:end)-CIpvsm(1,ii:end))-delta_dFF_for_sig>(meanpvsp(ii:end)+CIpvsp(2,ii:end)),1,'first');

    if (~isempty(delta_ii_sp_larger_than_sm))&(~isempty(delta_ii_sm_larger_than_sp))
        if delta_ii_sp_larger_than_sm<=delta_ii_sm_larger_than_sp
            if ii+delta_ii_sp_larger_than_sm-1+delta_ii_required<length(time_span)
                delta_ii_sp_larger_than_sm_last=find((meanpvsp(ii+delta_ii_sp_larger_than_sm-1:end)...
                    -safety_factor*CIpvsp(1,ii+delta_ii_sp_larger_than_sm-1:end)-delta_dFF_for_sig)>(meanpvsm(ii+delta_ii_sp_larger_than_sm-1:end)...
                    +safety_factor*CIpvsm(2,ii+delta_ii_sp_larger_than_sm-1:end)),1,'last');
                if isempty(delta_ii_sp_larger_than_sm_last)
                    all_done=1;
                    ii_div=ii+delta_ii_sp_larger_than_sm-1;
                    sig_div=1;
                else
                    if delta_ii_sp_larger_than_sm_last>=delta_ii_required
                        all_done=1;
                        ii_div=ii+delta_ii_sp_larger_than_sm-1;
                        sig_div=1;
                    end
                end
                ii=ii+delta_ii_sp_larger_than_sm+delta_ii_sp_larger_than_sm_last;
            else
                all_done=1;
                ii_div=length(time_span);
            end
        else
            if ii+delta_ii_sm_larger_than_sp-1+delta_ii_required<length(time_span)
                delta_ii_sm_larger_than_sp_last=find((meanpvsm(ii+delta_ii_sm_larger_than_sp-1:end)...
                    -safety_factor*CIpvsm(1,ii+delta_ii_sm_larger_than_sp-1:end)-delta_dFF_for_sig)>(meanpvsp(ii+delta_ii_sm_larger_than_sp-1:end)...
                    +safety_factor*CIpvsp(2,ii+delta_ii_sm_larger_than_sp-1:end)),1,'last');
                if isempty(delta_ii_sm_larger_than_sp_last)
                    all_done=1;
                    ii_div=ii+delta_ii_sm_larger_than_sp-1;
                    sig_div=1;
                else
                    if delta_ii_sm_larger_than_sp_last>=delta_ii_required
                        all_done=1;
                        ii_div=ii+delta_ii_sm_larger_than_sp-1;
                        sig_div=1;
                    end
                end
                ii=ii+delta_ii_sm_larger_than_sp+delta_ii_sm_larger_than_sp_last;
            else
                all_done=1;
                ii_div=length(time_span);
            end
        end

    else
        if ~isempty(delta_ii_sp_larger_than_sm)
            if ii+delta_ii_sp_larger_than_sm-1+delta_ii_required<length(time_span)
                delta_ii_sp_larger_than_sm_last=find((meanpvsp(ii+delta_ii_sp_larger_than_sm-1:end)...
                    -safety_factor*CIpvsp(1,ii+delta_ii_sp_larger_than_sm-1:end)-delta_dFF_for_sig)>(meanpvsm(ii+delta_ii_sp_larger_than_sm-1:end)...
                    +safety_factor*CIpvsm(2,ii+delta_ii_sp_larger_than_sm-1:end)),1,'last');
                if isempty(delta_ii_sp_larger_than_sm_last)
                    all_done=1;
                    ii_div=ii+delta_ii_sp_larger_than_sm-1;
                    sig_div=1;
                else
                    if delta_ii_sp_larger_than_sm_last>=delta_ii_required
                        all_done=1;
                        ii_div=ii+delta_ii_sp_larger_than_sm-1;
                        sig_div=1;
                    end
                end
                ii=ii+delta_ii_sp_larger_than_sm+delta_ii_sp_larger_than_sm_last;
            else
                all_done=1;
                ii_div=length(time_span);
            end
        end

        if ~isempty(delta_ii_sm_larger_than_sp)
            if ii+delta_ii_sm_larger_than_sp-1+delta_ii_required<length(time_span)
                delta_ii_sm_larger_than_sp_last=find((meanpvsm(ii+delta_ii_sm_larger_than_sp-1:end)...
                    -safety_factor*CIpvsm(1,ii+delta_ii_sm_larger_than_sp-1:end)-delta_dFF_for_sig)>(meanpvsp(ii+delta_ii_sm_larger_than_sp-1:end)...
                    +safety_factor*CIpvsp(2,ii+delta_ii_sm_larger_than_sp-1:end)),1,'last');
                if isempty(delta_ii_sm_larger_than_sp_last)
                    all_done=1;
                    ii_div=ii+delta_ii_sm_larger_than_sp-1;
                    sig_div=1;
                else
                    if delta_ii_sm_larger_than_sp_last>=delta_ii_required
                        all_done=1;
                        ii_div=ii+delta_ii_sm_larger_than_sp-1;
                        sig_div=1;
                    end
                end
                ii=ii+delta_ii_sm_larger_than_sp+delta_ii_sm_larger_than_sp_last;
            else
                all_done=1;
                ii_div=length(time_span);
            end
        end
    end

    if (isempty(delta_ii_sp_larger_than_sm))&(isempty(delta_ii_sm_larger_than_sp))
        all_done=1;
        ii_div=length(time_span);
    end
end

if sig_div==1
    div_t=time_span(ii_div);

    if show_figures==1
        try
        close(1)
    catch
    end
    hFig=figure(1);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .4 .4])


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

    plot([div_t div_t],this_ylim,'-b','LineWidth',2)

    xlim([-7 15])

    xlabel('Time(sec)')
    ylabel('dFF')

        
    end
 
    pffft=1;
end
