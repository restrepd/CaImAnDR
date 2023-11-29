function [sig_div,div_t]=drgCaImAn_bootstrap_mean_resp_dFF_div(handles_in)

%First time for divergence
all_done=0;
time_span=handles_in.time_span;
t_start=handles_in.t_start;
t_end=handles_in.t_end;
min_bins=handles_in.min_bins;
dt_per_bin=handles_in.dt_per_bin;
pre_start=handles_in.pre_start;
pre_end=handles_in.pre_end;
% dt_required=handles_in.dt_required;
dFF=handles_in.dFF;
% delta_dFF_for_sig=0.3; %0.3 works well
% safety_factor=1; %1.5 with delta_dFF_for_sig=0 too lenient

show_figures=0;

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;


% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

myMeanDiff = @(x1,x2) mean(x1)-mean(x2);
nReps=1000;
alpha=0.05;

%Do a bootstrap for each 0.5 sec segments
sig_div_dt=zeros(1,round((t_end-t_start)/dt_per_bin));
ii=0;
nSp=size(dFFsplus,2);
nSm=size(dFFsminus,2);

for t=t_start:dt_per_bin:t_end
    these_dFFsp=zeros(1,nSp);
    these_dFFsp(1,:)=mean(dFFsplus((time_span>=t)&(time_span<t+dt_per_bin),:));
    these_dFFsm=zeros(1,nSm);
    these_dFFsm(1,:)=mean(dFFsminus((time_span>=t)&(time_span<t+dt_per_bin),:));
    bootstrapStat = zeros(nReps,1);
    for i=1:nReps
        sampX1 = these_dFFsp(ceil(rand(nSp,1)*nSp));
        sampX2 = these_dFFsm(ceil(rand(nSm,1)*nSm));
        bootstrapStat(i) = myMeanDiff(sampX1,sampX2);
    end

    %Calculate the confidence interval (I could make a function out of this...)
    CI = prctile(bootstrapStat,[100*alpha/2,100*(1-alpha/2)]);
    %Hypothesis test: Does the confidence interval cover zero?
    %H=0 the samples do not differ
    %H=1 the sample means do differ
    ii=ii+1;
    sig_div_dt(ii) = CI(1)>0 | CI(2)<0;
end

if sum(sig_div_dt==1)>=min_bins
    sig_div=1;
else
    sig_div=0;
end

ii_div_first=find(sig_div_dt==1,1,"first");
if ~isempty(ii_div_first)
    div_t=(ii_div_first-1)*dt_per_bin+(dt_per_bin/2);
else
    div_t=t_end;
end

