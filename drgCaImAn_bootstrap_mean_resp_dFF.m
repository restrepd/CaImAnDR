function [sig_resp,n]=drgCaImAn_bootstrap_mean_resp_dFF(handles_in)

%First time for divergence
all_done=0;
time_span=handles_in.time_span;
t_start=handles_in.t_start;
t_end=handles_in.t_end;
min_bins=handles_in.min_bins;
dt_per_bin=handles_in.dt_per_bin;
pre_t_start=handles_in.pre_t_start;
pre_t_end=handles_in.pre_t_end;
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
sig_resp_dt=zeros(1,round((t_end-t_start)/dt_per_bin));
ii=0;
n=size(dFF,2);

these_pre_dFF=[];
for t=pre_t_start:dt_per_bin:pre_t_end
    for trNo=1:n
        these_pre_dFF=[these_pre_dFF mean(dFF((time_span>=t)&(time_span<t+dt_per_bin),trNo))];
    end
end
n_pre=length(these_pre_dFF);

for t=t_start:dt_per_bin:t_end

    these_dFF_bin=zeros(1,n);
    these_dFF_bin(1,:)=mean(dFF((time_span>=t)&(time_span<t+dt_per_bin),:));
    bootstrapStat = zeros(nReps,1);
    for i=1:nReps
        sampX1 = these_pre_dFF(ceil(rand(n_pre,1)*n_pre));
        sampX2 = these_dFF_bin(ceil(rand(n,1)*n));
        bootstrapStat(i) = myMeanDiff(sampX1,sampX2);
    end

    %Calculate the confidence interval (I could make a function out of this...)
    CI = prctile(bootstrapStat,[100*alpha/2,100*(1-alpha/2)]);
    %Hypothesis test: Does the confidence interval cover zero?
    %H=0 the samples do not differ
    %H=1 the sample means do differ
    ii=ii+1;
    sig_resp_dt(ii) = CI(1)>0 | CI(2)<0;
end

if sum(sig_resp_dt==1)>=min_bins
    sig_resp=1;
else
    sig_resp=0;
end


