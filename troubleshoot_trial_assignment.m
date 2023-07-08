%troubleshoot_trial_assignment
%Used to troubleshoot trial assignment in
%drgCaImAn_SVZ_entire_session_shufflingv3
%and drgCaImAn_parse_out_trialsv2

%Place abreak after drgCaImAn_parse_out_trialsv2 in 
%drgCaImAn_SVZ_entire_session_shufflingv3

%Number and identity of trials using events
% 6=Hit (on for the duration of odor on)
% 7=Miss
% 8=FA
% 9=CR
at_end=0;
epochEvents=handles.dropcData.epochEvent;
this_ii=0;
evhit_per_trial=[];
evmiss_per_trial=[];
evcr_per_trial=[];
evfa_per_trial=[];
while at_end==0
    dii_next_trial=find((epochEvents(this_ii+1:end)==6)|(epochEvents(this_ii+1:end)==7)|(epochEvents(this_ii+1:end)==8)|(epochEvents(this_ii+1:end)==9),1,'first');
    if isempty(dii_next_trial)
        at_end=1;
    else
        switch epochEvents(this_ii+dii_next_trial)
            case 6
                % Hit
                evhit_per_trial=[evhit_per_trial 1];
                evmiss_per_trial=[evmiss_per_trial 0];
                evcr_per_trial=[evcr_per_trial 0];
                evfa_per_trial=[evfa_per_trial 0];
            case 7
                % Miss
                evhit_per_trial=[evhit_per_trial 0];
                evmiss_per_trial=[evmiss_per_trial 1];
                evcr_per_trial=[evcr_per_trial 0];
                evfa_per_trial=[evfa_per_trial 0];
            case 8
                % FA
                evhit_per_trial=[evhit_per_trial 0];
                evmiss_per_trial=[evmiss_per_trial 0];
                evcr_per_trial=[evcr_per_trial 0];
                evfa_per_trial=[evfa_per_trial 1];
            case 9
                % CR
                evhit_per_trial=[evhit_per_trial 0];
                evmiss_per_trial=[evmiss_per_trial 0];
                evcr_per_trial=[evcr_per_trial 1];
                evfa_per_trial=[evfa_per_trial 0];
        end
        this_ii=this_ii+dii_next_trial;
    end

end

no_trials=length(evhit_per_trial)
no_splus_trials=sum(evhit_per_trial)+sum(evmiss_per_trial)
no_sminus_trials=sum(evcr_per_trial)+sum(evfa_per_trial)

%Number reported by drgCaImAn_parse_out_trialsv2
no_trials=length(hit_per_trial)
no_splus_trials=sum(hit_per_trial)+sum(miss_per_trial)
no_sminus_trials=sum(cr_per_trial)+sum(fa_per_trial)

%for  Slide1.sld
%-20220131-Frame2-spm_E_10_Iter_3500_output_ncorr_extract_batchv2_pre_per
%drgCaImAn_parse_out_trialsv2 appears to skip the first two trials and the last trial