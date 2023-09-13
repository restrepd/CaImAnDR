function [hit_per_trial,cr_per_trial,dFFs_sp_per_trial_per_ROI_post_shifted,...
    dFFs_sm_per_trial_per_ROI_post_shifted,dFF_per_trial_sm,dFF_per_trial_sp,training_decisions_post...
    ,which_model_for_traces_loo,decisions_per_trial,...
    ii_pointer_to_td,epochs_sp_post,measurements_post,...
    ii_post,trial_no,epochs_sm_post,miss_per_trial,fa_per_trial] = ...
    drgCaImAn_parse_out_trials_and_licks(dt, dt_span,epochs,no_points_post_shift,no_points_post,traces,ii_p_threshold,no_odor_trials,...
            time,trimmed_licks)
%This function parses out the trials from the epochs vector

%Note that this line is here because some experiments were run with slow
%6Hz acquisition rate
if no_points_post==0
    no_points_post=1;
end

%training_decisions is 1 if S+ and 2 if S-
%epochs has masks for the following epochs
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR

%Do both S+ and S-
at_end=0;
this_ii=0;
ii_post=0;
ii_pre=0;
ii=0;
trial_no=0;
ii_sp_post=0;
ii_sm_post=0;
ii_which_model=0;
dt_post_which_model=floor(20/dt); %Points that model will be used beyond the training period
ii_span=ceil(dt_span/dt);
dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI_post_shifted=[];
dFFs_sm_per_trial_per_ROI_post_shifted=[];
hit_per_trial=[];
miss_per_trial=[];
cr_per_trial=[];
fa_per_trial=[];
training_decisions_post=[];
decisions_per_trial=[];
which_model_for_traces_loo=no_odor_trials*ones(1,size(traces,2));
measurements_post=[];
ii_pointer_to_td=[];
% measurements_pre=[];


%training_decisions is 1 if S+ and 2 if S-
%epochs has masks for the following epochs
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR
while (at_end==0)
    %6 is hit, 7 is miss, 8 FA and 9 CR
    next_ii_sp=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    delta_ii_sp=find(epochs(this_ii+next_ii_sp:end)~=epochs(this_ii+next_ii_sp),1,'first');

    next_ii_sm=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    delta_ii_sm=find(epochs(this_ii+next_ii_sm:end)~=epochs(this_ii+next_ii_sm),1,'first');

    if (isempty(next_ii_sp))&(isempty(next_ii_sm))
        at_end=1;
    else

        if isempty(next_ii_sm)
            %This is S+

            next_ii=next_ii_sp;
            if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii>0)...
                    &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)&...
                    (this_ii+next_ii+ii_span<length(time))&(this_ii+next_ii-ii_span>=1)&(no_points_post_shift+this_ii+next_ii+ii_p_threshold<=size(traces,2))
                if epochs(this_ii+1+next_ii_sp)==6
                    hit_per_trial=[hit_per_trial 1];
                    miss_per_trial=[miss_per_trial 0];
                else
                    hit_per_trial=[hit_per_trial 0];
                    miss_per_trial=[miss_per_trial 1];
                end
                cr_per_trial=[cr_per_trial 0];
                fa_per_trial=[fa_per_trial 0];
                measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                training_decisions_post(ii_post+1:ii_post+no_points_post)=trimmed_licks(no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1);
                ii_sp_post=ii_sp_post+1;
                dFFs_sp_per_trial_per_ROI_post_shifted(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                dFFs_this_trial=traces(:,this_ii+next_ii-ii_span:this_ii+next_ii+ii_span);
                dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;
                trial_no=trial_no+1;
                decisions_per_trial(trial_no)=1;
                which_model_for_traces_loo(1,ii_which_model+1:this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                ii_which_model=this_ii+next_ii+no_points_post-1+dt_post_which_model;
                ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=this_ii+next_ii:this_ii+next_ii+no_points_post-1;
                epochs_sp_post(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
                %                 measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)';
                %                 epochs_sp_pre(1,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)=1;
                this_ii=this_ii+next_ii+delta_ii_sp-1;
                ii_post=ii_post+no_points_post;
                %                 ii_pre=ii_pre+no_points_pre;
            else
                if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||((no_points_post_shift+this_ii+next_ii-ii_span>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs)))||(no_points_post_shift+this_ii+next_ii+ii_p_threshold>size(traces,2))
                    at_end=1;
                else
                    this_ii=this_ii+next_ii+delta_ii_sp-1;
                end
            end
        end

        if isempty(next_ii_sp)
            %This is S-


            next_ii=next_ii_sm;
            if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii>0)...
                    &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)&...
                    (this_ii+next_ii+ii_span<length(time))&(this_ii+next_ii-ii_span>=1)&(no_points_post_shift+this_ii+next_ii+ii_p_threshold<=size(traces,2))
                if epochs(this_ii+1+next_ii_sm)==9
                    cr_per_trial=[cr_per_trial 1];
                    fa_per_trial=[fa_per_trial 0];
                else
                    cr_per_trial=[cr_per_trial 0];
                    fa_per_trial=[fa_per_trial 1];
                end
                hit_per_trial=[hit_per_trial 0];
                miss_per_trial=[miss_per_trial 0];
                measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                training_decisions_post(ii_post+1:ii_post+no_points_post)=trimmed_licks(no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1);

                ii_sm_post=ii_sm_post+1;
                dFFs_sm_per_trial_per_ROI_post_shifted(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                dFFs_this_trial=traces(:,this_ii+next_ii-ii_span:this_ii+next_ii+ii_span);
                dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;
                trial_no=trial_no+1;
                decisions_per_trial(trial_no)=0;
                which_model_for_traces_loo(1,ii_which_model+1:this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                ii_which_model=this_ii+next_ii+no_points_post-1+dt_post_which_model;
                ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=this_ii+next_ii:this_ii+next_ii+no_points_post-1;
                epochs_sm_post(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
                %                 measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)';
                %                 epochs_sm_pre(1,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)=1;
                this_ii=this_ii+next_ii+delta_ii_sm-1;
                ii_post=ii_post+no_points_post;
                %                 ii_pre=ii_pre+no_points_pre;
            else
                if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||((no_points_post_shift+this_ii+next_ii-ii_span>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs)))||(no_points_post_shift+this_ii+next_ii+ii_p_threshold>size(traces,2))
                    at_end=1;
                else
                    this_ii=this_ii+next_ii+delta_ii_sm-1;
                end
            end
        end

        if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
            if next_ii_sm<next_ii_sp
                %This is S-
                next_ii=next_ii_sm;



                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)&...
                        (this_ii+next_ii-ii_span>=1)&(this_ii+next_ii+ii_span<length(time))&(no_points_post_shift+this_ii+next_ii+ii_p_threshold<=size(traces,2))
                    if epochs(this_ii+1+next_ii_sm)==9
                        cr_per_trial=[cr_per_trial 1];
                        fa_per_trial=[fa_per_trial 0];
                    else
                        cr_per_trial=[cr_per_trial 0];
                        fa_per_trial=[fa_per_trial 1];
                    end
                    hit_per_trial=[hit_per_trial 0];
                    miss_per_trial=[miss_per_trial 0];
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post)=trimmed_licks(no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1);
                    ii_sm_post=ii_sm_post+1;
                    dFFs_sm_per_trial_per_ROI_post_shifted(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,this_ii+next_ii-ii_span:this_ii+next_ii+ii_span);
                    dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;
                    trial_no=trial_no+1;
                    decisions_per_trial(trial_no)=0;
                    which_model_for_traces_loo(1,ii_which_model+1:this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                    ii_which_model=this_ii+next_ii+no_points_post-1+dt_post_which_model;
                    ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=this_ii+next_ii:this_ii+next_ii+no_points_post-1;
                    epochs_sm_post(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
                    %                     measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)';
                    %                     epochs_sm_pre(1,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)=1;
                    this_ii=this_ii+next_ii+delta_ii_sm-1;
                    ii_post=ii_post+no_points_post;
                    %                     ii_pre=ii_pre+no_points_pre;
                else
                    if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||((no_points_post_shift+this_ii+next_ii-ii_span>0)...
                            &(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs)))||(no_points_post_shift+this_ii+next_ii+ii_p_threshold>size(traces,2))
                        at_end=1;
                    else
                        this_ii=this_ii+next_ii+delta_ii_sm-1;
                    end
                end
            else
                %This is S+
                next_ii=next_ii_sp;



                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii>0)...
                        &(no_points_post_shift+this_ii+next_ii+ii_span<=length(epochs))&(no_points_post_shift+this_ii+next_ii-ii_span>0)&...
                        (this_ii+next_ii-ii_span>=1)&(this_ii+next_ii+ii_span<length(time))&(no_points_post_shift+this_ii+next_ii+ii_p_threshold<=size(traces,2))
                    if epochs(this_ii+1+next_ii_sp)==6
                        hit_per_trial=[hit_per_trial 1];
                        miss_per_trial=[miss_per_trial 0];
                    else
                        hit_per_trial=[hit_per_trial 0];
                        miss_per_trial=[miss_per_trial 1];
                    end
                    cr_per_trial=[cr_per_trial 0];
                    fa_per_trial=[fa_per_trial 0];
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post)=trimmed_licks(no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1);
                    ii_sp_post=ii_sp_post+1;
                    dFFs_sp_per_trial_per_ROI_post_shifted(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,this_ii+next_ii-ii_span:this_ii+next_ii+ii_span);
                    dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;
                    trial_no=trial_no+1;
                    decisions_per_trial(trial_no)=1;
                    which_model_for_traces_loo(1,ii_which_model+1:this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                    ii_which_model=this_ii+next_ii+no_points_post-1+dt_post_which_model;
                    ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=this_ii+next_ii:this_ii+next_ii+no_points_post-1;
                    epochs_sp_post(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
                    %                     measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)';
                    %                     epochs_sp_pre(1,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)=1;
                    this_ii=this_ii+next_ii+delta_ii_sp-1;
                    ii_post=ii_post+no_points_post;
                    %                     ii_pre=ii_pre+no_points_pre;
                else
                    if (no_points_post_shift+this_ii+next_ii+no_points_post>length(epochs))||((no_points_post_shift+this_ii+next_ii-ii_span>0)...
                            &(no_points_post_shift+this_ii+next_ii+ii_span>length(epochs)))||(no_points_post_shift+this_ii+next_ii+ii_p_threshold>size(traces,2))
                        at_end=1;
                    else
                        this_ii=this_ii+next_ii+delta_ii_sp-1;
                    end
                end
            end
        end


    end

end



pffft=1;
