function [percent_correct] = drgCaImAnFindPercentCorrect(pre_per_PathName,pre_per_FileName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load([pre_per_PathName pre_per_FileName])

hit_per_trial=[];
miss_per_trial=[];
cr_per_trial=[];
fa_per_trial=[];
at_end=0;
this_ii=0;


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
    next_ii_sm=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    
    if (isempty(next_ii_sp))&(isempty(next_ii_sm))
        at_end=1;
    else
        
        if isempty(next_ii_sm)
            %This is S+
            if epochs(this_ii+1+next_ii_sp)==6
                hit_per_trial=[hit_per_trial 1];
                miss_per_trial=[miss_per_trial 0];
            else
                hit_per_trial=[hit_per_trial 0];
                miss_per_trial=[miss_per_trial 1];
            end
            cr_per_trial=[cr_per_trial 0];
            fa_per_trial=[fa_per_trial 0];
            
            delta_next_ii_sp=find(epochs(this_ii+next_ii_sp:end)~=epochs(this_ii+1+next_ii_sp),1,'first');
            this_ii=this_ii+next_ii_sp+delta_next_ii_sp;
            
        end
        
        if isempty(next_ii_sp)
            %This is S-
            if epochs(this_ii+1+next_ii_sm)==9
                cr_per_trial=[cr_per_trial 1];
                fa_per_trial=[fa_per_trial 0];
            else
                cr_per_trial=[cr_per_trial 0];
                fa_per_trial=[fa_per_trial 1];
            end
            hit_per_trial=[hit_per_trial 0];
            miss_per_trial=[miss_per_trial 0];
            
            delta_next_ii_sm=find(epochs(this_ii+next_ii_sm:end)~=epochs(this_ii+1+next_ii_sm),1,'first');
            this_ii=this_ii+next_ii_sm+delta_next_ii_sm;
            
            
        end
        
        if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
            if next_ii_sm<next_ii_sp
                %This is S-
                next_ii=next_ii_sm;
                
                if epochs(this_ii+1+next_ii_sm)==9
                    cr_per_trial=[cr_per_trial 1];
                    fa_per_trial=[fa_per_trial 0];
                else
                    cr_per_trial=[cr_per_trial 0];
                    fa_per_trial=[fa_per_trial 1];
                end
                hit_per_trial=[hit_per_trial 0];
                miss_per_trial=[miss_per_trial 0];
                
                delta_next_ii_sm=find(epochs(this_ii+next_ii_sm:end)~=epochs(this_ii+1+next_ii_sm),1,'first');
                this_ii=this_ii+next_ii_sm+delta_next_ii_sm;
                
                
            else
                %This is S+
                next_ii=next_ii_sp;
                
                if epochs(this_ii+1+next_ii_sp)==6
                    hit_per_trial=[hit_per_trial 1];
                    miss_per_trial=[miss_per_trial 0];
                else
                    hit_per_trial=[hit_per_trial 0];
                    miss_per_trial=[miss_per_trial 1];
                end
                cr_per_trial=[cr_per_trial 0];
                fa_per_trial=[fa_per_trial 0];
                
                delta_next_ii_sp=find(epochs(this_ii+next_ii_sp:end)~=epochs(this_ii+1+next_ii_sp),1,'first');
                this_ii=this_ii+next_ii_sp+delta_next_ii_sp;
                
            end
        end
        
        
    end
    
end

%Calculate percent correct
percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);

end

