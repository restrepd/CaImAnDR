function epochs=drgCaImAn_find_trials_in_eopchs(epochs, ii_shuffled)
%This finds the Sp and Sm and changes the epochs to a different point in the fellows series 
ii=1;
not_at_end=1;
spsm_trials=[];
trial_from=[];
trial_to=[];
ii_trials=0;
while not_at_end==1
    this_ii_sp=find((epochs(ii:end)==6)|(epochs(ii:end)==7),1,'first');
    this_ii_sm=find((epochs(ii:end)==8)|(epochs(ii:end)==9),1,'first');


    %Neither is empty
    if ((~isempty(this_ii_sp)&(~isempty(this_ii_sm))))

        %Next is an Sp
        if (this_ii_sp<this_ii_sm)|((isempty(this_ii_sm)&(~isempty(this_ii_sp))))
            ii_trials=ii_trials+1;
            spsm_trials(ii_trials)=1;
            trial_from(ii_trials)=ii+this_ii_sp-1;
            this_epoch=epochs(ii+this_ii_sp-1);
            this_ii_sp_off=find((epochs(ii+this_ii_sp-1:end)~=this_epoch),1,'first');
            trial_to(ii_trials)=ii+this_ii_sp-1+this_ii_sp_off-2;
            ii=ii+this_ii_sp-1+this_ii_sp_off-1;
        end

        %Next is an Sm
        if (this_ii_sp>this_ii_sm)
            ii_trials=ii_trials+1;
            spsm_trials(ii_trials)=0;
            trial_from(ii_trials)=ii+this_ii_sm-1;
            this_epoch=epochs(ii+this_ii_sm-1);
            this_ii_off=find((epochs(ii+this_ii_sm-1:end)~=this_epoch),1,'first');
            trial_to(ii_trials)=ii+this_ii_sm-1+this_ii_off-2;
            ii=ii+this_ii_sm-1+this_ii_off-1;
        end
    end

    %Sp is empty
    if ((isempty(this_ii_sp)&(~isempty(this_ii_sm))))
        ii_trials=ii_trials+1;
        spsm_trials(ii_trials)=0;
        trial_from(ii_trials)=ii+this_ii_sm-1;
        this_epoch=epochs(ii+this_ii_sm-1);
        this_ii_off=find((epochs(ii+this_ii_sm-1:end)~=this_epoch),1,'first');
        trial_to(ii_trials)=ii+this_ii_sm-1+this_ii_off-2;
        ii=ii+this_ii_sm-1+this_ii_off-1;
    end

    %Sm is empty
    if ((isempty(this_ii_sm)&(~isempty(this_ii_sp))))
        ii_trials=ii_trials+1;
        spsm_trials(ii_trials)=1;
        trial_from(ii_trials)=ii+this_ii_sp-1;
        this_epoch=epochs(ii+this_ii_sp-1);
        this_ii_sp_off=find((epochs(ii+this_ii_sp-1:end)~=this_epoch),1,'first');
        trial_to(ii_trials)=ii+this_ii_sp-1+this_ii_sp_off-2;
        ii=ii+this_ii_sp-1+this_ii_sp_off-1;
    end


    if isempty(this_ii_sp)&isempty(this_ii_sm)
        not_at_end=0;
    end
end

switch ii_shuffled
    case 1
        [randomFellows randomOpto]=dropcGetSlotnickOdorList();
        ii_start=11;
    case 2
        rng;
        these_rand=rand(1,100);
        randomFellows=these_rand>0.5;
        ii_start=1;
end


for ii_trial=1:ii_trials
    for ii=trial_from(ii_trial):trial_to(ii_trial)
        if randomFellows(ii_start)==0
            epochs(ii)=7;
        else
            epochs(ii)=9;
        end
    end
    ii_start=ii_start+1;
end