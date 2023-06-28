function handles_out2=drgCaImAn_SVZ_entire_session_selectROI(handles_choices)
%This program trains several decoding algorithms with the post odorant and then determines what happens throughout the entire timecouse
%The user enters the choices entered under exist('handles_choices')==0
%
% processing_algorithm= 1 and 2 were used for troublehsooting and do not
% produce reliable results because of overtraining, use
% processing_algoritm=3, that was vetted for our manuscript for regular
% decoding
%
%
% the input is a pre_per file version 2

if exist('handles_choices')==0
    clear all
    close all

    %Load file
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the *pre_per.mat');

    [pre_per_rdecFileName,pre_per_rdecPathName] = uigetfile({'*pre_per_rdec.mat'},'Select the *pre_per_rdec.mat');

    process_low=0; %1= process with ROIs with accuracy<=0.35 0= process with accuracy >=0.65
    processing_algorithm=3; %Use 3 for manuscript (trained with all points in the training window)
    k_fold=5; %Only used for processing_algorithm=2,
    post_time=5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on, 5
    post_shift=0.5; %Set to 0 if you want to train with odor on points, 0
    pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
    MLalgo_to_use=6; %Vector with the decoding algorithms you want to use
    ii_cost=3;
    p_threshold=1.1; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
    dt_p_threshold=5; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures

    no_ROI_draws=1000; %Number of times that decoding is calculated for each set of no_ROIs
    no_ROIs=1; %Number of ROIs used in the decoding (sampled randomly from the total number of ROIs)
    fileNo=1;

    handles_choices.pre_perFileName=pre_perFileName;
    handles_choices.pre_perPathName=pre_perPathName;
    handles_choices.processing_algorithm=processing_algorithm;
    handles_choices.post_time=post_time;
    handles_choices.k_fold=k_fold;
    handles_choices.post_shift=post_shift;
    handles_choices.MLalgo_to_use=MLalgo_to_use;
    handles_choices.pre_time=pre_time;
    handles_choices.p_threshold=p_threshold;
    handles_choices.dt_p_threshold=dt_p_threshold;
    handles_choices.show_figures=show_figures;
    handles_choices.ii_cost=ii_cost;
    handles_choices.no_ROI_draws=no_ROI_draws;
    handles_choices.process_low=process_low;

else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
    pre_per_rdecPathName=pre_perPathName;
    pre_per_rdecFileName=[pre_perFileName(1:end-4) '_rdec.mat'];
    processing_algorithm=handles_choices.processing_algorithm;
    post_time=handles_choices.post_time;
    k_fold=handles_choices.k_fold;
    post_shift=handles_choices.post_shift;
    MLalgo_to_use=handles_choices.MLalgo_to_use;
    pre_time=handles_choices.pre_time;
    p_threshold=handles_choices.p_threshold;
    dt_p_threshold=handles_choices.dt_p_threshold;
    show_figures=handles_choices.show_figures;
    ii_cost=handles_choices.ii_cost;
    no_ROI_draws=handles_choices.no_ROI_draws;
    no_ROIs=handles_choices.no_ROIs;
    process_low=handles_choices.process_low;
    fileNo=handles_choices.fileNo;
end

w = warning ('off','all');
handles_out2.handles_choices=handles_choices;

% warning('off')

tic

%Restart random seeds
rng('shuffle');


convert_z=1; %Convert dFF traces to z
dt_span=15; %Seconds shown before and after odor on in the figures, 40 was used before
moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
no_shuffles=3; %Number of shuffles for per trial shuffling
window_no=2;
time_windows=[-1 0;
    3.1 4.1];

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;
  
handles_out2.failed=0;

%Load the rdec file and choose the ROIs for decoding
if exist([pre_per_rdecPathName pre_per_rdecFileName])==0
    handles_out2.failed=1;
    fprintf(1, ['File no %d does not exist\n'],fileNo);
else
    load([pre_per_rdecPathName pre_per_rdecFileName])
    ii_out=1;
    handles_out2rdec=handles_out.ii_out(ii_out).handles_out;
    if ~isfield(handles_out2rdec,'no_ROI_draws')
        handles_out2.failed=1;
        fprintf(1, ['File no %d must be re-run\n'],fileNo);
    end
    if show_figures==1
        fprintf(1, ['\ndrgCaImAnInspectMultiROI run for ' pre_perFileName '\n\n']);
    end  
end


if handles_out2.failed==0
    handles_out2.failed=0;
    no_ROI_draws=handles_out2rdec.no_ROI_draws;
    rdecMLalgo=6;
    rdectime_windows=[3.1 4.1];
    rdectime_span=handles_out2rdec.time_span;

    rdec_accuracy_per_ROI=[];

    for iiROI=1:no_ROI_draws
        rdec_accuracy_per_ROI=[rdec_accuracy_per_ROI mean(mean(handles_out2rdec.ROI(iiROI).MLalgo(rdecMLalgo).this_correct_predict(:,(rdectime_span>=rdectime_windows(1))&(rdectime_span<=rdectime_windows(2))),2))];
    end

    no_ROI_draws=1; %Number of times that decoding is calculated for each set of no_ROIs
    no_ROIs=10000; %Number of ROIs used in the decoding (sampled randomly from the total number of ROIs)

    load([pre_perPathName pre_perFileName])
    if show_figures==1
        fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
        fprintf(1, 'post_time = %d, p_threshold= %d, post_shift= %d, cost %d\n',post_time,p_threshold,post_shift, ii_cost);
    end

    if convert_z==1
        for trace_no=1:size(traces,1)
            traces(trace_no,:)=traces(trace_no,:)/std(traces(trace_no,:));
        end
    end

    classifier_names{1}='Linear Discriminant';
    classifier_names{2}='Support Vector Machine';
    classifier_names{3}='Naive Bayes Classifier';
    classifier_names{4}='Neural Network';
    classifier_names{5}='Decision tree';
    classifier_names{6}='Binomial glm';

    handles_out2.handles_choices=handles_choices;
    handles_out2.pre_perFileName=pre_perFileName;
    handles_out2.pre_perPathName=pre_perPathName;
    handles_out2.post_time=post_time;
    handles_out2.k_fold=k_fold;
    handles_out2.post_shift=post_shift;
    handles_out2.pre_time=pre_time;
    handles_not_out.MLalgo_to_use=MLalgo_to_use;

    figNo=0;

    %time has the time for the dF/F traces(ROI,time)
    if show_figures==1
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        hold on

        % Determine the y spacing of the traces
        y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));

        %Plot the traces and do z normalization
        %For S+ and S- plot odor on and reinforcement
        for epoch=1:handles.dropcData.epochIndex
            %Epoch 2 is odor on, 3 is odor off
            plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
            if plot_epoch
                if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                    plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                        '-r','LineWidth',1)
                else
                    plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                        '-b','LineWidth',1)
                end
            end
        end

        for trNo=1:no_traces
            % for trNo=1:20
            plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        end

        ylim([-y_shift*0.2 (no_traces+2)*y_shift])
        xlabel('time(sec)')
        title(['All dFF timecourses ' num2str(size(traces,1)) ' ROIs'])

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        hold on

        % Determine the y spacing of the traces
        if process_low==1
            these_traces=traces(rdec_accuracy_per_ROI<=0.35,:);
        else
            these_traces=traces(rdec_accuracy_per_ROI>=0.65,:);
        end
        y_shift=1.2*(prctile(traces(:),95)-prctile(these_traces(:),5));

        %Plot the traces and do z normalization
        %For S+ and S- plot odor on and reinforcement
        for epoch=1:handles.dropcData.epochIndex
            %Epoch 2 is odor on, 3 is odor off
            plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
            if plot_epoch
                if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                    plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                        '-r','LineWidth',1)
                else
                    plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                        '-b','LineWidth',1)
                end
            end
        end

        trNo_included=0;
        for trNo=1:no_traces
            if process_low==1
                if rdec_accuracy_per_ROI(trNo)<=0.35
                    trNo_included=trNo_included+1;
                    plot(time,traces(trNo,:)+y_shift*trNo_included,'-k','LineWidth',1)
                end
            else
                if rdec_accuracy_per_ROI(trNo)>=0.65
                    trNo_included=trNo_included+1;
                    plot(time,traces(trNo,:)+y_shift*trNo_included,'-k','LineWidth',1)
                end
            end
        end
        no_these_traces=size(these_traces,1);
        ylim([-y_shift*0.2 (no_these_traces+2)*y_shift])
        xlabel('time(sec)')
        if process_low==1
            title(['dFF timecourses for ROIs with accuracy<=0.35'])
        else
            title(['dFF timecourses for ROIs with accuracy>=0.65'])
        end
    end

    %if the number of ROIs is larger than the number available decrease it to
    %the number of ROIs
    if no_ROIs>no_traces
        no_ROIs=no_traces;
        no_ROI_draws=1;
    end

    if no_ROIs==1
        no_ROI_draws=no_traces;
    end
    %epochs is a vector of the length of time that gives information on
    %behavior
    % 1=Final Valve
    % 6=Hit (on for the duration of odor on)
    % 7=Miss
    % 8=FA
    % 9=CR

    %For example Hit||Miss shows S+ odor application times (red)
    %and FA||CR gives S- (blue)

    %Post points
    Nall=size(traces,1);
    dt=median(time(2:end)-time(1:end-1));
    ii_p_threshold=ceil(dt_p_threshold/dt);
    no_points_post=floor(post_time/dt);
    no_points_post_shift=floor(post_shift/dt);
    no_points_pre=floor(pre_time/dt);


    measurements_post=[];
    measurements_pre=[];
    epochs_sp_post=zeros(1,length(time));
    epochs_sp_pre=zeros(1,length(time));
    epochs_sm_post=zeros(1,length(time));
    epochs_sm_pre=zeros(1,length(time));
    which_model_for_traces_loo=no_odor_trials*ones(1,size(traces,2));



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
    dFFs_sp_per_trial_per_ROI=[];
    dFFs_sm_per_trial_per_ROI=[];
    hit_per_trial=[];
    miss_per_trial=[];
    cr_per_trial=[];
    fa_per_trial=[];

    [hit_per_trial,cr_per_trial,dFFs_sp_per_trial_per_ROI,...
        dFFs_sm_per_trial_per_ROI,dFF_per_trial_sm,dFF_per_trial_sp,training_decisions_post,...
        which_model_for_traces_loo,decisions_per_trial,...
        ii_pointer_to_td,epochs_sp_post,measurements_post,...
        measurements_pre,epochs_sp_pre,ii_post,trial_no...
        ,epochs_sm_post,epochs_sm_pre] = ...
        drgCaImAn_parse_out_trials(dt, dt_span,epochs,no_points_post_shift,no_points_post,no_points_pre,traces,ii_p_threshold,no_odor_trials);


    %Calculate percent correct
    handles_out.percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);
    if show_figures==1
        fprintf(1,'Percent corect %d\n',handles_out.percent_correct)
    end

    which_model_for_traces_loo(which_model_for_traces_loo>trial_no)=trial_no;
    %
    % %Do S+
    % at_end=0;
    % this_ii=0;
    % ii_post=0;
    % ii_pre=0;
    % ii=0;
    % trial_no=0;
    %
    % while (at_end==0)
    %     next_ii=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    %     if ~isempty(next_ii)
    %         if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
    %             measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
    %             trial_no=trial_no+1;
    %             which_model_for_traces_loo(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=trial_no;
    %             ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
    %             epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
    %             measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
    %             epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
    %             this_ii=this_ii+next_ii+no_points_post;
    %             ii_post=ii_post+no_points_post;
    %             ii_pre=ii_pre+no_points_pre;
    %         else
    %             at_end=1;
    %         end
    %     else
    %         at_end=1;
    %     end
    % end
    %
    % training_decisions_post=ones(1,size(measurements_post,1));
    % ii_sp_post=size(measurements_post,1);
    %
    %
    % %Do S-
    % at_end=0;
    % this_ii=0;
    % ii=0;
    %
    %
    % while (at_end==0)
    %     next_ii=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    %     if ~isempty(next_ii)
    %         if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
    %             measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
    %             trial_no=trial_no+1;
    %             which_model_for_traces_loo(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=trial_no;
    %             ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
    %             epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
    %             measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
    %             epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
    %             this_ii=this_ii+next_ii+no_points_post;
    %             ii_post=ii_post+no_points_post;
    %             ii_pre=ii_pre+no_points_pre;
    %         else
    %             at_end=1;
    %         end
    %     else
    %         at_end=1;
    %     end
    % end
    %
    % ii_sm_post=ii_post-ii_sp_post;
    %
    % training_decisions_post=[training_decisions_post zeros(1,size(measurements_post,1)-ii_sp_post)];

    %Now let's limit the ROIs to those below p_threshold

    p_values=ones(1,size(dFFs_sm_per_trial_per_ROI,2));
    for iiROI=1:size(dFFs_sm_per_trial_per_ROI,2)
        dFF_sm=zeros(size(dFFs_sm_per_trial_per_ROI,1),size(dFFs_sm_per_trial_per_ROI,3));
        dFF_sm(:,:)=dFFs_sm_per_trial_per_ROI(:,iiROI,:);
        dFF_sp=zeros(size(dFFs_sp_per_trial_per_ROI,1),size(dFFs_sp_per_trial_per_ROI,3));
        dFF_sp(:,:)=dFFs_sp_per_trial_per_ROI(:,iiROI,:);

        [h,p_values(iiROI)]=ttest2(mean(dFF_sp,2),mean(dFF_sm,2));
    end

    p_value_masks=[];
    for iiROI=1:no_ROI_draws

        %     if no_ROIs~=1
        %         found_it=0;
        %         while found_it==0
        %             these_mask_ROIs=randperm(length(p_values));
        %             these_mask_ROIs=these_mask_ROIs(1:no_ROIs);
        %             p_value_mask=logical(zeros(1,length(p_values)));
        %             p_value_mask(these_mask_ROIs)=1;
        %             found_it=1;
        %             for jjROI=1:size(p_value_masks,1)
        %                 if sum(p_value_masks(jjROI,:)==p_value_mask)==length(p_values)
        %                     found_it=0;
        %                 end
        %             end
        %         end
        %     else
        %         p_value_mask=logical(zeros(1,length(p_values)));
        %         p_value_mask(1,iiROI)=logical(1);
        %     end

        if process_low==1
            p_value_mask=rdec_accuracy_per_ROI<=0.35;
        else
            p_value_mask=rdec_accuracy_per_ROI>=0.65;
        end

        p_value_masks(iiROI,:)=p_value_mask;
        %     end

        handles_out2.p_value_masks=p_value_masks;

        %Trim the number of ROIs in all matrices
        noROIs_before_trimming=size(measurements_post,2);
        %     dFF_per_trial_sp=dFF_per_trial_sp(:,p_value_mask,:);
        %     dFF_per_trial_sm=dFF_per_trial_sm(:,p_value_mask,:);
        measurements_post_trimmed=measurements_post(:,p_value_mask);
        %     measurements_pre_trimmed=measurements_pre(:,p_value_mask);
        traces_trimmed=traces(p_value_mask,:);
        %     no_traces=size(traces,1);

        %Save odor times
        handles_out2.sp_times=[];
        handles_out2.sp_times_ii=0;
        handles_out2.sm_times=[];
        handles_out2.sm_times_ii=0;

        %For S+ and S- plot odor on and reinforcement
        for epoch=1:handles.dropcData.epochIndex
            %Epoch 2 is odor on, 3 is odor off
            plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
            if plot_epoch
                if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                    handles_out2.sp_times_ii=handles_out2.sp_times_ii+1;
                    handles_out2.sp_times(handles_out2.sp_times_ii)=handles.dropcData.epochTime(epoch);
                    %             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces_trimmed+2)*y_shift],...
                    %                 '-r','LineWidth',1)
                else
                    handles_out2.sm_times_ii=handles_out2.sm_times_ii+1;
                    handles_out2.sm_times(handles_out2.sm_times_ii)=handles.dropcData.epochTime(epoch);
                    %             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces_trimmed+2)*y_shift],...
                    %                 '-b','LineWidth',1)
                end
            end
        end


        training_decisions_post_sh=zeros(1,length(training_decisions_post));
        training_decisions_post_sh(1,:)=training_decisions_post(randperm(length(training_decisions_post)));

        for ii=1:no_shuffles
            these_shuffled_decisions_per_trial=decisions_per_trial(randperm(length(decisions_per_trial)));
            ww=0;
            for jj=1:length(decisions_per_trial)
                training_decisions_post_sh2(ii,ww+1:ww+no_points_post)=these_shuffled_decisions_per_trial(jj)*ones(1,no_points_post);
                ww=ww+no_points_post;
            end
        end


        handles_out2.dt=dt;

        fprintf(1, ['Training run number %d with %d ROIs (original no ROIs %d)...\n'],iiROI,size(measurements_post_trimmed,2),noROIs_before_trimming);


        for MLalgo=MLalgo_to_use

            this_cost=[0 ii_cost;ii_cost 0];
            labels=[];
            timepoint_processed=[];
            correct_predict=[];
            correct_predict_shuffled=[];



            Nall_post=size(measurements_post_trimmed,1);



            switch processing_algorithm
                case 1
                    %Train with all the post data
                    %Store the training data in a table.
                    tblTrn=[];
                    tblTrn = array2table(measurements_post_trimmed);

                    %Store the decisions in Y
                    Y=training_decisions_post;

                    switch MLalgo
                        case 1
                            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        case 2
                            Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                        case 4
                            Mdl = fitcnet(tblTrn,Y);
                        case 5
                            Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                    end

                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason tisdid not work for net when I did this:
                    % [label_traces_trimmed,score] = predict(Mdl,traces_trimmed');
                    % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                    %I had to resolrt to the for loop:
                    switch MLalgo
                        case 3,4
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label_traces(ii),score] = predict(Mdl,this_time_point);
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label_post(ii),score] = predict(Mdl,this_time_point);
                            end
                        otherwise
                            [label_traces,score] = predict(Mdl,traces_trimmed');
                            [label_post,score] = predict(Mdl,measurements_post_trimmed);
                    end

                case 2
                    %k-fold cross validation evaluating performance with left-out data
                    %Store the training data in a table.

                    handles_not_out.MLalgo(MLalgo).models=[];
                    handles_not_out.MLalgo(MLalgo).processed_succesfully=1;
                    points_masked=floor(no_points_post/k_fold);
                    which_model=ones(1,size(measurements_post_trimmed,1));
                    for kk=1:k_fold


                        training_mask=ones(size(measurements_post_trimmed,1),1);
                        at_end=0;
                        ii=0;

                        while at_end==0
                            if ii+no_points_post<=size(measurements_post_trimmed,1)
                                %This is what I originally used, these samples are adjacent
                                %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                                %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                                %
                                %I then changed it to this, spaced by delta_ii
                                delta_ii=floor(no_points_post/points_masked);
                                for ii_pm=1:points_masked
                                    training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                    which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                                end

                                ii=ii+no_points_post;
                            else
                                at_end=1;
                            end
                        end

                        these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                        these_training_decisions=zeros(1,sum(training_mask));


                        jj=0;
                        for ii=1:size(measurements_post_trimmed,1)
                            if training_mask(ii)==1
                                jj=jj+1;
                                these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                                these_training_decisions(jj)=training_decisions_post(ii);
                            end
                        end

                        tblTrn=[];
                        tblTrn = array2table(these_training_measurements);

                        %Store the decisions in Y
                        Y=these_training_decisions;

                        switch MLalgo
                            case 1
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                            case 2
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                            case 3
                                %The try catch was entered here because of this
                                %error
                                % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                % A normal distribution cannot be fit for the combination of class 1 and predictor
                                % these_training_measurements107. The data has zero variance.
                                try
                                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                catch
                                    handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                                end
                            case 4
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                            case 5
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                            case 6
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                        end
                    end

                    if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
                        %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                        %Note: for some reason tisdid not work for net when I did this:
                        % [label_traces,score] = predict(Mdl,traces');
                        % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                        %I had to resolrt to the for loop:

                        which_model_for_traces=randi(k_fold,1,size(traces_trimmed,2));

                        for ii=1:length(which_model)
                            which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                        end

                        if MLalgo==6
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_traces(ii)=1;
                                else
                                    label_traces(ii)=0;
                                end
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_post(ii)=1;
                                else
                                    label_post(ii)=0;
                                end
                            end
                        else
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label_traces(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces(ii)).Mdl,this_time_point);
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label_post(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                            end
                        end


                        handles_not_out.MLalgo(MLalgo).label_post=label_post;
                        handles_not_out.MLalgo(MLalgo).label_traces=label_traces;

                        %Now do prediction with shuffled training decisions
                        %k-fold cross validation evaluating performance with left-out data
                        %Store the training data in a table.

                        handles_not_out.MLalgo(MLalgo).sh_models=[];
                        points_masked=floor(no_points_post/k_fold);
                        which_model=ones(1,size(measurements_post_trimmed,1));
                        for kk=1:k_fold


                            training_mask=ones(size(measurements_post_trimmed,1),1);
                            at_end=0;
                            ii=0;

                            while at_end==0
                                if ii+no_points_post<=size(measurements_post_trimmed,1)
                                    %This is what I originally used, these samples are adjacent
                                    %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                                    %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                                    %
                                    %I then changed it to this, spaced by delta_ii
                                    delta_ii=floor(no_points_post/points_masked);
                                    for ii_pm=1:points_masked
                                        training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                        which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                                    end

                                    ii=ii+no_points_post;
                                else
                                    at_end=1;
                                end
                            end

                            these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                            these_training_decisions=zeros(1,sum(training_mask));


                            jj=0;
                            for ii=1:size(measurements_post_trimmed,1)
                                if training_mask(ii)==1
                                    jj=jj+1;
                                    these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                                    these_training_decisions(jj)=training_decisions_post_sh(ii);
                                end
                            end

                            tblTrn=[];
                            tblTrn = array2table(these_training_measurements);

                            %Store the decisions in Y
                            Y=these_training_decisions;

                            handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;

                            switch MLalgo
                                case 1
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                case 2
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                case 3
                                    %The try catch was entered here because of this
                                    %error
                                    % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                    % A normal distribution cannot be fit for the combination of class 1 and predictor
                                    % these_training_measurements107. The data has zero variance.
                                    try
                                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                                    end
                                case 4
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                                case 5
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                case 6
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                            end
                        end

                        %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                        %Note: for some reason this did not work for net when I did this:
                        % [label_traces,score] = predict(Mdl,traces');
                        % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                        %I had to resort to the for loop:

                        which_model_for_traces=randi(k_fold,1,size(traces_trimmed,2));

                        for ii=1:length(which_model)
                            which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                        end

                        if MLalgo==6

                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_traces_sh(ii)=1;
                                else
                                    label_traces_sh(ii)=0;
                                end
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_post_sh(ii)=1;
                                else
                                    label_post_sh(ii)=0;
                                end
                            end
                        else
                            label_traces_sh=zeros(1,size(traces_trimmed,2));
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label_traces_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces(ii)).Mdl,this_time_point);
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label_post_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                            end
                        end

                        handles_not_out.MLalgo(MLalgo).label_post=label_post;
                        handles_not_out.MLalgo(MLalgo).label_traces=label_traces;

                        handles_not_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
                        handles_not_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;

                        %Now do shuffling on a per trial basis sh2
                        for ii_sh=1:no_shuffles


                            %Now do prediction with shuffled training decisions
                            %k-fold cross validation evaluating performance with left-out data
                            %Store the training data in a table.
                            %                         which_model
                            handles_not_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                            points_masked=floor(no_points_post/k_fold);
                            which_model=ones(1,size(measurements_post_trimmed,1));
                            for kk=1:k_fold


                                training_mask=ones(size(measurements_post_trimmed,1),1);
                                at_end=0;
                                ii=0;

                                while at_end==0
                                    if ii+no_points_post<=size(measurements_post_trimmed,1)
                                        %This is what I originally used, these samples are adjacent
                                        %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                                        %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                                        %
                                        %I then changed it to this, spaced by delta_ii
                                        delta_ii=floor(no_points_post/points_masked);
                                        for ii_pm=1:points_masked
                                            training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                            which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                                        end

                                        ii=ii+no_points_post;
                                    else
                                        at_end=1;
                                    end
                                end

                                these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                                these_training_decisions=zeros(1,sum(training_mask));


                                jj=0;
                                for ii=1:size(measurements_post_trimmed,1)
                                    if training_mask(ii)==1
                                        jj=jj+1;
                                        these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                                        these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                                    end
                                end

                                tblTrn=[];
                                tblTrn = array2table(these_training_measurements);

                                %Store the decisions in Y
                                Y=these_training_decisions;

                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;

                                switch MLalgo
                                    case 1
                                        try
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                        catch
                                            % Error using ClassificationDiscriminant (line 380)
                                            % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                            % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                        end
                                    case 2
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                    case 3
                                        %The try catch was entered here because of this
                                        %error
                                        % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                        % A normal distribution cannot be fit for the combination of class 1 and predictor
                                        % these_training_measurements107. The data has zero variance.
                                        try
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                        catch
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                                        end
                                    case 4
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                                    case 5
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                    case 6
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                                end
                            end

                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                            %Note: for some reason this did not work for net when I did this:
                            % [label_traces,score] = predict(Mdl,traces');
                            % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                            %I had to resort to the for loop:

                            which_model_for_traces=randi(k_fold,1,size(traces_trimmed,2));

                            for ii=1:length(which_model)
                                which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                            end

                            if MLalgo==6



                                for ii=1:size(measurements_post_trimmed,1)
                                    this_time_point=zeros(1,size(traces_trimmed,1));
                                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                                    if label>0.5
                                        label_post_sh2(ii_sh,ii)=1;
                                    else
                                        label_post_sh2(ii_sh,ii)=0;
                                    end
                                end
                            else


                                for ii=1:size(measurements_post_trimmed,1)
                                    this_time_point=zeros(1,size(traces_trimmed,1));
                                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                    [label_post_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                                end
                            end

                            handles_not_out.MLalgo(MLalgo).label_post=label_post_sh2;

                        end


                    end
                case 3
                    %leave one trial out
                    %Store the training data in a table.


                    handles_not_out.MLalgo(MLalgo).models=[];
                    handles_not_out.MLalgo(MLalgo).processed_succesfully=1;
                    points_masked=floor(no_points_post/k_fold);
                    which_model=ones(1,size(measurements_post_trimmed,1));
                    no_trials=ii_post/no_points_post;
                    for kk=1:no_trials


                        training_mask=ones(size(measurements_post_trimmed,1),1);
                        training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                        which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;


                        these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                        these_training_decisions=zeros(1,sum(training_mask));

                        jj=0;
                        for ii=1:size(measurements_post_trimmed,1)
                            if training_mask(ii)==1
                                jj=jj+1;
                                these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                                these_training_decisions(jj)=training_decisions_post(ii);
                            end
                        end

                        tblTrn=[];
                        tblTrn = array2table(these_training_measurements);

                        %Store the decisions in Y
                        Y=these_training_decisions;

                        switch MLalgo
                            case 1
                                try
                                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                catch
                                    % Error using ClassificationDiscriminant (line 380)
                                    % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                    % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                end
                            case 2
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                            case 3
                                %The try catch was entered here because of this
                                %error
                                % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                % A normal distribution cannot be fit for the combination of class 1 and predictor
                                % these_training_measurements107. The data has zero variance.
                                try
                                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                catch
                                    handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                                end
                            case 4
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                            case 5
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                            case 6
                                handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                        end
                    end

                    if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
                        %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                        %Note: for some reason this did not work for net when I did this:
                        % [label_traces,score] = predict(Mdl,traces');
                        % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                        %I had to resort to the for loop:


                        if MLalgo==6
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                scores(ii,:)=score;
                                if label>0.5
                                    label_traces(ii)=1;
                                else
                                    label_traces(ii)=0;
                                end
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                                scores_post(ii,:)=score;
                                if label>0.5
                                    label_post(ii)=1;
                                else
                                    label_post(ii)=0;
                                end
                            end
                        else
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                try
                                    [label_traces(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                    scores(ii,:)=score;
                                catch
                                    scores(ii,:)=rand(1,2);
                                    if rand(1)>0.5
                                        label_traces(ii)=1;
                                    else
                                        label_traces(ii)=0;
                                    end
                                end

                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                try
                                    [label_post(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                                    scores_post(ii,:)=score;
                                catch
                                    scores_post(ii,:)=rand(1,2);
                                    if rand(1)>0.5
                                        label_post(ii)=1;
                                    else
                                        label_post(ii)=0;
                                    end
                                end

                            end
                        end


                        handles_not_out.MLalgo(MLalgo).label_post=label_post;
                        handles_out.MLalgo(MLalgo).label_traces=label_traces;
                        handles_out.MLalgo(MLalgo).scores=scores;
                        handles_not_out.MLalgo(MLalgo).scores_post=scores_post;

                        %Now do prediction with shuffled training decisions
                        %k-fold cross validation evaluating performance with left-out data
                        %Store the training data in a table.

                        handles_not_out.MLalgo(MLalgo).sh_models=[];
                        points_masked=floor(no_points_post/k_fold);
                        which_model=ones(1,size(measurements_post_trimmed,1));
                        for kk=1:no_trials


                            training_mask=ones(size(measurements_post_trimmed,1),1);
                            training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                            which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;


                            these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                            these_training_decisions=zeros(1,sum(training_mask));


                            jj=0;
                            for ii=1:size(measurements_post_trimmed,1)
                                if training_mask(ii)==1
                                    jj=jj+1;
                                    these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                                    these_training_decisions(jj)=training_decisions_post_sh(ii);
                                end
                            end

                            tblTrn=[];
                            tblTrn = array2table(these_training_measurements);

                            %Store the decisions in Y
                            Y=these_training_decisions;

                            handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;

                            switch MLalgo
                                case 1
                                    try
                                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        % Error using ClassificationDiscriminant (line 380)
                                        % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                        % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                    end
                                case 2
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                case 3

                                    %The try catch was entered here because of this
                                    %error
                                    % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                    % A normal distribution cannot be fit for the combination of class 1 and predictor
                                    % these_training_measurements107. The data has zero variance.
                                    try
                                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                                    end
                                case 4
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                                case 5
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                case 6
                                    handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                            end
                        end

                        %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                        %Note: for some reason this did not work for net when I did this:
                        % [label_traces,score] = predict(Mdl,traces');
                        % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                        %I had to resort to the for loop:



                        if MLalgo==6

                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                scores_sh(ii,:)=score;
                                if label>0.5
                                    label_traces_sh(ii)=1;
                                else
                                    label_traces_sh(ii)=0;
                                end
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                                scores_post_sh(ii,:)=score;
                                if label>0.5
                                    label_post_sh(ii)=1;
                                else
                                    label_post_sh(ii)=0;
                                end
                            end
                        else
                            label_traces_sh=zeros(1,size(traces_trimmed,2));
                            for ii=1:size(traces_trimmed,2)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=traces_trimmed(:,ii);
                                [label_traces_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                scores_sh(ii,:)=score;
                            end

                            for ii=1:size(measurements_post_trimmed,1)
                                this_time_point=zeros(1,size(traces_trimmed,1));
                                this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                [label_post_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                                scores_post_sh(ii,:)=score;
                            end
                        end



                        handles_not_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
                        handles_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
                        handles_out.MLalgo(MLalgo).scores_sh=scores_sh;
                        handles_not_out.MLalgo(MLalgo).scores_post_sh=scores_post_sh;

                        %Now do shuffling on a per trial basis
                        for ii_sh=1:no_shuffles


                            %Now do prediction with shuffled training decisions
                            %k-fold cross validation evaluating performance with left-out data
                            %Store the training data in a table.

                            handles_not_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                            points_masked=floor(no_points_post/k_fold);
                            which_model=ones(1,size(measurements_post_trimmed,1));
                            for kk=1:no_trials


                                training_mask=ones(size(measurements_post_trimmed,1),1);
                                training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                                which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;

                                these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                                these_training_decisions=zeros(1,sum(training_mask));


                                jj=0;
                                for ii=1:size(measurements_post_trimmed,1)
                                    if training_mask(ii)==1
                                        jj=jj+1;
                                        these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                                        these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                                    end
                                end

                                tblTrn=[];
                                tblTrn = array2table(these_training_measurements);

                                %Store the decisions in Y
                                Y=these_training_decisions;

                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;

                                switch MLalgo
                                    case 1
                                        try
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                        catch
                                            % Error using ClassificationDiscriminant (line 380)
                                            % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                            % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl  = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                        end
                                    case 2
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                    case 3
                                        %The try catch was entered here because of this
                                        %error
                                        % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                        % A normal distribution cannot be fit for the combination of class 1 and predictor
                                        % these_training_measurements107. The data has zero variance.
                                        try
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                        catch
                                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                                        end
                                    case 4
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                                    case 5
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                    case 6
                                        handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                                end
                            end

                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                            %Note: for some reason this did not work for net when I did this:
                            % [label_traces,score] = predict(Mdl,traces');
                            % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                            %I had to resort to the for loop:


                            if MLalgo==6
                                for ii=1:size(traces_trimmed,2)
                                    this_time_point=zeros(1,size(traces_trimmed,1));
                                    this_time_point(1,:)=traces_trimmed(:,ii);
                                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                    scores_sh2(ii_sh,ii,:)=score;
                                    if label>0.5
                                        label_traces_sh2(ii_sh,ii)=1;
                                    else
                                        label_traces_sh2(ii_sh,ii)=0;
                                    end
                                end


                                for ii=1:size(measurements_post_trimmed,1)
                                    this_time_point=zeros(1,size(traces_trimmed,1));
                                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                                    scores_post_sh2(ii_sh,ii,:)=score;
                                    if label>0.5
                                        label_post_sh2(ii_sh,ii)=1;
                                    else
                                        label_post_sh2(ii_sh,ii)=0;
                                    end
                                end
                            else

                                for ii=1:size(traces_trimmed,2)
                                    this_time_point=zeros(1,size(traces_trimmed,1));
                                    this_time_point(1,:)=traces_trimmed(:,ii);
                                    [label_traces_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                    scores_sh2(ii_sh,ii,:)=score;
                                end

                                for ii=1:size(measurements_post_trimmed,1)
                                    this_time_point=zeros(1,size(traces_trimmed,1));
                                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                                    [label_post_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                                    scores_post_sh2(ii_sh,ii,:)=score;
                                end
                            end

                            handles_out.MLalgo(MLalgo).label_traces_sh2=label_traces_sh2;
                            handles_not_out.MLalgo(MLalgo).label_post_sh2=label_post_sh2;
                            handles_out.MLalgo(MLalgo).scores_sh2=scores_sh2;
                            handles_not_out.MLalgo(MLalgo).scores_post_sh2=scores_post_sh2;

                        end


                    end

                case 4
                    %leave one trial out
                    %Train each time point separately
                    gcp

                    handles_not_out.MLalgo(MLalgo).models=[];
                    handles_not_out.MLalgo(MLalgo).processed_succesfully=1;

                    sp_trials=size(dFF_per_trial_sp,1);
                    sm_trials=size(dFF_per_trial_sm,1);
                    no_trials=sp_trials+sm_trials;

                    no_time_points=size(dFF_per_trial_sm,3);
                    this_no_ROIs=sum(p_value_mask);
                    dFF_per_trial_all=zeros(sm_trials+sp_trials,this_no_ROIs,no_time_points);
                    dFF_per_trial_all(1:sm_trials,:,:)=dFF_per_trial_sm(:,p_value_mask,:);
                    dFF_per_trial_all(sm_trials+1:end,:,:)=dFF_per_trial_sp(:,p_value_mask,:);
                    decisions_per_trial4=zeros(1,sp_trials+sm_trials);
                    decisions_per_trial4(1,sm_trials+1:end)=1;

                    label_traces=zeros(no_trials,no_time_points);
                    scores=zeros(no_trials,no_time_points,2);

                    time_span=[0:dt:dt*no_time_points]-dt_span+dt;
                    time_span=time_span(1:end-1)+post_shift;

                    %Do decoding per time point
                    for ii_t=1:no_time_points
                        parfor kk=1:no_trials


                            training_mask=ones(size(dFF_per_trial_all,1),1);
                            training_mask(kk,1)=0;

                            these_training_measurements=zeros(sum(training_mask),this_no_ROIs);
                            these_training_decisions=zeros(1,sum(training_mask));


                            jj=0;
                            for ii=1:no_trials
                                if training_mask(ii)==1
                                    jj=jj+1;
                                    these_training_measurements(jj,:)=dFF_per_trial_all(ii,:,ii_t);
                                    these_training_decisions(jj)=decisions_per_trial4(ii);
                                end
                            end

                            tblTrn=[];
                            tblTrn = array2table(these_training_measurements);

                            %Store the decisions in Y
                            Y=these_training_decisions;

                            switch MLalgo
                                case 1
                                    try
                                        Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        % Error using ClassificationDiscriminant (line 380)
                                        % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                        % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                        Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                    end
                                case 2
                                    Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                case 3
                                    %The try catch was entered here because of this
                                    %error
                                    % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                    % A normal distribution cannot be fit for the combination of class 1 and predictor
                                    % these_training_measurements107. The data has zero variance.
                                    try
                                        Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        %                                     handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                                    end
                                case 4
                                    Mdl = fitcnet(tblTrn,Y);
                                case 5
                                    Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                case 6
                                    Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                            end

                            these_predict_measurements=zeros(1,this_no_ROIs);
                            these_predict_measurements(1,:)=dFF_per_trial_all(kk,:,ii_t);




                            [label,score] = predict(Mdl,these_predict_measurements);
                            scores(kk,ii_t,:)=score;
                            if label>0.5
                                label_traces(kk,ii_t)=1;
                            else
                                label_traces(kk,ii_t)=0;
                            end


                        end
                        %                     if time_span(ii_t)>-6.5
                        %                         pffft=1;
                        %                     end
                    end

                    %Now do shuffled decoding
                    label_traces_sh=zeros(no_shuffles,no_trials,no_time_points);
                    scores_sh=zeros(no_shuffles,no_trials,no_time_points,2);
                    randperms=[];

                    for ii_sh=1:no_shuffles

                        if ii_sh==1
                            this_randperm=randperm(no_trials);
                            randperms(1,:)=this_randperm;
                        else
                            not_found=1;
                            while not_found==1
                                for rp_ii=1:size(randperms,1)
                                    this_randperm=randperm(no_trials);
                                    if sum(randperms(rp_ii,:)==this_randperm)~=no_trials
                                        rp_ii=rp_ii+1;
                                        randperms(rp_ii,:)=this_randperm;
                                        not_found=0;
                                    end
                                end
                            end
                        end

                        for ii_t=1:no_time_points
                            parfor kk=1:no_trials


                                training_mask=ones(size(dFF_per_trial_all,1),1);
                                training_mask(kk,1)=0;

                                these_training_measurements=zeros(sum(training_mask),this_no_ROIs);
                                these_training_decisions=zeros(1,sum(training_mask));

                                shuffled_decisions=decisions_per_trial4(this_randperm);

                                jj=0;
                                for ii=1:no_trials
                                    if training_mask(ii)==1
                                        jj=jj+1;
                                        these_training_measurements(jj,:)=dFF_per_trial_all(ii,:,ii_t);
                                        these_training_decisions(jj)=shuffled_decisions(ii);
                                    end
                                end

                                tblTrn=[];
                                tblTrn = array2table(these_training_measurements);

                                %Store the decisions in Y
                                Y=these_training_decisions;

                                switch MLalgo
                                    case 1
                                        try
                                            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                        catch
                                            % Error using ClassificationDiscriminant (line 380)
                                            % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                            % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                                        end
                                    case 2
                                        Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                    case 3
                                        %The try catch was entered here because of this
                                        %error
                                        % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                        % A normal distribution cannot be fit for the combination of class 1 and predictor
                                        % these_training_measurements107. The data has zero variance.
                                        try
                                            Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                        catch
                                            %                                         handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                                        end
                                    case 4
                                        Mdl = fitcnet(tblTrn,Y);
                                    case 5
                                        Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                    case 6
                                        Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                                end

                                these_predict_measurements=zeros(1,this_no_ROIs);
                                these_predict_measurements(1,:)=dFF_per_trial_all(kk,:,ii_t);




                                [label,score] = predict(Mdl,these_predict_measurements);
                                scores_sh(ii_sh,kk,ii_t,:)=score;
                                if label>0.5
                                    label_traces_sh(ii_sh,kk,ii_t)=1;
                                else
                                    label_traces_sh(ii_sh,kk,ii_t)=0;
                                end





                            end
                        end
                    end

            end

            switch processing_algorithm
                case 3
                    if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
                        %label is the predicted label, and score is the predicted class
                        %posterior probability
                        for ii=1:length(training_decisions_post)
                            if label_post(ii)==training_decisions_post(ii)
                                correct_predict_tr(ii)=1;
                            else
                                correct_predict_tr(ii)=0;
                            end
                        end


                        %Calculate wta for windows of no_points_post
                        correct_predict_tr_wta=zeros(1,length(training_decisions_post));
                        for ii=1:length(training_decisions_post)-no_points_post
                            this_correct=zeros(1,no_points_post);
                            for jj=1:no_points_post
                                if label_post(ii+jj-1)==training_decisions_post(ii+jj-1)
                                    this_correct(jj)=1;
                                end
                            end
                            if sum(this_correct)>(no_points_post/2)
                                correct_predict_tr_wta(ii+floor(no_points_post/2))=1;
                            else
                                correct_predict_tr_wta(ii+floor(no_points_post/2))=0;
                            end
                        end

                        %Now do shuffled
                        %label is the predicted label, and score is the predicted class
                        %posterior probability
                        for ii=1:length(training_decisions_post_sh)
                            if label_post_sh(ii)==training_decisions_post_sh(ii)
                                correct_predict_tr_sh(ii)=1;
                            else
                                correct_predict_tr_sh(ii)=0;
                            end
                        end

                        %Calculate wta for windows of no_points_post
                        correct_predict_tr_wta_sh=zeros(1,length(training_decisions_post_sh));
                        for ii=1:length(training_decisions_post_sh)-no_points_post
                            this_correct=zeros(1,no_points_post);
                            for jj=1:no_points_post
                                if label_post_sh(ii+jj-1)==training_decisions_post_sh(ii+jj-1)
                                    this_correct(jj)=1;
                                end
                            end
                            if sum(this_correct)>(no_points_post/2)
                                correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=1;
                            else
                                correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=0;
                            end
                        end

                        %Now do shuffled per trial
                        correct_predict_tr_wta_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
                        correct_predict_tr_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
                        for ii_sh=1:no_shuffles
                            %label is the predicted label, and score is the predicted class
                            %posterior probability
                            for ii=1:size(training_decisions_post_sh2,2)
                                if label_post_sh2(ii_sh,ii)==training_decisions_post_sh2(ii_sh,ii)
                                    correct_predict_tr_sh2(ii_sh,ii)=1;
                                else
                                    correct_predict_tr_sh2(ii_sh,ii)=0;
                                end
                            end

                            %Calculate wta for windows of no_points_post

                            for ii=1:size(training_decisions_post_sh2,2)-no_points_post
                                this_correct=zeros(1,no_points_post);
                                for jj=1:no_points_post
                                    if label_post_sh2(ii_sh,ii+jj-1)==training_decisions_post_sh2(ii_sh,ii+jj-1)
                                        this_correct(jj)=1;
                                    end
                                end
                                if sum(this_correct)>(no_points_post/2)
                                    correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=1;
                                else
                                    correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=0;
                                end
                            end

                        end

                        if show_figures==1
                            fprintf(1, ['Training accuracy for ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
                                ,ii_cost, sum(correct_predict_tr)/length(correct_predict_tr),sum(correct_predict_tr_wta)/length(correct_predict_tr_wta));
                            fprintf(1, ['Shuffled training accuracy for ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
                                ,ii_cost, sum(correct_predict_tr_sh)/length(correct_predict_tr_sh),sum(correct_predict_tr_wta_sh)/length(correct_predict_tr_wta_sh));
                            fprintf(1, ['Training accuracy for shuffled per trial ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
                                ,ii_cost, sum(correct_predict_tr_sh2(:))/length(correct_predict_tr_sh2(:)),sum(correct_predict_tr_wta_sh2(:))/length(correct_predict_tr_wta_sh2(:)));
                            fprintf(1, ['Mean label trace %d, variance %d\n'],mean(label_traces),var(label_traces));
                        end

                        handles_not_out.MLalgo(MLalgo).correct_predict_tr=correct_predict_tr;
                        handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta=correct_predict_tr_wta;
                        handles_out2.ROI(iiROI).MLalgo(MLalgo).accuracy_tr=sum(correct_predict_tr)/length(correct_predict_tr);
                        handles_out2.ROI(iiROI).MLalgo(MLalgo).accuracy_tr_wta=sum(correct_predict_tr_wta)/length(correct_predict_tr_wta);

                        handles_not_out.MLalgo(MLalgo).correct_predict_tr_sh=correct_predict_tr_sh;
                        handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta_sh=correct_predict_tr_wta_sh;
                        handles_out.MLalgo(MLalgo).accuracy_tr_sh=sum(correct_predict_tr_sh)/length(correct_predict_tr_sh);
                        handles_out.MLalgo(MLalgo).accuracy_tr_wta_sh=sum(correct_predict_tr_wta_sh)/length(correct_predict_tr_wta_sh);

                        handles_not_out.MLalgo(MLalgo).correct_predict_tr_sh2=correct_predict_tr_sh2;
                        handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta_sh2=correct_predict_tr_wta_sh2;
                        handles_out.MLalgo(MLalgo).accuracy_tr_sh2=sum(correct_predict_tr_sh2(:))/length(correct_predict_tr_sh2(:));
                        handles_out.MLalgo(MLalgo).accuracy_tr_wta_sh2=sum(correct_predict_tr_wta_sh2(:))/length(correct_predict_tr_wta_sh2(:));

                        handles_not_out.MLalgo(MLalgo).mean_label_traces=mean(label_traces);
                        handles_not_out.MLalgo(MLalgo).var_label_traces=var(label_traces);

                        moving_mean_label_traces = movmean(label_traces,moving_mean_n);
                        handles_not_out.MLalgo(MLalgo).label_traces=label_traces;

                        %Now let's do the carpentry
                        %             moving_mean_label_traces_sh = movmean(mean(label_traces_sh2),moving_mean_n);

                        moving_mean_label_traces_sh2 = movmean(label_traces_sh2,moving_mean_n);
                        handles_not_out.MLalgo(MLalgo).label_traces_sh2=label_traces_sh2;
                        moving_mean_label_traces_sh = movmean(label_traces_sh,moving_mean_n);
                        handles_not_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;

                        %             if show_figures==1
                        %                 figNo=figNo+1;
                        %                 try
                        %                     close(figNo)
                        %                 catch
                        %                 end
                        %
                        %                 hFig = figure(figNo);
                        %
                        %                 set(hFig, 'units','normalized','position',[.05 .1 .85 .3])
                        %
                        %                 hold on
                        %
                        %
                        %                 %                 CIsm = bootci(1000, @mean, moving_mean_label_traces_sh2);
                        %                 %                 meansm=mean(moving_mean_label_traces_sh2,1);
                        %                 %                 CIsm(1,:)=meansm-CIsm(1,:);
                        %                 %                 CIsm(2,:)=CIsm(2,:)-meansm;
                        %                 %
                        %                 %                 %S- Proficient
                        %                 %                 [hlsm, hpsm] = boundedline(time',mean(moving_mean_label_traces_sh2,1)', CIsm', 'cmap',[80/255 194/255 255/255]);
                        %                 %
                        %
                        %                 per95=prctile(moving_mean_label_traces_sh2(:),95);
                        %                 per5=prctile(moving_mean_label_traces_sh2(:),5);
                        %                 CIsh=[mean(moving_mean_label_traces_sh2(:))-per5 per95-mean(moving_mean_label_traces_sh2(:))]';
                        %                 [hlCR, hpCR] = boundedline([time(1) time(end)],[mean(moving_mean_label_traces_sh2(:)) mean(moving_mean_label_traces_sh2(:))], CIsh', 'cmap',[80/255 194/255 255/255]);
                        %
                        %
                        %                 %                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
                        %                 plot(time,moving_mean_label_traces,'-k','LineWidth',1)
                        %                 plot(time,1.1*((epochs==8)+(epochs==9)),'-b')
                        %                 plot(time,1.1*((epochs==6)+(epochs==7)),'-r')
                        %
                        %                 %                 plot(time,moving_mean_label_traces_sh(1,:),'-b')
                        %
                        %                 ylim([-0.2 1.2])
                        %                 title(['Label prediction for entire session for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
                        %
                        %             end

                        handles_out2.time=time;

                        %Calculate Shanon enthropy
                        pone=sum(label_traces==1)/length(label_traces);
                        pzero=sum(label_traces==0)/length(label_traces);
                        handles_out.MLalgo(MLalgo).shannon_e=-pzero*log2(pzero) - pone*log2(pone);

                        pone=sum(label_traces_sh==1)/length(label_traces_sh);
                        pzero=sum(label_traces_sh==0)/length(label_traces_sh);
                        handles_out.MLalgo(MLalgo).shannon_e_sh=-pzero*log2(pzero) - pone*log2(pone);

                        for ii=1:no_shuffles
                            pone=sum(label_traces_sh2(ii,:)==1)/size(label_traces_sh2,2);
                            pzero=sum(label_traces_sh2(ii,:)==0)/size(label_traces_sh2,2);
                            handles_out.MLalgo(MLalgo).shannon_e_sh2(ii)=-pzero*log2(pzero) - pone*log2(pone);
                        end
                        if show_figures==1
                            fprintf(1, ['Shannon entropy %d, shuffled 1 %d, shuffled 2 %d\n\n'],handles_out.MLalgo(MLalgo).shannon_e...
                                ,handles_out.MLalgo(MLalgo).shannon_e_sh,mean(handles_out.MLalgo(MLalgo).shannon_e_sh2));
                        end

                        %Now let's do accounting and show it in a bar graph

                        %post Splus
                        post_label_sp=label_traces(logical(epochs_sp_post));
                        points_per_cut=no_points_post;
                        no_cuts=floor(length(post_label_sp)/points_per_cut);
                        mean_post_label_sp=zeros(1,no_cuts);
                        for ii=1:no_cuts
                            mean_post_label_sp(ii)=mean(post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        end
                        %     sh_post_label_sp=sum(sh_label_traces,1)/size(sh_label_traces,1);
                        %     sh_post_label_sp=sh_post_label_sp(logical(epochs_sp_post));
                        %     mean_sh_post_label_sp=zeros(1,no_cuts);
                        %     for ii=1:no_cuts
                        %         mean_sh_post_label_sp(ii)=mean(sh_post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        %     end

                        %post Sminus
                        post_label_sm=label_traces(logical(epochs_sm_post));
                        points_per_cut=no_points_post;
                        no_cuts=floor(length(post_label_sm)/points_per_cut);
                        mean_post_label_sm=zeros(1,no_cuts);
                        for ii=1:no_cuts
                            mean_post_label_sm(ii)=mean(post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        end
                        %     sh_post_label_sm=sum(sh_label_traces,1)/size(sh_label_traces,1);
                        %     sh_post_label_sm=sh_post_label_sm(logical(epochs_sm_post));
                        %     mean_sh_post_label_sm=zeros(1,no_cuts);
                        %     for ii=1:no_cuts
                        %         mean_sh_post_label_sm(ii)=mean(sh_post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        %     end


                        %pre Sminus
                        pre_label_sm=label_traces(logical(epochs_sm_pre));
                        points_per_cut=no_points_pre;
                        no_cuts=floor(length(pre_label_sm)/points_per_cut);
                        mean_pre_label_sm=zeros(1,no_cuts);
                        for ii=1:no_cuts
                            mean_pre_label_sm(ii)=mean(pre_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        end


                        %pre Splus
                        pre_label_sp=label_traces(logical(epochs_sp_pre));
                        points_per_cut=no_points_pre;
                        no_cuts=floor(length(pre_label_sp)/points_per_cut);
                        mean_pre_label_sp=zeros(1,no_cuts);
                        for ii=1:no_cuts
                            mean_pre_label_sp(ii)=mean(pre_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        end

                        %     sh_pre_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
                        %     sh_pre_label=sh_pre_label(logical(epochs_sm_pre+epochs_sp_pre));
                        %     mean_sh_pre_label=zeros(1,no_cuts);
                        %     for ii=1:no_cuts
                        %         mean_sh_pre_label(ii)=mean(sh_pre_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        %     end

                        %all
                        all_label=label_traces;
                        %points_per_cut=no_points_pre;
                        points_per_cut=no_points_post;
                        no_cuts=floor(length(all_label)/points_per_cut);
                        mean_all_label=zeros(1,no_cuts);
                        for ii=1:no_cuts
                            mean_all_label(ii)=mean(all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        end
                        %     sh_all_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
                        %     mean_sh_all_label=zeros(1,no_cuts);
                        %     for ii=1:no_cuts
                        %         mean_sh_all_label(ii)=mean(sh_all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
                        %     end

                        %             if show_figures==1
                        %                 %Note that pre and post are refrenced to the start of the training period
                        %                 edges=[0:0.033:1.2];
                        %                 rand_offset=0.8;
                        %
                        %
                        %                 figNo=figNo+1;
                        %                 try
                        %                     close(figNo)
                        %                 catch
                        %                 end
                        %
                        %                 hFig = figure(figNo);
                        %
                        %                 set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
                        %
                        %                 hold on
                        %
                        %                 bar_offset=1;
                        %
                        %                 %                 %S- pre
                        %                 %                 bar_offset=1;
                        %                 %
                        %                 %                 bar(bar_offset,mean(mean_pre_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                        %                 %
                        %                 %                 %Violin plot
                        %                 %                 [mean_out, CIout]=drgViolinPoint(mean_pre_label_sm...
                        %                 %                     ,edges,bar_offset,rand_offset,'k','k',3);
                        %                 %
                        %                 %                 bar_offset=bar_offset+1;
                        %                 %
                        %                 %                 %S+ pre
                        %                 %                 bar(bar_offset,mean(mean_pre_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                        %                 %
                        %                 %                 %Violin plot
                        %                 %                 [mean_out, CIout]=drgViolinPoint(mean_pre_label_sp...
                        %                 %                     ,edges,bar_offset,rand_offset,'k','k',3);
                        %                 %
                        %                 %                 bar_offset=bar_offset+2;
                        %
                        %                 %S- post
                        %                 bar(bar_offset,mean(mean_post_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                        %
                        %                 %Violin plot
                        %                 [mean_out, CIout]=drgViolinPoint(mean_post_label_sm...
                        %                     ,edges,bar_offset,rand_offset,'k','k',3);
                        %                 bar_offset=bar_offset+1;
                        %
                        %                 %S+ post
                        %                 bar(bar_offset,mean(mean_post_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                        %
                        %                 %Violin plot
                        %                 [mean_out, CIout]=drgViolinPoint(mean_post_label_sp...
                        %                     ,edges,bar_offset,rand_offset,'k','k',3);
                        %
                        %                 bar_offset=bar_offset+2;
                        %                 bar(bar_offset,mean(mean_all_label),'LineWidth', 3,'EdgeColor','none','FaceColor','m')
                        %
                        %                 %Violin plot
                        %                 [mean_out, CIout]=drgViolinPoint(mean_all_label...
                        %                     ,edges,bar_offset,rand_offset,'k','k',3);
                        %
                        %                 xticks([1 2 4 5 7])
                        %                 xticklabels({'S- post', 'S+ post','Entire session'})
                        %                 title(['Label prediction for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
                        %
                        %             end

                        %Now show average labels for S+, S-, etc for 30 sec
                        handles_not_out.MLalgo(MLalgo).mean_all_label=mean_all_label;
                        handles_not_out.MLalgo(MLalgo).mean_pre_label_sm=mean_pre_label_sm;
                        handles_not_out.MLalgo(MLalgo).mean_pre_label_sp=mean_pre_label_sp;
                        handles_not_out.MLalgo(MLalgo).mean_post_label_sm=mean_post_label_sm;
                        handles_not_out.MLalgo(MLalgo).mean_post_label_sp=mean_post_label_sp;
                        handles_not_out.MLalgo(MLalgo).mean_all_label=mean_all_label;

                        at_end=0;
                        ii=1;
                        sp_ii=0;
                        per_trial_sp_timecourse=[];
                        epoch_before_sp=[];
                        sm_ii=0;
                        per_trial_sm_timecourse=[];
                        epoch_before_sm=[];
                        per_trial_scores_sp=[];
                        per_trial_scores_sm=[];

                        last_sp_sm=-1;


                        while at_end==0
                            next_ii=[];
                            next_ii_sp=find(epochs_sp_post(ii:end)==1,1,'first');
                            next_ii_sm=find(epochs_sm_post(ii:end)==1,1,'first');
                            if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
                                if next_ii_sp<next_ii_sm
                                    next_ii=next_ii_sp;
                                    this_sp_sm=1;
                                else
                                    next_ii=next_ii_sm;
                                    this_sp_sm=0;
                                end
                            else
                                if ~isempty(next_ii_sp)
                                    next_ii=next_ii_sp;
                                    this_sp_sm=1;
                                end
                                if ~isempty(next_ii_sm)
                                    next_ii=next_ii_sm;
                                    this_sp_sm=0;
                                end
                            end

                            if ~isempty(next_ii)
                                if ((ii+next_ii-ii_span)>0)&((ii+next_ii+ii_span<length(label_traces)))
                                    if this_sp_sm==1
                                        sp_ii=sp_ii+1;
                                        per_trial_sp_timecourse(sp_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                                        per_trial_scores_sp(sp_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
                                        epoch_before_sp(sp_ii)=last_sp_sm;
                                        last_sp_sm=1;
                                        ii_next_post=find(epochs_sp_post(ii+next_ii:end)==0,1,'first');
                                        ii=ii+next_ii+ii_next_post;
                                    else
                                        sm_ii=sm_ii+1;
                                        per_trial_sm_timecourse(sm_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                                        per_trial_scores_sm(sm_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
                                        epoch_before_sm(sm_ii)=last_sp_sm;
                                        last_sp_sm=0;
                                        ii_next_post=find(epochs_sm_post(ii+next_ii:end)==0,1,'first');
                                        ii=ii+next_ii+ii_next_post;
                                    end
                                else
                                    if  ((ii+next_ii+ii_span>length(label_traces)))
                                        at_end=1;
                                    else
                                        ii=ii+next_ii;
                                    end
                                end
                            else
                                at_end=1;
                            end

                        end

                        handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse=per_trial_sp_timecourse;
                        handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse=per_trial_sm_timecourse;



                        this_moving_mean_n=10;
                        moving_mean_per_trial_sp_timecourse = movmean(per_trial_sp_timecourse',this_moving_mean_n)';
                        moving_mean_per_trial_sm_timecourse = movmean(per_trial_sm_timecourse',this_moving_mean_n)';

                        handles_out2.ROI(iiROI).MLalgo(MLalgo).moving_mean_per_trial_sp_timecourse=moving_mean_per_trial_sp_timecourse;
                        handles_out2.ROI(iiROI).MLalgo(MLalgo).moving_mean_per_trial_sp_timecourse=moving_mean_per_trial_sp_timecourse;

                        time_span=[0:dt:dt*size(per_trial_sp_timecourse,2)]-dt_span+dt;
                        time_span=time_span(1:end-1)+post_shift;

                        %Calculate correct predict
                        for ii_tr=1:size(per_trial_sp_timecourse,1)
                            for ii_time=1:size(per_trial_sp_timecourse,2)
                                if per_trial_sp_timecourse(ii_tr,ii_time)==1
                                    this_correct_predict(ii_tr,ii_time)=1;
                                else
                                    this_correct_predict(ii_tr,ii_time)=0;
                                end
                            end
                        end

                        for ii_tr=1:size(per_trial_sm_timecourse,1)
                            for ii_time=1:size(per_trial_sm_timecourse,2)
                                if per_trial_sm_timecourse(ii_tr,ii_time)==0
                                    this_correct_predict(ii_tr+sp_ii,ii_time)=1;
                                else
                                    this_correct_predict(ii_tr+sp_ii,ii_time)=0;
                                end
                            end
                        end


                        %Calculate correct predict shuffled
                        ii_plus=0;
                        for ww=1:10
                            for ii_tr=1:size(per_trial_sp_timecourse,1)
                                rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                                for ii_time=1:size(per_trial_sp_timecourse,2)
                                    if per_trial_sp_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                                    else
                                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                                    end
                                end
                            end

                            ii_plus=ii_plus+sp_ii;

                            for ii_tr=1:size(per_trial_sm_timecourse,1)
                                rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                                for ii_time=1:size(per_trial_sm_timecourse,2)
                                    if per_trial_sm_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                                    else
                                        this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                                    end
                                end
                            end
                            ii_plus=ii_plus+sm_ii;
                        end

                        handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict=this_correct_predict;
                        handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh=this_correct_predict_sh;

                        %             fprintf(1,'Accuracy for ROI No %d is %d\n', iiROI,mean(this_correct_predict(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2)))))
                        %             fprintf(1,'Shuffled trial accuracy for ROI No %d is %d\n', iiROI,mean(this_correct_predict_sh(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2)))))

                        %Calculate correct predict for the odor window


                        %             if show_figures==1
                        %                 %Plot the prediction for S+ and S- and the correct predictions
                        %                 figNo=figNo+1;
                        %                 try
                        %                     close(figNo)
                        %                 catch
                        %                 end
                        %
                        %                 hFig = figure(figNo);
                        %
                        %                 set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
                        %
                        %                 subplot(2,1,1)
                        %                 hold on
                        %
                        %                 CIsm = bootci(1000, @mean, moving_mean_per_trial_sm_timecourse);
                        %                 meansm=mean(moving_mean_per_trial_sm_timecourse,1);
                        %                 CIsm(1,:)=meansm-CIsm(1,:);
                        %                 CIsm(2,:)=CIsm(2,:)-meansm;
                        %
                        %                 [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sm_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
                        %
                        %                 CIsp = bootci(1000, @mean, moving_mean_per_trial_sp_timecourse);
                        %                 meansp=mean(moving_mean_per_trial_sp_timecourse,1);
                        %                 CIsp(1,:)=meansp-CIsp(1,:);
                        %                 CIsp(2,:)=CIsp(2,:)-meansp;
                        %
                        %
                        %                 [hlsp, hpsp] = boundedline(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', CIsp', 'cmap',[0 114/255 178/255]);
                        %
                        %                 plot(time_span',mean(moving_mean_per_trial_sm_timecourse,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                        %                 plot(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
                        %
                        %                 text(30,0.75,'S-','Color',[158/255 31/255 99/255])
                        %                 text(30,0.65,'S+','Color',[0 114/255 178/255])
                        %
                        %                 ylim([0 1])
                        %                 title(['Label prediction per trial for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
                        %                 xlabel('Time(sec)')
                        %                 ylabel('Label prediction, S+=1, S-=0')
                        %
                        %                 subplot(2,1,2)
                        %                 hold on
                        %
                        %                 CIsm = bootci(1000, @mean, this_correct_predict_sh);
                        %                 meansm=mean(this_correct_predict_sh,1);
                        %                 CIsm(1,:)=meansm-CIsm(1,:);
                        %                 CIsm(2,:)=CIsm(2,:)-meansm;
                        %
                        %                 [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');
                        %
                        %
                        %                 CIsm = bootci(1000, @mean, this_correct_predict);
                        %                 meansm=mean(this_correct_predict,1);
                        %                 CIsm(1,:)=meansm-CIsm(1,:);
                        %                 CIsm(2,:)=CIsm(2,:)-meansm;
                        %
                        %                 [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);
                        %
                        %                 plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
                        %                 plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);
                        %
                        %                 text(30,0.75,'Shuffled','Color','k')
                        %                 text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])
                        %
                        %                 ylim([0 1])
                        %                 title(['Percent correct accuracy per trial for ' classifier_names{MLalgo} ' and cost ' num2str(ii_cost)])
                        %                 xlabel('Time(sec)')
                        %                 ylabel('Accuracy')
                        %
                        %                 %Plot the posterior probabilities for Sp (scores,:,2) and Sm (scores(:,1))
                        %                 figNo=figNo+1;
                        %                 try
                        %                     close(figNo)
                        %                 catch
                        %                 end
                        %
                        %                 hFig = figure(figNo);
                        %
                        %                 set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
                        %
                        %                 subplot(2,1,1)
                        %                 hold on
                        %
                        %                 these_scores_sm=zeros(size(per_trial_scores_sp,1),size(per_trial_scores_sp,3));
                        %                 these_scores_sm(:,:)=per_trial_scores_sp(:,1,:);
                        %                 CIsm = bootci(1000, @mean, these_scores_sm);
                        %                 meansm=mean(these_scores_sm,1);
                        %                 CIsm(1,:)=meansm-CIsm(1,:);
                        %                 CIsm(2,:)=CIsm(2,:)-meansm;
                        %
                        %                 [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
                        %
                        %                 these_scores_sp=zeros(size(per_trial_scores_sp,1),size(per_trial_scores_sp,3));
                        %                 these_scores_sp(:,:)=per_trial_scores_sp(:,2,:);
                        %                 CIsp = bootci(1000, @mean, these_scores_sp);
                        %                 meansp=mean(these_scores_sp,1);
                        %                 CIsp(1,:)=meansp-CIsp(1,:);
                        %                 CIsp(2,:)=CIsp(2,:)-meansp;
                        %
                        %
                        %                 [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);
                        %
                        %                 plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                        %                 plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
                        %
                        %                 text(30,0.75,'S-','Color',[158/255 31/255 99/255])
                        %                 text(30,0.65,'S+','Color',[0 114/255 178/255])
                        %
                        %                 ylim([0 1])
                        %                 title(['Posterior probability for S+ or S- prediction for S+ trials ' classifier_names{MLalgo} ' ' num2str(p_threshold)])
                        %                 xlabel('Time(sec)')
                        %                 ylabel('Posterior probability')
                        %
                        %                 subplot(2,1,2)
                        %                 hold on
                        %
                        %                 these_scores_sm=zeros(size(per_trial_scores_sm,1),size(per_trial_scores_sm,3));
                        %                 these_scores_sm(:,:)=per_trial_scores_sm(:,1,:);
                        %                 CIsm = bootci(1000, @mean, these_scores_sm);
                        %                 meansm=mean(these_scores_sm,1);
                        %                 CIsm(1,:)=meansm-CIsm(1,:);
                        %                 CIsm(2,:)=CIsm(2,:)-meansm;
                        %
                        %                 [hlsm, hpsm] = boundedline(time_span',mean(these_scores_sm,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
                        %
                        %                 these_scores_sp=zeros(size(per_trial_scores_sm,1),size(per_trial_scores_sm,3));
                        %                 these_scores_sp(:,:)=per_trial_scores_sm(:,2,:);
                        %                 CIsp = bootci(1000, @mean, these_scores_sp);
                        %                 meansp=mean(these_scores_sp,1);
                        %                 CIsp(1,:)=meansp-CIsp(1,:);
                        %                 CIsp(2,:)=CIsp(2,:)-meansp;
                        %
                        %
                        %                 [hlsp, hpsp] = boundedline(time_span',mean(these_scores_sp,1)', CIsp', 'cmap',[0 114/255 178/255]);
                        %
                        %                 plot(time_span',mean(these_scores_sm,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                        %                 plot(time_span',mean(these_scores_sp,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
                        %
                        %                 text(30,0.75,'S-','Color',[158/255 31/255 99/255])
                        %                 text(30,0.65,'S+','Color',[0 114/255 178/255])
                        %
                        %                 ylim([0 1])
                        %                 title(['Posterior probability for S+ or S- prediction for S- trials ' classifier_names{MLalgo} ' ' num2str(p_threshold)])
                        %                 xlabel('Time(sec)')
                        %                 ylabel('Posterior probability')
                        %
                        %             end
                    else
                        fprintf(1, [classifier_names{MLalgo} ' was not processed succesfully\n']);
                    end

                case 4

                    %get the per trial sp and sm timecourse for predicted
                    %labels

                    per_trial_sp_timecourse=zeros(sum(decisions_per_trial4==1),size(label_traces,2));
                    per_trial_sp_timecourse_sh=zeros(sum(decisions_per_trial4==1),size(label_traces,2));
                    ii_sp=0;

                    per_trial_sm_timecourse=zeros(sum(decisions_per_trial4==0),size(label_traces,2));
                    per_trial_sm_timecourse_sh=zeros(sum(decisions_per_trial4==0),size(label_traces,2));
                    ii_sm=0;


                    for trNo=1:no_trials
                        if  decisions_per_trial4(trNo)==1
                            ii_sp=ii_sp+1;
                            per_trial_sp_timecourse(ii_sp,:)=label_traces(trNo,:);

                            this_mean_label_trace=zeros(1,size(label_traces,2));
                            this_mean_label_trace(:,:)=mean(label_traces_sh(:,trNo,:),1);
                            per_trial_sp_timecourse_sh(ii_sp,:)=this_mean_label_trace;
                        else
                            ii_sm=ii_sm+1;
                            per_trial_sm_timecourse(ii_sm,:)=label_traces(trNo,:);

                            this_mean_label_trace=zeros(1,size(label_traces,2));
                            this_mean_label_trace(:,:)=mean(label_traces_sh(:,trNo,:),1);
                            per_trial_sm_timecourse_sh(ii_sm,:)=this_mean_label_trace;
                        end
                    end

                    handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse=per_trial_sp_timecourse;
                    handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse=per_trial_sm_timecourse;

                    time_span=[0:dt:dt*size(per_trial_sp_timecourse,2)]-dt_span+dt;
                    time_span=time_span(1:end-1)+post_shift;

                    this_correct_predict=zeros(no_trials,size(per_trial_sp_timecourse,2));
                    this_correct_predict_sh=zeros(no_trials,size(per_trial_sp_timecourse,2));

                    %Calculate correct predict
                    for ii_tr=1:size(per_trial_sp_timecourse,1)
                        for ii_time=1:size(per_trial_sp_timecourse,2)
                            if per_trial_sp_timecourse(ii_tr,ii_time)==1
                                this_correct_predict(ii_tr,ii_time)=1;
                            else
                                this_correct_predict(ii_tr,ii_time)=0;
                            end
                        end
                    end

                    sp_ii=size(per_trial_sp_timecourse,1);
                    for ii_tr=1:size(per_trial_sm_timecourse,1)
                        for ii_time=1:size(per_trial_sm_timecourse,2)
                            if per_trial_sm_timecourse(ii_tr,ii_time)==0
                                this_correct_predict(ii_tr+sp_ii,ii_time)=1;
                            else
                                this_correct_predict(ii_tr+sp_ii,ii_time)=0;
                            end
                        end
                    end


                    %Calculate correct predict shuffled
                    sm_ii=size(per_trial_sm_timecourse,1);
                    this_correct_predict_sh=[];
                    ii_plus=0;
                    for ww=1:10
                        for ii_tr=1:size(per_trial_sp_timecourse,1)
                            rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                            for ii_time=1:size(per_trial_sp_timecourse,2)
                                if per_trial_sp_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                                    this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                                else
                                    this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                                end
                            end
                        end

                        ii_plus=ii_plus+sp_ii;

                        for ii_tr=1:size(per_trial_sm_timecourse,1)
                            rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                            for ii_time=1:size(per_trial_sm_timecourse,2)
                                if per_trial_sm_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                                    this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                                else
                                    this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                                end
                            end
                        end
                        ii_plus=ii_plus+sm_ii;
                    end

                    handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict=this_correct_predict;
                    handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh=this_correct_predict_sh;

            end


        end

        if show_figures==1
            %Plot the accuracy
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);

            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

            hold on

            CIsm = bootci(1000, @mean, this_correct_predict_sh);
            meansm=mean(this_correct_predict_sh,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');


            CIsm = bootci(1000, @mean, this_correct_predict);
            meansm=mean(this_correct_predict,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);

            plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
            plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);

            text(30,0.75,'Shuffled','Color','k')
            text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])


            ylim([-0.1 1.1])
            xlim([-7 15])
            this_ylim=ylim;

            %FV
            rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
            plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

            %Odor on markers
            plot([0 0],this_ylim,'-k')
            odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

            title(['Accuracy for ii_ROI '  num2str(iiROI)])
            xlabel('Time(sec)')
            ylabel('Accuracy')


            %Plot S+ and S- predictions
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end

            hFig = figure(figNo);

            set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

            hold on

            this_per_trial_sp_timecourse=handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse;

            CIsp = bootci(1000, @mean, this_per_trial_sp_timecourse);
            meansp=mean(this_per_trial_sp_timecourse,1);
            CIsp(1,:)=meansp-CIsp(1,:);
            CIsp(2,:)=CIsp(2,:)-meansp;

            [hlsp, hpsp] = boundedline(time_span',mean(this_per_trial_sp_timecourse,1)', CIsp','cmap',[0 114/255 178/255]);


            this_per_trial_sm_timecourse=handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse+1.2;

            CIsm = bootci(1000, @mean, this_per_trial_sm_timecourse);
            meansm=mean(this_per_trial_sm_timecourse,1);
            CIsm(1,:)=meansm-CIsm(1,:);
            CIsm(2,:)=CIsm(2,:)-meansm;

            [hlsm, hpsm] = boundedline(time_span',mean(this_per_trial_sm_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);

            plot(time_span',mean(this_per_trial_sp_timecourse,1)','-','Color',[0 114/255 178/255])
            plot(time_span',mean(this_per_trial_sm_timecourse,1)', '-','Color',[158/255 31/255 99/255]);

            text(30,0.75,'S+','Color',[0 114/255 178/255])
            text(30,0.65,'S-','Color',[158/255 31/255 99/255])


            ylim([-0.1 2.2])
            xlim([-7 15])
            this_ylim=ylim;

            %FV
            rectangle(Position=[-1.5,this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)),1.5,0.03*(this_ylim(2)-this_ylim(1))], FaceColor=[0.9 0.9 0.9], EdgeColor=[0.9 0.9 0.9])
            plot([-1.5 -1.5],[this_ylim],'-','Color',[0.5 0.5 0.5])

            %Odor on markers
            plot([0 0],this_ylim,'-k')
            odorhl=plot([0 mean(delta_odor)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],this_ylim,'-k')

            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],this_ylim,'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on+delta_reinf)],[this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1)) this_ylim(1)+0.1*(this_ylim(2)-this_ylim(1))],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on+delta_reinf) mean(delta_odor_on_reinf_on+delta_reinf)],this_ylim,'-r')

            title(['Prediction for ii_ROI '  num2str(iiROI)])
            xlabel('Time(sec)')
            ylabel('Prediction')
            pffft=1;
        end
    end
    %
    % %Plot the correct predictions
    % if show_figures==1
    %     for iiMLalgo=MLalgo_to_use
    %
    %         figNo=figNo+1;
    %         try
    %             close(figNo)
    %         catch
    %         end
    %
    %         hFig = figure(figNo);
    %
    %         set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    %
    %         hold on
    %         accuracy_tr=handles_out2.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
    %         accuracy_tr_other_ROIs=[];
    %         for iiROI=2:no_ROI_draws
    %             accuracy_tr_other_ROIs=[accuracy_tr_other_ROIs handles_out2.ROI(iiROI).MLalgo(iiMLalgo).accuracy_tr];
    %         end
    %
    %         histogram(accuracy_tr_other_ROIs)
    %         yl=ylim;
    %         ylim([yl(1) yl(2)+0.1*(yl(2)-yl(1))])
    %         plot([accuracy_tr accuracy_tr],yl,'-k','LineWidth',3)
    %         title(['Histogram of accuracy for p <= ' num2str(p_threshold) ' vs. other ROI sets'])
    %         xlabel('Accuracy')
    %     end
    % end


    %Save other variables to handles_out2
    handles_out2.time_span=time_span;
    handles_out2.MLalgo=MLalgo;
    handles_out2.time_windows=time_windows;
    handles_out2.wondow_no=window_no;
    handles_out2.no_ROI_draws=no_ROI_draws;
    handles_out2.dFF_per_trial_sm=dFF_per_trial_sm;
    handles_out2.dFF_per_trial_sp=dFF_per_trial_sp;

    %Plot the correct predictions for the odor window
    if show_figures==1
        %     for iiMLalgo=MLalgo_to_use
        %
        %         figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %
        %         hFig = figure(figNo);
        %
        %         set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
        %
        %         hold on
        % %         accuracy_tr=handles_out2.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
        %         accuracy_tr_ROIs=[];
        %         for iiROI=1:no_ROI_draws
        %             accuracy_tr_ROIs=[accuracy_tr_ROIs handles_out2.ROI(iiROI).MLalgo(iiMLalgo).accuracy_tr];
        %         end
        %
        %         histogram(accuracy_tr_ROIs)
        %         yl=ylim;
        %         ylim([yl(1) yl(2)+0.1*(yl(2)-yl(1))])
        % %         plot([accuracy_tr accuracy_tr],yl,'-k','LineWidth',3)
        %         title(['Histogram of accuracy for p <= ' num2str(p_threshold)])
        %         xlabel('Accuracy')
        %     end

        %Plot accuracy
        %     for iiMLalgo=MLalgo_to_use
        %This is odor window
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

        hold on
        edges=[0:0.05:1];



        %         accuracy_tr=handles_out2.ROI(1).MLalgo(iiMLalgo).accuracy_tr;
        accuracy_per_ROI=[];
        for iiROI=1:no_ROI_draws
            accuracy_per_ROI=[accuracy_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),2))];
        end

        histogram(accuracy_per_ROI,edges,'Normalization','Probability')



        accuracy_per_ROI_sh=[];
        for iiROI=1:no_ROI_draws
            accuracy_per_ROI_sh=[accuracy_per_ROI_sh mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),2))];
        end

        histogram(accuracy_per_ROI_sh,edges,'Normalization','Probability')


        title(['Histogram of accuracy for p <= ' num2str(p_threshold)])
        xlabel('Accuracy')

        fprintf(1,'Mean accuracy all runs= %d\n',mean(accuracy_per_ROI))
        fprintf(1,'Shuffled trial mean accuracy all runs= %d\n',mean(accuracy_per_ROI_sh))
        %     end

        %Plot pseudocolor accuracy vs dFF fraction positive
        %This is odor window
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

        hold on
        colors=colormap(jet);

        ii_t_start=find(time_span>=time_windows(window_no,1),1,'first');
        ii_t_end=find(time_span<=time_windows(window_no,2),1,'last');

        for iiROI=1:no_ROI_draws

            for ii_t=ii_t_start:ii_t_end
                dFF_sm_nonzero=0;
                n_dFF_sm=0;
                for trNo=1:size(dFF_per_trial_sm,1)
                    n_dFF_sm=n_dFF_sm+1;
                    if dFF_per_trial_sm(trNo,iiROI,ii_t)>0.05
                        dFF_sm_nonzero=dFF_sm_nonzero+1;
                    end
                end
                dFF_sm_nonzero_f=dFF_sm_nonzero/n_dFF_sm;

                dFF_sp_nonzero=0;
                n_dFF_sp=0;
                for trNo=1:size(dFF_per_trial_sp,1)
                    n_dFF_sp=n_dFF_sp+1;
                    if dFF_per_trial_sp(trNo,iiROI,ii_t)>0.05
                        dFF_sp_nonzero=dFF_sp_nonzero+1;
                    end
                end
                dFF_sp_nonzero_f=dFF_sp_nonzero/n_dFF_sp;
                if ((accuracy_per_ROI(iiROI)<0.3)||(accuracy_per_ROI(iiROI)>0.7))
                    ii_color=1+ceil(accuracy_per_ROI(iiROI)*255);
                    plot(dFF_sm_nonzero_f,dFF_sp_nonzero_f,'o','MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                        'MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                        'MarkerSize',5)
                end
            end
        end



        %Plot prediction

        %     for iiMLalgo=MLalgo_to_use
        %This is odor window

        %Do Sp first
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

        hold on

        edges=[0:0.1:1];

        %S+
        prediction_sp_per_ROI=[];
        for iiROI=1:no_ROI_draws
            prediction_sp_per_ROI=[prediction_sp_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),2))];
        end

        histogram(prediction_sp_per_ROI, edges)

        %S-
        prediction_sm_per_ROI=[];
        for iiROI=1:no_ROI_draws
            prediction_sm_per_ROI=[prediction_sm_per_ROI mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),2))];
        end

        histogram(prediction_sm_per_ROI, edges)


        title(['Histogram of label predictions for p <= ' num2str(p_threshold)])
        xlabel('Prediction')
        legend('S+','S-')





        %     end

        %Plot prediction S+ vs S-
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

        hold on


        for iiROI=1:no_ROI_draws
            ii_color=1+ceil(accuracy_per_ROI(iiROI)*255);
            plot(prediction_sm_per_ROI(iiROI), prediction_sp_per_ROI(iiROI),'o',...
                'MarkerFaceColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                'MarkerEdgeColor',[colors(ii_color,1) colors(ii_color,2) colors(ii_color,3)],...
                'MarkerSize',5)
        end


        title(['Prediction S+ vs. S-, color accuracy'])
        ylabel('Prediction S+')
        xlabel('Prediction S-')



        %Plot prediction d prime

        %     for iiMLalgo=MLalgo_to_use
        %This is odor window

        %Do Sp first
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

        hold on

        edges=[0:0.1:1];

        %Calculate d prime
        d_prime=[];
        for iiROI=1:no_ROI_draws
            these_Sp_predictions=mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),2);
            these_Sm_predictions=mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2))),2);
            this_d_prime=(mean(these_Sp_predictions)-mean(these_Sm_predictions))/sqrt((std(these_Sp_predictions)^2 + std(these_Sm_predictions)^2 )/2);

            d_prime=[d_prime this_d_prime];
        end

        histogram(d_prime)


        title(['d prime for label predictions for p <= ' num2str(p_threshold)])
        xlabel('d''')

        %     end




    end
    if show_figures==1
        fprintf(1,'Elapsed time (hr) %d\n',toc/(60*60))
    end
end
pffft=1;