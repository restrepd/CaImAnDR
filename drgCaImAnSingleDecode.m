function [outputArg1,outputArg2] = drgCaImAnSingleDecode(inputArg1,inputArg2)
%This function performs leave one out decoding case 3
%wand ouptuts the data for a few points

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
