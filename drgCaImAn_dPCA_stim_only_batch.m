function handles_out=drgCaImAn_dPCA_stim_only_batch(handles_choices)
%This program performs a demixed principal component analysis (dPCA) 
%according to Kobak et al eLife 2016 DOI: 10.7554/eLife.10989
%
% the input is either no input in which case you will be asked to enter a pre_per file 
% or a handles_choices variable with the pre_per file location and all
% other choices
%
% To find where you setup the dPCA search for "dPCA_input=="




if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dPCA_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_LDA_fsdz run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.no_files;

if isfield(handles,'processing_algo')
    processing_algo=handles.processing_algo;
else
    processing_algo=1;
end

if isfield(handles,'suffix_out')
    suffix_out=handles.suffix_out;
else
    suffix_out='_dec.mat';
end

if isfield(handles,'first_file')
    first_file=handles.first_file;
else
    first_file=1;
end

%Choices for simulations and how files are grouped
process_pcorr=0; %If this is 1 then each group, percent correct is processed separately
do_dPCA=0; %dPCA is not performed, all the function does is to process the pre_per data to use in the dPCA code below
dPCA_input=1; %1 processes the pre_per file, 2 simulates Romo, 3 simulates go-no go


fr_per_names{1}='Fwd <40%';
fr_per_names{2}='Fwd 40-65%';
fr_per_names{3}='Fwd 65-80%';
fr_per_names{4}='Fwd >=80%';
fr_per_names{5}='Rev <40%';
fr_per_names{6}='Rev 40-65%';
fr_per_names{7}='Rev 65-80%';
fr_per_names{8}='Rev >=80%';

per_names{1}='40-65%';
per_names{2}='65-80%';
per_names{3}='>=80%';

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

%Time tinterval for shifting time base due to slow olfactometer computer
t_shift=0.61;


%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:handles.no_files


    %Make sure that all the files exist
    pre_per_FileName=handles.FileName_pre_per{filNum};
    if iscell(handles.PathName_pre_per)
        pre_per_PathName=handles.PathName_pre_per{filNum};
    else
        pre_per_PathName=handles.PathName_pre_per;
    end

    if exist([pre_per_PathName pre_per_FileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' pre_per_FileName ' does not exist\n'],filNum);
        all_files_present=0;
    end

end


% if exist([handles.PathName_out handles.FileName_out])==0
%     handles_out=[];
%     ii_out=0;
% else
%     load([handles.PathName_out handles.FileName_out])
%     ii_out=handles_out.last_ii_out;
%     first_file=handles_out.last_file_processed+1;
% end

figNo=0;
show_figures=1;

if all_files_present==1


    %Find the group per file
    no_pcorr=4;
    all_handles=[];
    all_Es=[];
    all_grNos=[];
    allPcorr_per_file=[];
    all_mice=[];
    all_Ns=[];
    all_t=[];
    these_groups=unique(handles.group);
    for fileNo=first_file:length(handles.FileName_pre_per)
        tic



        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};

        handles_choices.pre_per_PathName=pre_per_PathName;
        handles_choices.pre_per_FileName=pre_per_FileName;

        handles_choices.show_figures=handles.show_figures;

        handles_choices.do_dPCA=do_dPCA; %dPCA is not performed, all the function does is to process the pre_per data to use in the dPCA code below
        handles_choices.dPCA_input=dPCA_input; %1 processes the pre_per file, 2 simulates Romo, 3 simulates go-no go
        handles_choices.convert_z=2; %This ensures that the components are the same before odor on

%         all_handles(fileNo).handles_out=drgCaImAn_dPCA_stim_and_dec(handles_choices);
        all_handles(fileNo).handles_out=drgCaImAn_get_dFF_per_trial(handles_choices);

        pCorr=all_handles(fileNo).handles_out.percent_correct;
        allPcorr_per_file=[allPcorr_per_file pCorr];
 
        if fileNo==1
            ii_mouse=1;
            mouse_names{1}=handles.mouse{1};
            all_mice(1)=1;
        else
            %Find whether this mouse is already in the list
            mouse_found=0;
            for this_ii=1:length(mouse_names)
                if strcmp(handles.mouse{fileNo},mouse_names{this_ii})
                    mouse_found=1;
                    mouse_found_ii=this_ii;
                end
            end
            if mouse_found==0
                %This mouse's name is not in the list
                ii_mouse=ii_mouse+1;
                mouse_names{ii_mouse}=handles.mouse{fileNo};
                all_mice(fileNo)=ii_mouse;
            else
                %This mouse's name is in the list
                all_mice(fileNo)=mouse_found_ii;
            end
        end

        ii_pCorr=1;
        if (pCorr>=40)&(pCorr<=65)
            ii_pCorr=2;
        else
            if (pCorr>65)&(pCorr<80)
                ii_pCorr=3;
            else
                if pCorr>=80
                    ii_pCorr=4;
                end
            end
        end
        this_group=handles.group(fileNo);
        this_grNo=find(these_groups==this_group);
        all_grNos=[all_grNos (this_grNo-1)*no_pcorr+ii_pCorr];

        all_Es=[all_Es all_handles(fileNo).handles_out.E];

        all_Ns(fileNo)=all_handles(fileNo).handles_out.N;
        S=all_handles(fileNo).handles_out.S;
        D=all_handles(fileNo).handles_out.D;
        all_t(fileNo).time=all_handles(fileNo).handles_out.time;
        all_t(fileNo).T=all_handles(fileNo).handles_out.T;
        all_t(fileNo).trialNum=all_handles(fileNo).handles_out.trialNum;
        t_from=all_handles(fileNo).handles_out.trial_time_from;
        t_to=all_handles(fileNo).handles_out.trial_time_to;

    end

if process_pcorr==0
    all_grNos=ones(1,length(all_grNos));
end

    for mouseNo=1:ii_mouse
        for grNo=1:max(all_grNos)
            if (sum((all_grNos==grNo)&(all_mice==mouseNo))>0)
                %Let's process this group/mouse
                maxE=max(all_Es((all_grNos==grNo)&(all_mice==mouseNo)));
                ii_included=0;
                allN=sum(all_Ns((all_grNos==grNo)&(all_mice==mouseNo)));
                pcorr_per_n=zeros(1,allN);

                %Note: These data were acquired at different rates. Here they are all
                %resampled to a dt of 0.03 sec
                dt_res=0.03;
                time_span=t_from:dt_res:t_to;
                T=length(time_span);
                firingRates=NaN(allN,S,T,maxE);
                trialNum=zeros(allN,S);
                this_n=0;


                for fileNo=first_file:length(handles.FileName_pre_per)
                    if (all_grNos(fileNo)==grNo)&(all_mice(fileNo)==mouseNo)
                        ii_included=ii_included+1;
                        
                        these_firingRates_d=all_handles(fileNo).handles_out.firingRates;
                        %Collapse decisions into stimuli
                        these_firingRates=NaN(size(these_firingRates_d,1),S,size(these_firingRates_d,4),size(these_firingRates_d,5));
                        
                        N=all_Ns(fileNo);
                        
                        for n=1:N
                            for s=1:S
                                for d=1:D
                                    for trNo=1:all_t(fileNo).trialNum(n,s,d)
                                        if d==2
                                            trNo_zero=all_t(fileNo).trialNum(n,s,1);
                                        else
                                            trNo_zero=0;
                                        end
                                        for ii_tsp=1:size(these_firingRates_d,4)
                                            if ~isnan(these_firingRates_d(n,s,d,ii_tsp,trNo))
                                                these_firingRates(n,s,ii_tsp,trNo+trNo_zero)=these_firingRates_d(n,s,d,ii_tsp,trNo);
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        
                        all_trialNum=zeros(N,S);
                        for n=1:N
                            for s=1:S
                                all_trialNum(n,s)=all_t(fileNo).trialNum(n,s,1)+all_t(fileNo).trialNum(n,s,2);
                            end
                        end

                        this_ii_tspan=all_t(fileNo).T;
                        this_time_span=all_t(fileNo).time;
                        N=all_Ns(fileNo);
                        
                        for n=1:N
                            for s=1:S
                                for trNo=1:all_trialNum(n,s)
                                    these_FR=zeros(1,length(time_span));

                                    for ii_tsp=1:length(time_span)
                                        if time_span(ii_tsp)<this_time_span(1)
                                            these_FR(ii_tsp)=these_firingRates(n,s,1,trNo);
                                        else
                                            if time_span(ii_tsp)>this_time_span(end)
                                                these_FR(ii_tsp)=these_firingRates(n,s,this_ii_tspan,trNo);
                                            else
                                                ii_0=find(this_time_span<=time_span(ii_tsp),1,'last');
                                                ii_1=find(this_time_span>time_span(ii_tsp),1,'first');
                                                these_FR(ii_tsp)=these_firingRates(n,s,ii_0,trNo)+...
                                                    (these_firingRates(n,s,ii_1,trNo)-these_firingRates(n,s,ii_0,trNo))*...
                                                    (time_span(ii_tsp)-this_time_span(ii_0))/(this_time_span(ii_1)-this_time_span(ii_0));
                                            end
                                        end
                                    end

                                    %Now do the time shift
                                    no_shift=ceil(t_shift/(time_span(2)-time_span(1)));

                                    if all_handles(fileNo).handles_out.shift_time==1
                                        if no_shift>0
                                            these_FR(no_shift+1:end)=these_FR(1:end-no_shift);
                                            these_FR(1:no_shift)=these_FR(no_shift);
                                        else
                                            these_FR(1:end-no_shift)=these_FR(no_shift+1:end);
                                            these_FR(1:end-no_shift+1:end)=these_FR(no_shift);
                                        end
                                    end


                                    firingRates(this_n+n,s,1:T,trNo)=these_FR;
                                    trialNum(this_n+n,s)=all_trialNum(n,s);
                                end
                            end
                            pcorr_per_n(this_n+n)=allPcorr_per_file(fileNo);
                        end
                        this_n=this_n+n;
                    end
                end


               

                %Now do dPCA
               
                % firingRatesAverage: N x S x T
                firingRatesAverage = nanmean(firingRates, 4);
                firingRatesAverage(isnan(firingRatesAverage))=0; %This is key when there are zero trials for some s,d combinations.

                numComp=10; %Components to be extracted from data by dpca
                lambda=0.00001; %Note: if I do not make this larger than thisLambda = 1e-10; in line 131 of dpca
                time=time_span;

                    %% Define parameter grouping

                    % *** Don't change this if you don't know what you are doing! ***
                    % firingRates array has [N S D T E] size; here we ignore the 1st dimension
                    % (neurons), i.e. we have the following parameters:
                    %    1 - stimulus
                    %    2 - decision
                    %    3 - time
                    % There are three pairwise interactions:
                    %    [1 3] - stimulus/time interaction
                    %    [2 3] - decision/time interaction
                    %    [1 2] - stimulus/decision interaction
                    % And one three-way interaction:
                    %    [1 2 3] - rest
                    % As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

                    % combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
                    % margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
                    % margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

                    % For two parameters (e.g. stimulus and time, but no decision), we would have
                    % firingRates array of [N S T E] size (one dimension less, and only the following
                    % possible marginalizations:
                    %    1 - stimulus
                    %    2 - time
                    %    [1 2] - stimulus/time interaction
                    % They could be grouped as follows:
                    %    combinedParams = {{1, [1 2]}, {2}};


                    combinedParams = {{1, [1 2]}, {2}};
                    margNames = {'Stimulus', 'Stimulus-independent'};
                    margColours = [23 100 171; 187 20 25]/256;

                    % Time events of interest (e.g. stimulus onset/offset, cues etc.)
                    % They are marked on the plots with vertical lines
                    timeEvents = [0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

                    % check consistency between trialNum and firingRates
                    for n = 1:size(firingRates,1)
                        for s = 1:size(firingRates,2)
                            assert(isempty(find(isnan(firingRates(n,s,:,1:trialNum(n,s))), 1)), 'Something is wrong!')
                        end
                    end

                    % dPCA without regularization and ignoring noise covariance

                    % This is the core function.
                    % W is the decoder, V is the encoder (ordered by explained variance),
                    % whichMarg is an array that tells you which component comes from which
                    % marginalization

                    tic
                    [W,V,whichMarg] = dpca(firingRatesAverage, numComp, ...
                        'combinedParams', combinedParams,'lambda',lambda);
                    toc

                    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
                        'combinedParams', combinedParams);

                    dpca_handles=dpca_plot_drg(firingRatesAverage, W, V, @dpca_plot_default_drg, ...
                        'explainedVar', explVar, ...
                        'marginalizationNames', margNames, ...
                        'marginalizationColours', margColours, ...
                        'whichMarg', whichMarg,                 ...
                        'time', time,                        ...
                        'timeEvents', timeEvents,               ...
                        'timeMarginalization', 3, ...
                        'legendSubplot', 16,...
                        'ylims',[20 20 20 20], ...
                        'xlims',[-5 15], ...
                        'displaySumStimComps',1);

                    if process_pcorr==1
                        t=sgtitle(['Demixed PCA for mouse No ' num2str(mouseNo) ' ' fr_per_names{grNo}]);
                    else
                        t=sgtitle(['Demixed PCA for mouse No ' num2str(mouseNo) ' entire dataset']);
                    end
                    t.HorizontalAlignment='left';

                    %Now calculate the components using the decoder matrix W
                    if process_pcorr==0
                        for grNo=1:3
                            
                            switch grNo
                                case 1
                                    ns_included=(pcorr_per_n>=45)&(pcorr_per_n<=65);
                                case 2
                                    ns_included=(pcorr_per_n>65)&(pcorr_per_n<80);
                                case 3
                                    ns_included=(pcorr_per_n>=80);
                            end

                            alln_ii=1:length(ns_included);
                            ii_ns_included=alln_ii(ns_included);
                            ns_excluded=~ns_included;
                            ii_ns_excluded=alln_ii(ns_excluded);

                            thesefiringRatesAverage=firingRatesAverage;
                            ii_inc=0;
                            for n=1:allN
                                if sum(n==ii_ns_excluded)>0
                                    ii_inc=ii_inc+1;
                                    if ii_inc>length(ii_ns_included)
                                        ii_inc=1;
                                    end
                                    for s=1:S
                                        for ii_tsp=1:length(time_span)
                                            thesefiringRatesAverage(n,s,ii_tsp)=firingRatesAverage(ii_inc,s,ii_tsp);
                                        end
                                    end
                                end
                            end


                            dpca_handles=dpca_plot_drg(thesefiringRatesAverage, W, V, @dpca_plot_default_drg, ...
                                'explainedVar', explVar, ...
                                'marginalizationNames', margNames, ...
                                'marginalizationColours', margColours, ...
                                'whichMarg', whichMarg,                 ...
                                'time', time,                        ...
                                'timeEvents', timeEvents,               ...
                                'timeMarginalization', 3, ...
                                'legendSubplot', 16,...
                                'ylims',[20 20 20 20], ...
                                'xlims',[-5 15], ...
                                'displaySumStimComps',1);


                            t=sgtitle(['Demixed PCA for mouse No ' num2str(mouseNo) ' ' per_names{grNo}]);

                            t.HorizontalAlignment='left';
                        end
                    end

                    %% Step 4: dPCA with regularization did not work because many of theneurons have zero trials for some conditions
    
% 
%                     % This function takes some minutes to run. It will save the computations
%                     % in a .mat file with a given name. Once computed, you can simply load
%                     % lambdas out of this file:
%                     %   load('tmp_optimalLambdas.mat', 'optimalLambda')
%   
%                     % Please note that this now includes noise covariance matrix Cnoise which
%                     % tends to provide substantial regularization by itself (even with lambda set
%                     % to zero).
%                     ifSimultaneousRecording = false;
%                     optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
%                         'combinedParams', combinedParams, ...
%                         'simultaneous', ifSimultaneousRecording, ...
%                         'numRep', 2, ...  % increase this number to ~10 for better accuracy
%                         'filename', 'tmp_optimalLambdas.mat');
% 
%                     Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
%                         firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);
% 
%                     [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
%                         'combinedParams', combinedParams, ...
%                         'lambda', optimalLambda, ...
%                         'Cnoise', Cnoise);
% 
%                     explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%                         'combinedParams', combinedParams);
% 
%                     dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%                         'explainedVar', explVar, ...
%                         'marginalizationNames', margNames, ...
%                         'marginalizationColours', margColours, ...
%                         'whichMarg', whichMarg,                 ...
%                         'time', time,                        ...
%                         'timeEvents', timeEvents,               ...
%                         'timeMarginalization', 3,           ...
%                         'legendSubplot', 16);
                    pffft=1;
            end

        end
    end
end
 

pffft=1;