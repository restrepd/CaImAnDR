function handles_out=drgCaImAn_pval_and_cc_batch(handles_choices)
%This program performs a demixed principal component analysis (dPCA)
%according to Kobak et al eLife 2016 DOI: 10.7554/eLife.10989
%
% the input is either no input in which case you will be asked to enter a pre_per file
% or a handles_choices variable with the pre_per file location and all
% other choices
%
% To find where you setup the dPCA search for "dPCA_input=="


handles_out=[];

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
process_pcorr=1; %If this is 1 then each group, percent correct is processed separately
do_dPCA=0; %dPCA is not performed, all the function does is to process the pre_per data to use in the p val code below
dPCA_input=1; %1 processes the pre_per file, 2 simulates Romo, 3 simulates go-no go


fr_per_names{1}='Fwd <40%%';
fr_per_names{2}='Fwd 40-65%%';
fr_per_names{3}='Fwd 65-80%%';
fr_per_names{4}='Fwd >=80%%';
fr_per_names{5}='Rev <40%%';
fr_per_names{6}='Rev 40-65%%';
fr_per_names{7}='Rev 65-80%%';
fr_per_names{8}='Rev >=80%%';

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
    allPcorr=[];
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

        all_handles(fileNo).handles_out=drgCaImAn_dPCA_stim_and_dec(handles_choices);

        pCorr=all_handles(fileNo).handles_out.percent_correct;
        allPcorr=[allPcorr pCorr];

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
        pCorr_per_file(fileNo)=pCorr;
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

    figureNo=0;
    cd(handles.choiceBatchPathName)
    mkdir figures

    for mouseNo=1:ii_mouse
        for grNo=1:max(all_grNos)
            if (sum((all_grNos==grNo)&(all_mice==mouseNo))>0)
                %Let's process this group/mouse
                maxE=max(all_Es((all_grNos==grNo)&(all_mice==mouseNo)));
                ii_included=0;
                allN=sum(all_Ns((all_grNos==grNo)&(all_mice==mouseNo)));

                %Note: These data were acquired at different rates. Here they are all
                %resampled to a dt of 0.03 sec
                dt_res=0.03;
                time_span=t_from:dt_res:t_to;
                T=length(time_span);



                handles_out.output_data_pre=[];
                handles_out.output_data_FV=[];
                handles_out.output_data_odor=[];
                handles_out.output_data_N=[];

                for fileNo=first_file:length(handles.FileName_pre_per)
                    if (all_grNos(fileNo)==grNo)&(all_mice(fileNo)==mouseNo)
                        handles_out.file(fileNo).mouseNo=mouseNo;
                        handles_out.file(fileNo).grNo=grNo;
                        handles_out.file(fileNo).pCorr=pCorr_per_file(fileNo);

                        N=all_Ns(fileNo);
                        firingRates=NaN(N,S,T,maxE);
                        trialNum=zeros(N,S);

                        ii_included=ii_included+1;

                        these_firingRates_d=all_handles(fileNo).handles_out.firingRates;
                        %Collapse decisions into stimuli
                        these_firingRates=NaN(size(these_firingRates_d,1),S,size(these_firingRates_d,4),size(these_firingRates_d,5));



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


                                    firingRates(n,s,1:T,trNo)=these_FR;
                                    trialNum(n,s)=all_trialNum(n,s);
                                end
                            end

                        end




                        %Now do p-vals
                        pre_t=[-3 -2];
                        fv_t=[-1 0];
                        odor_t=[2 3];

                        %do t-test or ranksum
                        handles_out.file(fileNo).no_neurons=N;


                        %Pre period
                        id_ii=0;
                        input_data=[];
                        for n=1:N

                            %S+
                            s=1;
                            id_ii=id_ii+1;
                            for trNo=1:trialNum(n,s)
                                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=pre_t(1))&(time_span<pre_t(2)),trNo),3);
                            end
                            input_data(id_ii).description=['Neuron ' num2str(n) ' S+'];

                            %S-
                            s=2;
                            id_ii=id_ii+1;
                            for trNo=1:trialNum(n,s)
                                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=pre_t(1))&(time_span<pre_t(2)),trNo),3);
                            end

                            input_data(id_ii).description=['Neuron ' num2str(n) ' S-'];
                        end

                        fprintf(1, ['Ranksum or t test for pre period for file No ' num2str(fileNo) '\n']);
                        fprintf(1, ['Percent correct behavior for this file is ' num2str(pCorr_per_file(fileNo)) '\n']);

                        [handles_out.file(fileNo).output_data_pre] = drgPairsRanksumorTtest(input_data);

                        %FV period
                        id_ii=0;
                        input_data=[];
                        for n=1:N

                            %S+
                            s=1;
                            id_ii=id_ii+1;
                            for trNo=1:trialNum(n,s)
                                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=fv_t(1))&(time_span<fv_t(2)),trNo),3);
                            end
                            input_data(id_ii).description=['Neuron ' num2str(n) ' S+'];

                            %S-
                            s=2;
                            id_ii=id_ii+1;
                            for trNo=1:trialNum(n,s)
                                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=fv_t(1))&(time_span<fv_t(2)),trNo),3);
                            end

                            input_data(id_ii).description=['Neuron ' num2str(n) ' S-'];
                        end

                        fprintf(1, ['Ranksum or t test for final valve period for file No ' num2str(fileNo) '\n']);
                        fprintf(1, ['Percent correct behavior for this file is ' num2str(pCorr_per_file(fileNo)) '\n']);

                        [handles_out.file(fileNo).output_data_FV] = drgPairsRanksumorTtest(input_data);

                        %Odor
                        id_ii=0;
                        input_data=[];
                        for n=1:N

                            %S+
                            s=1;
                            id_ii=id_ii+1;
                            for trNo=1:trialNum(n,s)
                                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=odor_t(1))&(time_span<odor_t(2)),trNo),3);
                            end
                            input_data(id_ii).description=['Neuron ' num2str(n) ' S+'];

                            %S-
                            s=2;
                            id_ii=id_ii+1;
                            for trNo=1:trialNum(n,s)
                                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=odor_t(1))&(time_span<odor_t(2)),trNo),3);
                            end

                            input_data(id_ii).description=['Neuron ' num2str(n) ' S-'];
                        end

                        fprintf(1, ['Ranksum or t test for odor period for file No ' num2str(fileNo) '\n']);
                        [handles_out.file(fileNo).output_data_odor] = drgPairsRanksumorTtest(input_data);

                        handles_out.output_data_pre=[handles_out.output_data_pre handles_out.file(fileNo).output_data_pre];
                        handles_out.output_data_FV=[handles_out.output_data_FV handles_out.file(fileNo).output_data_FV];
                        handles_out.output_data_odor=[handles_out.output_data_odor handles_out.file(fileNo).output_data_odor];
                        handles_out.output_data_N=[handles_out.output_data_N N];


                        %Plot the timecourse for the lowest p value for the
                        %odor period
                        figureNo = figureNo + 1;
                        try
                            close(figureNo)
                        catch
                        end
                        hFig=figure(figureNo);
                        hold on

                        %get the dF/F
                        [min_p min_ii]=min(handles_out.file(fileNo).output_data_odor.p);
                        dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                        dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                        dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                        dFFsminus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,2));

                        try
                            %S-
                            CIpv = bootci(1000, @mean, dFFsminus');
                            meanpv=mean(dFFsminus',1);
                            CIpv(1,:)=meanpv-CIpv(1,:);
                            CIpv(2,:)=CIpv(2,:)-meanpv;


                            [hlpvl, hppvl] = boundedline(time_span,mean(dFFsminus'), CIpv','cmap',[158/255 31/255 99/255]);

                            %S+
                            CIpv = bootci(1000, @mean, dFFsplus');
                            meanpv=mean(dFFsplus',1);
                            CIpv(1,:)=meanpv-CIpv(1,:);
                            CIpv(2,:)=CIpv(2,:)-meanpv;


                            [hlpvl, hppvl] = boundedline(time_span, mean(dFFsplus'), CIpv','cmap',[0 114/255 178/255]);


                            plot(time_span',mean(dFFsminus')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
                            plot(time_span',mean(dFFsplus')','Color',[0 114/255 178/255],'LineWidth',1.5);
                        catch
                        end

                        this_ylim=ylim;
                        plot([0 0],this_ylim,'-k')
                        xlim([-7 15])

                        xlabel('Time(sec)')
                        ylabel('dFF')

                        title(['dFF for neuron ' num2str(min_ii), ' for file No ' num2str(fileNo)])

                        savefig([handles.choiceBatchPathName 'figures\' 'fil' num2str(fileNo) 'neuron' num2str(min_ii)])
                    end
                end


                pffft=1;
            end

        end
    end
end


pffft=1;