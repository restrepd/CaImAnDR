function handles_out=drgCaImAn_pval_batch(handles_choices)
%This program performs a demixed principal component analysis (dPCA)
%according to Kobak et al eLife 2016 DOI: 10.7554/eLife.10989
%
% the input is either no input in which case you will be asked to enter a pre_per file
% or a handles_choices variable with the pre_per file location and all
% other choices
%
% To find where you setup the dPCA search for "dPCA_input=="

close all
clear all

handles_out=[];

show_figures=0;

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


fr_per_names{1}='Fwd <40%';
fr_per_names{2}='Fwd 40-65%';
fr_per_names{3}='Fwd 65-80%';
fr_per_names{4}='Fwd >=80%';
fr_per_names{5}='Rev <40%';
fr_per_names{6}='Rev 40-65%';
fr_per_names{7}='Rev 65-80%';
fr_per_names{8}='Rev >=80%';

per_names{2}='40-65%';
per_names{3}='65-80%';
per_names{4}='>=80%';

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

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

    first_file=handles.first_file;
    first_file_processed=0;

    all_dFFsplus=[];
    all_dFFsminus=[];
    all_dFFspm=[];
    ii_dFF=0;

    for fileNo=first_file:length(handles.FileName_pre_per)
        tic

        fprintf(1, ['Processing file No  ' num2str(fileNo) '\n']);

        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};

        handles_choices.pre_per_PathName=pre_per_PathName;
        handles_choices.pre_per_FileName=pre_per_FileName;

        handles_choices.show_figures=handles.show_figures;

        handles_choices.do_dPCA=do_dPCA; %dPCA is not performed, all the function does is to process the pre_per data to use in the dPCA code below
        handles_choices.dPCA_input=dPCA_input; %1 processes the pre_per file, 2 simulates Romo, 3 simulates go-no go
        handles_choices.convert_z=2; %This ensures that the components are the same before odor on
 
        all_handles(fileNo).handles_out=drgCaImAn_get_dFF_per_trial(handles_choices);


        pCorr=all_handles(fileNo).handles_out.percent_correct;
        allPcorr=[allPcorr pCorr];

        if first_file_processed==0
            ii_mouse=1;
            mouse_names{1}=handles.mouse{1};
            all_mice(1)=1;
            first_file_processed=1;
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

    %Let's process all files
    figureNo=0;
    cd(handles.choiceBatchPathName)

    try
        rmdir('figures','s')
    catch
    end
    mkdir figures

    ii_included=0;

    %Note: These data were acquired at different rates. Here they are all
    %resampled to a dt of 0.03 sec
    dt_res=0.03;
    time_span=t_from:dt_res:t_to;
    T=length(time_span);

    handles_out.output_data_pre=[];
    handles_out.output_data_FV=[];
    handles_out.output_data_odor=[];
    handles_out.output_data_N=[];
    handles_out.mouseNo=[];
    handles_out.perCorr=[];
    handles_out.grNo=[];
    handles_out.output_data_Sp=[];
    handles_out.output_data_Sm=[];
    handles_out.output_data_SporSm=[];

    for fileNo=first_file:length(handles.FileName_pre_per)

        mouseNo=all_mice(fileNo);
        handles_out.file(fileNo).mouseNo=mouseNo;
        handles_out.file(fileNo).pCorr=pCorr_per_file(fileNo);
        handles_out.file(fileNo).grNo=handles.group(fileNo);

        handles_out.grNo=[handles_out.grNo handles.group(fileNo)];
        handles_out.perCorr=[handles_out.perCorr pCorr_per_file(fileNo)];
        handles_out.mouseNo=[handles_out.mouseNo mouseNo];



        N=all_Ns(fileNo);
        firingRates=NaN(N,S,T,all_Es(fileNo));
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

        handles_out.output_data_pre=[handles_out.output_data_pre handles_out.file(fileNo).output_data_pre.no_significant_pFDR];
        handles_out.output_data_FV=[handles_out.output_data_FV handles_out.file(fileNo).output_data_FV.no_significant_pFDR];
        handles_out.output_data_odor=[handles_out.output_data_odor handles_out.file(fileNo).output_data_odor.no_significant_pFDR];
        handles_out.output_data_N=[handles_out.output_data_N N];



        %Show the figure for divergent dFF in the odor period if p<pFDR
        %Save the dFFSplus and dFFSminus
      
        for ii_p=1:length(handles_out.file(fileNo).output_data_odor.p)
            this_p=handles_out.file(fileNo).output_data_odor.p(ii_p);
            min_ii=ii_p;
            if this_p<handles_out.file(fileNo).output_data_odor.pFDR
                if show_figures==1
                    figureNo = figureNo + 1;
                    try
                        close(figureNo)
                    catch
                    end
                    hFig=figure(figureNo);

                    ax=gca;ax.LineWidth=3;
                    set(hFig, 'units','normalized','position',[.2 .2 .25 .25])


                    hold on
                end

                %get the dF/F

                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

                try
                    %S-
                    CIpv = bootci(1000, @mean, dFFsminus');
                    meanpv=mean(dFFsminus',1);
                    CIpv(1,:)=meanpv-CIpv(1,:);
                    CIpv(2,:)=CIpv(2,:)-meanpv;

                    if show_figures==1
                        [hlpvl, hppvl] = boundedline(time_span,mean(dFFsminus'), CIpv','cmap',[158/255 31/255 99/255]);
                    end

                    %S+
                    CIpv = bootci(1000, @mean, dFFsplus');
                    meanpv=mean(dFFsplus',1);
                    CIpv(1,:)=meanpv-CIpv(1,:);
                    CIpv(2,:)=CIpv(2,:)-meanpv;

                    if show_figures==1
                        [hlpvl, hppvl] = boundedline(time_span, mean(dFFsplus'), CIpv','cmap',[0 114/255 178/255]);


                        plot(time_span',mean(dFFsminus')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
                        plot(time_span',mean(dFFsplus')','Color',[0 114/255 178/255],'LineWidth',1.5);
                    end

                    ii_dFF=ii_dFF+1;
                    all_dFFspm_group(ii_dFF)=handles.group(fileNo);
                    all_dFFspm_pcorr(ii_dFF)=handles_out.perCorr(fileNo);
                    all_dFFsplus(ii_dFF,1:length(time_span))=mean(dFFsplus')';
                    all_dFFsminus(ii_dFF,1:length(time_span))=mean(dFFsminus')';
                    all_dFFspm(1:length(time_span),ii_dFF)=mean(dFFsplus')';
                    all_dFFspm(length(time_span)+1:2*length(time_span),ii_dFF)=mean(dFFsminus')';

                catch
                end

                if show_figures==1
                    this_ylim=ylim;
                    for ii_te=1:length(timeEvents)
                        plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
                    end

                    text(-5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
                    text(-5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)

                    xlim([-7 15])

                    xlabel('Time(sec)')
                    ylabel('dFF')

                    if all_handles(fileNo).handles_out.shift_time==1
                        title(['dFF odor n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo)) ' tsh'])
                    else
                        title(['dFF odor n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo))])
                    end

                    savefig([handles.choiceBatchPathName 'figures\' 'odor_f' num2str(fileNo) 'n' num2str(min_ii)])
                end
            end
        end

        %Plot the timecourse for the lowest p value for the
        %FV period
        [min_p min_ii]=min(handles_out.file(fileNo).output_data_FV.p);

        %Show the figure only if min_p<pFDR
        if min_p<handles_out.file(fileNo).output_data_FV.pFDR
            if show_figures==1
                figureNo = figureNo + 1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);

                ax=gca;ax.LineWidth=3;
                set(hFig, 'units','normalized','position',[.2 .2 .25 .25])


                hold on

                %get the dF/F

                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

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
                for ii_te=1:length(timeEvents)
                    plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
                end

                text(-5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
                text(-5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)

                xlim([-7 15])


                xlabel('Time(sec)')
                ylabel('dFF')

                if all_handles(fileNo).handles_out.shift_time==1
                    title(['dFF FV n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo)) ' tsh'])
                else
                    title(['dFF FV n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo))])
                end

                savefig([handles.choiceBatchPathName 'figures\' 'odor_f' num2str(fileNo) 'n' num2str(min_ii)])
            end
        end

        %Is there an odorant response (regardless of divergence)
        %Pre period

        %S+ responses
        id_ii=0;
        input_data=[];
        for n=1:N

            %S+ pre
            s=1;
            id_ii=id_ii+1;
            for trNo=1:trialNum(n,s)
                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=pre_t(1))&(time_span<pre_t(2)),trNo),3);
            end
            input_data(id_ii).description=['Neuron ' num2str(n) ' S+'];

            %S+ odor
            s=1;
            id_ii=id_ii+1;
            for trNo=1:trialNum(n,s)
                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=odor_t(1))&(time_span<odor_t(2)),trNo),3);
            end

            input_data(id_ii).description=['Neuron ' num2str(n) ' S-'];
        end

        fprintf(1, ['Ranksum or t test for S+ odor response for file No ' num2str(fileNo) '\n']);
        fprintf(1, ['Percent correct behavior for this file is ' num2str(pCorr_per_file(fileNo)) '\n']);

        [handles_out.file(fileNo).output_data_odor_Sp] = drgPairsRanksumorTtest(input_data);

        %Plot the timecourse if min_p<pFDR

        [min_p min_ii]=min(handles_out.file(fileNo).output_data_odor_Sp.p);

        %         for ii_p=1:length(handles_out.file(fileNo).output_data_odor_Sp.p)
        %             this_p=handles_out.file(fileNo).output_data_odor_Sp.p(ii_p);
        %             min_ii=ii_p;
        %             if this_p<handles_out.file(fileNo).output_data_odor_Sp.pFDR
        if min_p<handles_out.file(fileNo).output_data_odor_Sp.pFDR
            if show_figures==1
                figureNo = figureNo + 1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);

                ax=gca;ax.LineWidth=3;
                set(hFig, 'units','normalized','position',[.3 .2 .25 .25])


                hold on

                %get the dF/F

                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

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
                for ii_te=1:length(timeEvents)
                    plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
                end

                text(-5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
                text(-5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)

                xlim([-7 15])


                xlabel('Time(sec)')
                ylabel('dFF')

                if all_handles(fileNo).handles_out.shift_time==1
                    title(['dFF S+ n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo)) ' tsh'])
                else
                    title(['dFF S+ n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo))])
                end

                savefig([handles.choiceBatchPathName 'figures\' 'Sp_f' num2str(fileNo) 'n' num2str(min_ii)])
            end
        end
        %         end

        %S- responses
        id_ii=0;
        input_data=[];
        for n=1:N

            %S- pre
            s=2;
            id_ii=id_ii+1;
            for trNo=1:trialNum(n,s)
                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=pre_t(1))&(time_span<pre_t(2)),trNo),3);
            end
            input_data(id_ii).description=['Neuron ' num2str(n) ' S+'];

            %S- odor
            s=2;
            id_ii=id_ii+1;
            for trNo=1:trialNum(n,s)
                input_data(id_ii).data(trNo)=mean(firingRates(n,s,(time_span>=odor_t(1))&(time_span<odor_t(2)),trNo),3);
            end

            input_data(id_ii).description=['Neuron ' num2str(n) ' S-'];
        end

        fprintf(1, ['Ranksum or t test for S- odor response for file No ' num2str(fileNo) '\n']);
        fprintf(1, ['Percent correct behavior for this file is ' num2str(pCorr_per_file(fileNo)) '\n']);

        [handles_out.file(fileNo).output_data_odor_Sm] = drgPairsRanksumorTtest(input_data);

        handles_out.output_data_Sp=[handles_out.output_data_Sp handles_out.file(fileNo).output_data_odor_Sp.no_significant_pFDR];
        handles_out.output_data_Sm=[handles_out.output_data_Sm handles_out.file(fileNo).output_data_odor_Sm.no_significant_pFDR];
        odor_responses=logical(handles_out.file(fileNo).output_data_odor_Sp.p<handles_out.file(fileNo).output_data_odor_Sp.pFDR)|...
            logical(handles_out.file(fileNo).output_data_odor_Sm.p<handles_out.file(fileNo).output_data_odor_Sm.pFDR);
        no_odor_responses=sum(odor_responses);
        handles_out.output_data_SporSm=[handles_out.output_data_SporSm no_odor_responses];



        %Plot the timecourse for p<pFDR for the
        %S- response
        %         for ii_p=1:length(handles_out.file(fileNo).output_data_odor_Sm.p)
        %             this_p=handles_out.file(fileNo).output_data_odor_Sm.p(ii_p);
        %             min_ii=ii_p;
        %             %Show the figure only if min_p<pFDR
        %             if this_p<handles_out.file(fileNo).output_data_odor_Sm.pFDR


        [min_p min_ii]=min(handles_out.file(fileNo).output_data_odor_Sm.p);


        if min_p<handles_out.file(fileNo).output_data_odor_Sm.pFDR

            if show_figures==1
                figureNo = figureNo + 1;
                try
                    close(figureNo)
                catch
                end
                hFig=figure(figureNo);

                ax=gca;ax.LineWidth=3;
                set(hFig, 'units','normalized','position',[.2 .4 .25 .25])


                hold on

                %get the dF/F

                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

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
                for ii_te=1:length(timeEvents)
                    plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
                end

                text(-5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
                text(-5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)

                xlim([-7 15])


                xlabel('Time(sec)')
                ylabel('dFF')

                if all_handles(fileNo).handles_out.shift_time==1
                    title(['dFF S- n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo)) ' tsh'])
                else
                    title(['dFF S- n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo))])
                end

                savefig([handles.choiceBatchPathName 'figures\' 'Sp_f' num2str(fileNo) 'n' num2str(min_ii)])
            end
        end

        %         end
    end

    %Plot a bar graph showing percent signifcant divergence
    glm_sig=[];
    glm_sig_ii=0;
    input_sig_data=[];
    id_sig_ii=0;

    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

    hold on
    bar_offset=0;
    edges=[0:0.2:10];
    rand_offset=0.8;

    for grNo=1:3
        these_pre=[];
        these_FV=[];
        these_odor=[];

        for mouseNo=unique(handles_out.mouseNo)
            switch grNo
                case 1
                    files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                case 2
                    files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                case 3
                    files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
            end

            if sum(files_included)>0
                n_pre=zeros(1,sum(files_included));
                n_all=zeros(1,sum(files_included));
                n_pre(1,:)=handles_out.output_data_pre(files_included);
                n_all(1,:)=handles_out.output_data_N(files_included);
                this_pre=n_pre./n_all;
                these_pre=[these_pre 100*mean(this_pre)];

                n_FV=zeros(1,sum(files_included));
                n_FV(1,:)=handles_out.output_data_pre(files_included);
                this_FV=n_FV./n_all;
                these_FV=[these_FV 100*mean(this_FV)];

                n_odor=zeros(1,sum(files_included));
                n_odor(1,:)=handles_out.output_data_odor(files_included);
                this_odor=n_odor./n_all;
                these_odor=[these_odor 100*mean(this_odor)];
            end

        end

        %         bar_offset=bar_offset+1;
        %         bar(bar_offset,mean(these_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        %         if length(these_pre)>2
        %             %Violin plot
        %             [mean_out, CIout]=drgViolinPoint(these_pre...
        %                 ,edges,bar_offset,rand_offset,'k','k',3);
        %         end
        %
        %         glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_pre))=these_pre;
        %         glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_pre))=grNo*ones(1,length(these_pre));
        %         glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_pre))=1*ones(1,length(these_pre));
        %         glm_sig_ii=glm_sig_ii+length(these_pre);
        %
        %         id_sig_ii=id_sig_ii+1;
        %         input_sig_data(id_sig_ii).data=these_pre;
        %         input_sig_data(id_sig_ii).description=['Pre ' fr_per_names{grNo}];
        %
        %                 bar_offset=bar_offset+1;
        %         bar(bar_offset,mean(these_FV),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        %         if length(these_FV)>2
        %             %Violin plot
        %             [mean_out, CIout]=drgViolinPoint(these_FV...
        %                 ,edges,bar_offset,rand_offset,'k','k',3);
        %         end
        %
        %              glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_FV))=these_FV;
        %         glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_FV))=grNo*ones(1,length(these_FV));
        %         glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_FV))=1*ones(1,length(these_FV));
        %         glm_sig_ii=glm_sig_ii+length(these_FV);
        %
        %         id_sig_ii=id_sig_ii+1;
        %         input_sig_data(id_sig_ii).data=these_FV;
        %         input_sig_data(id_sig_ii).description=['FV ' fr_per_names{grNo}];


        bar_offset=bar_offset+1;
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',3);
        end

        glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
        glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
        glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
        glm_sig_ii=glm_sig_ii+length(these_odor);

        id_sig_ii=id_sig_ii+1;
        input_sig_data(id_sig_ii).data=these_odor;
        input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];

        bar_offset=bar_offset+1;
    end

    xticks([1 3 5])
    labels='xticklabels({';
    for ii_label=2:4
            labels=[labels '''' per_names{ii_label} ''', '];
    end
    labels=[labels(1:end-2) '})'];
    eval(labels)

    ylabel('% divergent')
    title('Percent divergent ROIs')

    %Perform the glm
    fprintf(1, ['\nglm for percent divergent\n'])


    tbl = table(glm_sig.data',glm_sig.grNo',...
        'VariableNames',{'percent_divergent','percent_correct'});
    mdl = fitglm(tbl,'percent_divergent~percent_correct'...
        ,'CategoricalVars',[2])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);


    %Plot a bar graph showing percent responding to S+
    glm_sig=[];
    glm_sig_ii=0;
    input_sig_data=[];
    id_sig_ii=0;

    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

    hold on
    bar_offset=0;
    edges=[0:0.2:10];
    rand_offset=0.8;

    for grNo=1:3
        these_pre=[];
        these_FV=[];
        these_odor=[];

        for mouseNo=unique(handles_out.mouseNo)
            switch grNo
                case 1
                    files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                case 2
                    files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                case 3
                    files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
            end

            if sum(files_included)>0
                n_all=zeros(1,sum(files_included));
                n_all(1,:)=handles_out.output_data_N(files_included);
                n_odor=zeros(1,sum(files_included));
                n_odor(1,:)=handles_out.output_data_Sp(files_included);
                this_odor=n_odor./n_all;
                these_odor=[these_odor 100*mean(this_odor)];
            end

        end

        %         bar_offset=bar_offset+1;
        %         bar(bar_offset,mean(these_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        %         if length(these_pre)>2
        %             %Violin plot
        %             [mean_out, CIout]=drgViolinPoint(these_pre...
        %                 ,edges,bar_offset,rand_offset,'k','k',3);
        %         end
        %
        %         glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_pre))=these_pre;
        %         glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_pre))=grNo*ones(1,length(these_pre));
        %         glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_pre))=1*ones(1,length(these_pre));
        %         glm_sig_ii=glm_sig_ii+length(these_pre);
        %
        %         id_sig_ii=id_sig_ii+1;
        %         input_sig_data(id_sig_ii).data=these_pre;
        %         input_sig_data(id_sig_ii).description=['Pre ' fr_per_names{grNo}];
        %
        %                 bar_offset=bar_offset+1;
        %         bar(bar_offset,mean(these_FV),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        %         if length(these_FV)>2
        %             %Violin plot
        %             [mean_out, CIout]=drgViolinPoint(these_FV...
        %                 ,edges,bar_offset,rand_offset,'k','k',3);
        %         end
        %
        %              glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_FV))=these_FV;
        %         glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_FV))=grNo*ones(1,length(these_FV));
        %         glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_FV))=1*ones(1,length(these_FV));
        %         glm_sig_ii=glm_sig_ii+length(these_FV);
        %
        %         id_sig_ii=id_sig_ii+1;
        %         input_sig_data(id_sig_ii).data=these_FV;
        %         input_sig_data(id_sig_ii).description=['FV ' fr_per_names{grNo}];


        bar_offset=bar_offset+1;
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',3);
        end

        glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
        glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
        glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
        glm_sig_ii=glm_sig_ii+length(these_odor);

        id_sig_ii=id_sig_ii+1;
        input_sig_data(id_sig_ii).data=these_odor;
        input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];

        bar_offset=bar_offset+1;
    end

    xticks([1 3 5])
    labels='xticklabels({';
    for ii_label=2:4
            labels=[labels '''' per_names{ii_label} ''', '];
    end
    labels=[labels(1:end-2) '})'];
    eval(labels)

    ylabel('% responsive')
    title('Percent ROIs responding to S+')

    %Perform the glm
    fprintf(1, ['\nglm for percent responsive to S+\n'])


    tbl = table(glm_sig.data',glm_sig.grNo',...
        'VariableNames',{'percent_divergent','percent_correct'});
    mdl = fitglm(tbl,'percent_divergent~percent_correct'...
        ,'CategoricalVars',[2])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);



    %Plot a bar graph showing percent responding to S-
    glm_sig=[];
    glm_sig_ii=0;
    input_sig_data=[];
    id_sig_ii=0;

    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

    hold on
    bar_offset=0;
    edges=[0:0.2:10];
    rand_offset=0.8;

    for grNo=1:3
        these_pre=[];
        these_FV=[];
        these_odor=[];

        for mouseNo=unique(handles_out.mouseNo)
            switch grNo
                case 1
                    files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                case 2
                    files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                case 3
                    files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
            end

            if sum(files_included)>0
                n_all=zeros(1,sum(files_included));
                n_all(1,:)=handles_out.output_data_N(files_included);
                n_odor=zeros(1,sum(files_included));
                n_odor(1,:)=handles_out.output_data_Sm(files_included);
                this_odor=n_odor./n_all;
                these_odor=[these_odor 100*mean(this_odor)];
            end

        end

        %         bar_offset=bar_offset+1;
        %         bar(bar_offset,mean(these_pre),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        %         if length(these_pre)>2
        %             %Violin plot
        %             [mean_out, CIout]=drgViolinPoint(these_pre...
        %                 ,edges,bar_offset,rand_offset,'k','k',3);
        %         end
        %
        %         glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_pre))=these_pre;
        %         glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_pre))=grNo*ones(1,length(these_pre));
        %         glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_pre))=1*ones(1,length(these_pre));
        %         glm_sig_ii=glm_sig_ii+length(these_pre);
        %
        %         id_sig_ii=id_sig_ii+1;
        %         input_sig_data(id_sig_ii).data=these_pre;
        %         input_sig_data(id_sig_ii).description=['Pre ' fr_per_names{grNo}];
        %
        %                 bar_offset=bar_offset+1;
        %         bar(bar_offset,mean(these_FV),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        %         if length(these_FV)>2
        %             %Violin plot
        %             [mean_out, CIout]=drgViolinPoint(these_FV...
        %                 ,edges,bar_offset,rand_offset,'k','k',3);
        %         end
        %
        %              glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_FV))=these_FV;
        %         glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_FV))=grNo*ones(1,length(these_FV));
        %         glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_FV))=1*ones(1,length(these_FV));
        %         glm_sig_ii=glm_sig_ii+length(these_FV);
        %
        %         id_sig_ii=id_sig_ii+1;
        %         input_sig_data(id_sig_ii).data=these_FV;
        %         input_sig_data(id_sig_ii).description=['FV ' fr_per_names{grNo}];


        bar_offset=bar_offset+1;
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',3);
        end

        glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
        glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
        glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
        glm_sig_ii=glm_sig_ii+length(these_odor);

        id_sig_ii=id_sig_ii+1;
        input_sig_data(id_sig_ii).data=these_odor;
        input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];

        bar_offset=bar_offset+1;
    end

    xticks([1 3 5])
    labels='xticklabels({';
    for ii_label=2:4
        labels=[labels '''' per_names{ii_label} ''', '];
    end
    labels=[labels(1:end-2) '})'];
    eval(labels)

    ylabel('% responsive')
    title('Percent ROIs responding to S-')

    %Perform the glm
    fprintf(1, ['\nglm for percent responsive to S-\n'])


    tbl = table(glm_sig.data',glm_sig.grNo',...
        'VariableNames',{'percent_divergent','percent_correct'});
    mdl = fitglm(tbl,'percent_divergent~percent_correct'...
        ,'CategoricalVars',[2])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    pffft=1;





    %Calculate the crosscorrelations
    croscorr_traces=abs(corrcoef(all_dFFspm)); %please note that I am using the absolute value

    %Set autocorrelations to zero
    for ii=1:size(croscorr_traces,1)
        croscorr_traces(ii,ii)=0;
    end
    Z = linkage(croscorr_traces,'complete','correlation');
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);
    [H,T,outperm]=dendrogram(Z,0,'Orientation','left');
    set(H,'LineWidth',2)
    hFig=figure(figureNo);
    set(hFig, 'units','normalized','position',[.05 .1 .12 .8])

    %re-sort the matrix
    for ii=1:size(croscorr_traces,1)
        for jj=1:size(croscorr_traces,1)
            perm_croscorr_traces(ii,jj)=croscorr_traces(outperm(ii),outperm(jj));
        end
    end



    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.15 .1 .6 .8])
    hold on
    pcolor(perm_croscorr_traces)
    colormap fire

    caxis([0    1])
    title(['Cross correlations for all ROIs'])
    xlim([1 ii_dFF])
    ylim([1 ii_dFF])

    %Plot rainbow
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


    prain=[0:0.6/99:0.6];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    %Now plot the timecourses
    for grNo=1:3

        switch grNo
            case 1
                ROIs_included=(all_dFFspm_pcorr>=45)&(all_dFFspm_pcorr<=65);
            case 2
                ROIs_included=(all_dFFspm_pcorr>65)&(all_dFFspm_pcorr<80);
            case 3
                ROIs_included=(all_dFFspm_pcorr>=80);
        end
        %S+
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
        hold on

        

        sorted_all_dFFsplus=[];
        ii_included=0;
        for ii=1:ii_dFF
            if ROIs_included(ii)
                ii_included=ii_included+1;
            sorted_all_dFFsplus(ii_included,:)=all_dFFsplus(outperm(ii),:);
            end
        end

        time_span_mat=repmat(time_span,ii_included,1);
        ROI_mat=repmat(1:ii_included,length(time_span),1)';

        pcolor(time_span_mat,ROI_mat,sorted_all_dFFsplus)
        colormap fire
        shading flat

        caxis([prctile(all_dFFspm(:),1) prctile(all_dFFspm(:),99)])

        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],[0 ii_dFF],'-r')
        end

        xlim([-7 15])
        ylim([1 ii_included])
        title(['S+ ' per_names{grNo+1}])
        xlabel('Time (sec)')
        ylabel('ROI number')

        %S-
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
        hold on

       

        sorted_all_dFFsminus=[];
        ii_included=0;
        for ii=1:ii_dFF
            if ROIs_included(ii)
                ii_included=ii_included+1;
            sorted_all_dFFsminus(ii_included,:)=all_dFFsminus(outperm(ii),:);
            end
        end

        time_span_mat=repmat(time_span,ii_included,1);
        ROI_mat=repmat(1:ii_included,length(time_span),1)';

        pcolor(time_span_mat,ROI_mat,sorted_all_dFFsminus)
        colormap fire
        shading flat

        caxis([prctile(all_dFFspm(:),1) prctile(all_dFFspm(:),99)])

        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],[0 ii_dFF],'-r')
        end

        xlim([-7 15])
        ylim([1 ii_included])
        title(['S- ' per_names{grNo+1}])
        xlabel('Time (sec)')
        ylabel('ROI number')

    end

    %Plot rainbow
    figureNo=figureNo+1;
    try
        close(figureNo)
    catch
    end

    hFig = figure(figureNo);

    set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


    prain=[prctile(all_dFFspm(:),1):(prctile(all_dFFspm(:),99)-prctile(all_dFFspm(:),1))/99:prctile(all_dFFspm(:),99)];
    pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
    %             colormap jet
    colormap fire
    shading interp
    ax=gca;
    set(ax,'XTickLabel','')

    pffft=1;

    %Uncomment this if you want to browse through the figures
    tic
    for figNo=1:figureNo
        figure(figNo)
        this_toc=toc;
        while toc-this_toc<4
        end
    end
    pffft=1;
end
