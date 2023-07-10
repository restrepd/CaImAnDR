function handles_out=drgCaImAn_pval_batch(handles_choices)
%This program performs a survey of ROI responsiveness to S+ and S-
% and divergent responses to S+ and S-
%
% if show_figures==1 the user can save the responsive/divergent 
% ROI timecourse figures when he/she presses "k" in the keyboard
% 
% the program surveys ROI cross correlations and evaluates 
% the hierarchical clustering of the ROIs using the maltab linkage function
% with unweighted average euclidean distance (UPGMA)
% the input a pval handles_choices file with the pre_per file location and all
% other choices
%

close all
clear all

handles_out=[];



if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_pval_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_pval_batch run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

new_no_files=handles.no_files;

if isfield(handles,'group_algo')
    group_algo=handles.group_algo;
else
    group_algo=1;
end

%handles.group_algo
% 1 Ming's spm with forward and reversed merged
% 2 Fabio's two odor application with spm sequence

if isfield(handles,'suffix_out')
    suffix_out=handles.suffix_out;
else
    suffix_out='_pval.mat';
end

if isfield(handles,'first_file')
    first_file=handles.first_file;
else
    first_file=1;
end

if isfield(handles,'show_figures')
    show_figures=handles.show_figures;
else
    show_figures=0;
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

pre_t=[-3 -2];
fv_t=[-1 0];
odor_t=[0.5 5.5];

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents = [0 delta_odor delta_odor_on_reinf_on delta_reinf+delta_odor_on_reinf_on];

%We performed controls and we do not need a time shift
% t_shift=0;


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

    handles_out.all_div_dFFsplus=[];
    handles_out.all_div_dFFsminus=[];
    handles_out.all_div_dFFspm=[];
    handles_out.all_div_ii_dFF=0;
    handles_out.all_div_t=[];
    handles_out.all_div_delta_dFFsplus=[];
    handles_out.all_div_delta_dFFsminus=[];
    handles_out.all_div_dFFspm_group=[];
    handles_out.all_div_dFFspm_pcorr=[];

    handles_out.all_spresp_ii_dFF=0;
    handles_out.all_spresp_dFFspm_group=[];
    handles_out.all_spresp_dFFspm_pcorr=[];
    handles_out.all_spresp_dFFsplus=[];
    handles_out.all_spresp_dFFsminus=[];
    handles_out.all_spresp_delta_dFFsplus=[];
    handles_out.all_spresp_delta_dFFsminus=[];
    handles_out.all_spresp_dFFspm=[];
    handles_out.all_spresp_dFFspm=[];

    handles_out.all_smresp_ii_dFF=0;
    handles_out.all_smresp_dFFspm_group=[];
    handles_out.all_smresp_dFFspm_pcorr=[];
    handles_out.all_smresp_dFFsplus=[];
    handles_out.all_smresp_dFFsminus=[];
    handles_out.all_smresp_delta_dFFsplus=[];
    handles_out.all_smresp_delta_dFFsminus=[];
    handles_out.all_smresp_dFFspm=[];
    handles_out.all_smresp_dFFspm=[];


    for fileNo=first_file:length(handles.FileName_pre_per)
        tic

        fprintf(1, ['Processing file No  ' num2str(fileNo) '\n']);

        pre_per_PathName=handles.PathName_pre_per{fileNo};
        pre_per_FileName=handles.FileName_pre_per{fileNo};

        handles_choices.pre_per_PathName=pre_per_PathName;
        handles_choices.pre_per_FileName=pre_per_FileName;

        handles_choices.show_figures=0; %Important, do not show the traces, otherwise the program gets bugged down....

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
    cd(handles.choiceBatchPathName)

    %If the user wants figures saved remove the figures directory and make a new figures directory
    if show_figures==1
        try
            rmdir('figures','s')
        catch
        end
        mkdir(['figures' suffix_out])
    end

    ii_included=0;

    %Note: These data were acquired at different rates. Here they are all
    %resampled to a dt of 0.03 sec
    dt_res=0.03;
    time_span=t_from:dt_res:t_to;
    handles_out.time_span=time_span;
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

%                     %Now do the time shift
%                     no_shift=floor(t_shift/(time_span(2)-time_span(1)));
% 
%                     if all_handles(fileNo).handles_out.shift_time==1
%                         if no_shift>0
%                             these_FR(no_shift+1:end)=these_FR(1:end-no_shift);
%                             these_FR(1:no_shift)=these_FR(no_shift);
%                         else
%                             these_FR(1:end-no_shift)=these_FR(no_shift+1:end);
%                             these_FR(1:end-no_shift+1:end)=these_FR(no_shift);
%                         end
%                     end


                    firingRates(n,s,1:T,trNo)=these_FR;
                    trialNum(n,s)=all_trialNum(n,s);
                end
            end

        end




        %Now do p-vals


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

        first_fig=1;
        

        %Show the figure for divergent dFF in the odor period if p<pFDR
        %Save the dFFSplus and dFFSminus

        for ii_p=1:length(handles_out.file(fileNo).output_data_odor.p)
            this_p=handles_out.file(fileNo).output_data_odor.p(ii_p);
            min_ii=ii_p;
            if this_p<handles_out.file(fileNo).output_data_odor.pFDR

                 %get the dF/F
                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

                if show_figures==1

                    try
                        close(first_fig)
                    catch
                    end
                    hFig=figure(first_fig);

                    ax=gca;ax.LineWidth=3;
                    set(hFig, 'units','normalized','position',[.2 .2 .25 .25])


                    hold on
                end

               

                try
                    %S-
                    CIpvsm = bootci(1000, @mean, dFFsminus');
                    meanpvsm=mean(dFFsminus',1);
                    CIpvsm(1,:)=meanpvsm-CIpvsm(1,:);
                    CIpvsm(2,:)=CIpvsm(2,:)-meanpvsm;

                    if show_figures==1
                        [hlpvl, hppvl] = boundedline(time_span,mean(dFFsminus'), CIpvsm','cmap',[158/255 31/255 99/255]);
                    end

                    %S+
                    CIpvsp = bootci(1000, @mean, dFFsplus');
                    meanpvsp=mean(dFFsplus',1);
                    CIpvsp(1,:)=meanpvsp-CIpvsp(1,:);
                    CIpvsp(2,:)=CIpvsp(2,:)-meanpvsp;

                    if show_figures==1
                        [hlpvl, hppvl] = boundedline(time_span, mean(dFFsplus'), CIpvsp','cmap',[0 114/255 178/255]);


                        plot(time_span',mean(dFFsminus')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
                        plot(time_span',mean(dFFsplus')','Color',[0 114/255 178/255],'LineWidth',1.5);

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


                        switch handles.save_all_figs
                            case 1
                                savefig([handles.choiceBatchPathName 'figures' suffix_out '\odor_f' num2str(fileNo) 'n' num2str(min_ii)])
                            case 2
                                answer = questdlg('Do you want to save the figure?', ...
                                    'Save figure', ...
                                    'Yes','No','No');
                                % Handle response
                                switch answer
                                    case 'Yes'
                                        savefig([handles.choiceBatchPathName 'figures' suffix_out '\odor_f' num2str(fileNo) 'n' num2str(min_ii)])
                                end
                        end
                    end

                    handles_out.all_div_ii_dFF=handles_out.all_div_ii_dFF+1;
                    handles_out.all_div_dFFspm_group(handles_out.all_div_ii_dFF)=handles.group(fileNo);
                    handles_out.all_div_dFFspm_pcorr(handles_out.all_div_ii_dFF)=handles_out.perCorr(fileNo);

                    handles_out.all_div_dFFsplus(handles_out.all_div_ii_dFF,1:length(time_span))=mean(dFFsplus')';
                    handles_out.all_div_dFFsminus(handles_out.all_div_ii_dFF,1:length(time_span))=mean(dFFsminus')';

                    handles_out.all_div_delta_dFFsplus(handles_out.all_div_ii_dFF)=mean(handles_out.all_div_dFFsplus(handles_out.all_div_ii_dFF,(time_span>=odor_t(1))&(time_span<=odor_t(2))))...
                        -mean(handles_out.all_div_dFFsplus(handles_out.all_div_ii_dFF,(time_span>=pre_t(1))&(time_span<=pre_t(2))));
                    handles_out.all_div_delta_dFFsminus(handles_out.all_div_ii_dFF)=mean(handles_out.all_div_dFFsminus(handles_out.all_div_ii_dFF,(time_span>=odor_t(1))&(time_span<=odor_t(2))))...
                        -mean(handles_out.all_div_dFFsminus(handles_out.all_div_ii_dFF,(time_span>=pre_t(1))&(time_span<=pre_t(2))));

                    handles_out.all_div_dFFspm(1:length(time_span),handles_out.all_div_ii_dFF)=mean(dFFsplus')';
                    handles_out.all_div_dFFspm(length(time_span)+1:2*length(time_span),handles_out.all_div_ii_dFF)=mean(dFFsminus')';

                    %First time for divergence
                    all_done=0;
                    ii=find(time_span>=pre_t(2),1,'first');
                    dt_required=0.5;
                    ii_required=ceil(dt_required/(time_span(2)-time_span(1)));
                    ii_included=0;
                    while all_done==0

                        if handles_out.all_div_delta_dFFsplus(handles_out.all_div_ii_dFF)-handles_out.all_div_delta_dFFsminus(handles_out.all_div_ii_dFF)>0
                            handles_out.all_div_ii_dFF_div=find((meanpvsp(ii:end)-CIpvsp(1,ii:end))>(meanpvsm(ii:end)+CIpvsm(2,ii:end)),1,'first');
                            if (~isempty(handles_out.all_div_ii_dFF_div))&(ii+handles_out.all_div_ii_dFF_div-1+ii_required<length(time_span))
                                if sum((meanpvsp(ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required)-CIpvsp(1,ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required))...
                                        >(meanpvsm(ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required)+CIpvsm(2,ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required)))...
                                        ==length(meanpvsm(ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required))
                                    all_done=1;

                                end
                                ii=ii+handles_out.all_div_ii_dFF_div;
                            else
                                all_done=1;
                                ii=NaN;
                            end

                        else
                            handles_out.all_div_ii_dFF_div=find((meanpvsm(ii:end)-CIpvsm(1,ii:end))>(meanpvsp(ii:end)+CIpvsp(2,ii:end)),1,'first');
                            if (~isempty(handles_out.all_div_ii_dFF_div))&(ii+handles_out.all_div_ii_dFF_div-1+ii_required<length(time_span))
                                if sum((meanpvsm(ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required)-CIpvsm(1,ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required))...
                                        >(meanpvsp(ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required)+CIpvsp(2,ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required)))...
                                        ==length(meanpvsm(ii+handles_out.all_div_ii_dFF_div-1:ii+handles_out.all_div_ii_dFF_div-1+ii_required))
                                    all_done=1;

                                end
                                ii=ii+handles_out.all_div_ii_dFF_div;
                            else
                                all_done=1;
                                ii=NaN;
                            end

                        end

                    end
                    handles_out.all_div_t(handles_out.all_div_ii_dFF)=time_span(ii);

                catch
                end


            end
        end

        %         %Plot the timecourse for the lowest p value for the
        %         %FV period
        %         [min_p min_ii]=min(handles_out.file(fileNo).output_data_FV.p);
        %
        %         %Show the figure only if min_p<pFDR
        %         if min_p<handles_out.file(fileNo).output_data_FV.pFDR
        %             if show_figures==1
        %                 figureNo = figureNo + 1;
        %                 try
        %                     close(figureNo)
        %                 catch
        %                 end
        %                 hFig=figure(figureNo);
        %
        %                 ax=gca;ax.LineWidth=3;
        %                 set(hFig, 'units','normalized','position',[.2 .2 .25 .25])
        %
        %
        %                 hold on
        %
        %                 %get the dF/F
        %
        %                 dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
        %                 dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
        %                 dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
        %                 dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));
        %
        %                 try
        %                     %S-
        %                     CIpv = bootci(1000, @mean, dFFsminus');
        %                     meanpv=mean(dFFsminus',1);
        %                     CIpv(1,:)=meanpv-CIpv(1,:);
        %                     CIpv(2,:)=CIpv(2,:)-meanpv;
        %
        %
        %                     [hlpvl, hppvl] = boundedline(time_span,mean(dFFsminus'), CIpv','cmap',[158/255 31/255 99/255]);
        %
        %                     %S+
        %                     CIpv = bootci(1000, @mean, dFFsplus');
        %                     meanpv=mean(dFFsplus',1);
        %                     CIpv(1,:)=meanpv-CIpv(1,:);
        %                     CIpv(2,:)=CIpv(2,:)-meanpv;
        %
        %
        %                     [hlpvl, hppvl] = boundedline(time_span, mean(dFFsplus'), CIpv','cmap',[0 114/255 178/255]);
        %
        %
        %                     plot(time_span',mean(dFFsminus')','Color',[158/255 31/255 99/255],'LineWidth',1.5);
        %                     plot(time_span',mean(dFFsplus')','Color',[0 114/255 178/255],'LineWidth',1.5);
        %
        %
        %
        %                 catch
        %                 end
        %
        %                 this_ylim=ylim;
        %                 for ii_te=1:length(timeEvents)
        %                     plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
        %                 end
        %
        %                 text(-5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
        %                 text(-5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)
        %
        %                 xlim([-7 15])
        %
        %
        %                 xlabel('Time(sec)')
        %                 ylabel('dFF')
        %
        %                 if all_handles(fileNo).handles_out.shift_time==1
        %                     title(['dFF FV n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo)) ' tsh'])
        %                 else
        %                     title(['dFF FV n= ' num2str(min_ii), ' fileNo= ' num2str(fileNo) ' pCorr= ' num2str(pCorr_per_file(fileNo))])
        %                 end
        %
        %                 savefig([handles.choiceBatchPathName 'figures\' 'odor_f' num2str(fileNo) 'n' num2str(min_ii)])
        %             end
        %         end

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

        %         [min_p min_ii]=min(handles_out.file(fileNo).output_data_odor_Sp.p);

        for ii_p=1:length(handles_out.file(fileNo).output_data_odor_Sp.p)
            this_p=handles_out.file(fileNo).output_data_odor_Sp.p(ii_p);
            min_ii=ii_p;
            if this_p<handles_out.file(fileNo).output_data_odor_Sp.pFDR
                %         if min_p<handles_out.file(fileNo).output_data_odor_Sp.pFDR
                %get the dF/F
                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

                if show_figures==1

                    try
                        close(first_fig+1)
                    catch
                    end
                    hFig=figure(first_fig+1);

                    ax=gca;ax.LineWidth=3;
                    set(hFig, 'units','normalized','position',[.3 .2 .25 .25])


                    hold on



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

                    switch handles.save_all_figs
                        case 1
                            savefig([handles.choiceBatchPathName 'figures' suffix_out '\Sp_f' num2str(fileNo) 'n' num2str(min_ii)])
                        case 2
                            answer = questdlg('Do you want to save the figure?', ...
                                'Save figure', ...
                                'Yes','No','No');
                            % Handle response
                            switch answer
                                case 'Yes'
                                    savefig([handles.choiceBatchPathName 'figures' suffix_out '\Sp_f' num2str(fileNo) 'n' num2str(min_ii)])
                            end
                    end


                end

                handles_out.all_spresp_ii_dFF=handles_out.all_spresp_ii_dFF+1;
                handles_out.all_spresp_dFFspm_group(handles_out.all_spresp_ii_dFF)=handles.group(fileNo);
                handles_out.all_spresp_dFFspm_pcorr(handles_out.all_spresp_ii_dFF)=handles_out.perCorr(fileNo);

                handles_out.all_spresp_dFFsplus(handles_out.all_spresp_ii_dFF,1:length(time_span))=mean(dFFsplus')';
                handles_out.all_spresp_dFFsminus(handles_out.all_spresp_ii_dFF,1:length(time_span))=mean(dFFsminus')';

                handles_out.all_spresp_delta_dFFsplus(handles_out.all_spresp_ii_dFF)=mean(handles_out.all_spresp_dFFsplus(handles_out.all_spresp_ii_dFF,(time_span>=odor_t(1))&(time_span<=odor_t(2))))...
                    -mean(handles_out.all_spresp_dFFsplus(handles_out.all_spresp_ii_dFF,(time_span>=pre_t(1))&(time_span<=pre_t(2))));
                handles_out.all_spresp_delta_dFFsminus(handles_out.all_spresp_ii_dFF)=mean(handles_out.all_spresp_dFFsminus(handles_out.all_spresp_ii_dFF,(time_span>=odor_t(1))&(time_span<=odor_t(2))))...
                    -mean(handles_out.all_spresp_dFFsminus(handles_out.all_spresp_ii_dFF,(time_span>=pre_t(1))&(time_span<=pre_t(2))));

                handles_out.all_spresp_dFFspm(1:length(time_span),handles_out.all_spresp_ii_dFF)=mean(dFFsplus')';
                handles_out.all_spresp_dFFspm(length(time_span)+1:2*length(time_span),handles_out.all_spresp_ii_dFF)=mean(dFFsminus')';

            end
        end

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
        for ii_p=1:length(handles_out.file(fileNo).output_data_odor_Sm.p)
            this_p=handles_out.file(fileNo).output_data_odor_Sm.p(ii_p);
            min_ii=ii_p;
            %Show the figure only if min_p<pFDR
            if this_p<handles_out.file(fileNo).output_data_odor_Sm.pFDR


                %get the dF/F

                dFFsplus=zeros(length(time_span),trialNum(min_ii,1));
                dFFsplus(:,:)=firingRates(min_ii,1,:,1:trialNum(min_ii,1));
                dFFsminus=zeros(length(time_span),trialNum(min_ii,2));
                dFFsminus(:,:)=firingRates(min_ii,2,:,1:trialNum(min_ii,2));

                if show_figures==1

                    try
                        close(first_fig+2)
                    catch
                    end
                    hFig=figure(first_fig+2);

                    ax=gca;ax.LineWidth=3;
                    set(hFig, 'units','normalized','position',[.2 .4 .25 .25])


                    hold on



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

                    switch handles.save_all_figs
                        case 1
                            savefig([handles.choiceBatchPathName 'figures' suffix_out '\Sm_f' num2str(fileNo) 'n' num2str(min_ii)])
                        case 2
                            answer = questdlg('Do you want to save the figure?', ...
                                'Save figure', ...
                                'Yes','No','No');
                            % Handle response
                            switch answer
                                case 'Yes'
                                    savefig([handles.choiceBatchPathName 'figures' suffix_out '\Sm_f' num2str(fileNo) 'n' num2str(min_ii)])
                            end
                    end


                end

                handles_out.all_smresp_ii_dFF=handles_out.all_smresp_ii_dFF+1;
                handles_out.all_smresp_dFFspm_group(handles_out.all_smresp_ii_dFF)=handles.group(fileNo);
                handles_out.all_smresp_dFFspm_pcorr(handles_out.all_smresp_ii_dFF)=handles_out.perCorr(fileNo);

                handles_out.all_smresp_dFFsplus(handles_out.all_smresp_ii_dFF,1:length(time_span))=mean(dFFsplus')';
                handles_out.all_smresp_dFFsminus(handles_out.all_smresp_ii_dFF,1:length(time_span))=mean(dFFsminus')';

                handles_out.all_smresp_delta_dFFsplus(handles_out.all_smresp_ii_dFF)=mean(handles_out.all_smresp_dFFsplus(handles_out.all_smresp_ii_dFF,(time_span>=odor_t(1))&(time_span<=odor_t(2))))...
                    -mean(handles_out.all_smresp_dFFsplus(handles_out.all_smresp_ii_dFF,(time_span>=pre_t(1))&(time_span<=pre_t(2))));
                handles_out.all_smresp_delta_dFFsminus(handles_out.all_smresp_ii_dFF)=mean(handles_out.all_smresp_dFFsminus(handles_out.all_smresp_ii_dFF,(time_span>=odor_t(1))&(time_span<=odor_t(2))))...
                    -mean(handles_out.all_smresp_dFFsminus(handles_out.all_smresp_ii_dFF,(time_span>=pre_t(1))&(time_span<=pre_t(2))));

                handles_out.all_smresp_dFFspm(1:length(time_span),handles_out.all_smresp_ii_dFF)=mean(dFFsplus')';
                handles_out.all_smresp_dFFspm(length(time_span)+1:2*length(time_span),handles_out.all_smresp_ii_dFF)=mean(dFFsminus')';
            end

        end
    end

    %Plot a bar graph showing percent signifcant divergence
    glm_sig=[];
    glm_sig_ii=0;
    input_sig_data=[];
    id_sig_ii=0;

    figureNo=0;
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    hold on
    bar_offset=0;
    edges=[0:0.2:10];
    rand_offset=0.8;

    switch handles.group_algo
        case 1
            %Ming
            groups=[1:3];
        case 2
            %Fabio
            groups=unique(handles.group);
    end

    for grNo=groups
        these_pre=[];
        these_FV=[];
        these_odor=[];

        for mouseNo=unique(handles_out.mouseNo)

            switch handles.group_algo
                case 1
                    %Ming
                    switch grNo
                        case 1
                            files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                        case 2
                            files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                        case 3
                            files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
                    end
                case 2
                    %Fabio
                    files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo);
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

        if ~isempty(these_odor)
%             bar_offset=bar_offset+1;
            bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            if length(these_odor)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(these_odor...
                    ,edges,bar_offset,rand_offset,'k','k',5);
            end

            glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
            glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
            glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
            glm_sig_ii=glm_sig_ii+length(these_odor);

            id_sig_ii=id_sig_ii+1;
            input_sig_data(id_sig_ii).data=these_odor;
            switch handles.group_algo
                case 1
                    input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];
                case 2
                    input_sig_data(id_sig_ii).description=['Odor ' handles.group_names{grNo}];
            end
        end

        bar_offset=bar_offset+1;
    end

    switch handles.group_algo
        case 1
            %Ming

            xticks([0 1 2])
            labels='xticklabels({';
            for ii_label=2:4
                labels=[labels '''' per_names{ii_label} ''', '];
            end
            labels=[labels(1:end-2) '})'];
            eval(labels)
            xlim([-0.7 2.7])
        case 2
            %Fabio
            xticks(2*groups-1)
            labels='xticklabels({';
            for ii_label=1:length(groups)
                labels=[labels '''' handles.group_names{ii_label} ''', '];
            end
            labels=[labels(1:end-2) '})'];
            eval(labels)
    end


    ylabel('% divergent')
    title('Percent divergent ROIs per mouse')

    %Perform the glm
    fprintf(1, ['\nglm for percent divergent per mouse\n'])


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
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    hold on
    bar_offset=0;
    edges=[0:0.2:10];
    rand_offset=0.8;

    for grNo=groups
        these_pre=[];
        these_FV=[];
        these_odor=[];

        for mouseNo=unique(handles_out.mouseNo)

            switch handles.group_algo
                case 1
                    %Ming
                    switch grNo
                        case 1
                            files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                        case 2
                            files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                        case 3
                            files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
                    end
                case 2
                    %Fabio
                    files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo);
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

        if ~isempty(these_odor)
%             bar_offset=bar_offset+1;
            bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
            if length(these_odor)>2
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(these_odor...
                    ,edges,bar_offset,rand_offset,'k','k',5);
            end

            glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
            glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
            glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
            glm_sig_ii=glm_sig_ii+length(these_odor);

            id_sig_ii=id_sig_ii+1;
            input_sig_data(id_sig_ii).data=these_odor;

            switch handles.group_algo
                case 1
                    input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];
                case 2
                    input_sig_data(id_sig_ii).description=['Odor ' handles.group_names{grNo}];
            end
        end

        bar_offset=bar_offset+1;
    end

    switch handles.group_algo
        case 1
            %Ming

            xticks([0 1 2])
            labels='xticklabels({';
            for ii_label=2:4
                labels=[labels '''' per_names{ii_label} ''', '];
            end
            labels=[labels(1:end-2) '})'];
            eval(labels)
            xlim([-0.7 2.7])
        case 2
            %Fabio
            xticks(2*groups-1)
            labels='xticklabels({';
            for ii_label=1:length(groups)
                labels=[labels '''' handles.group_names{ii_label} ''', '];
            end
            labels=[labels(1:end-2) '})'];
            eval(labels)
    end

    ylabel('% responsive')
    title('Percent ROIs responding to S+ per mouse')

    %Perform the glm
    fprintf(1, ['\nglm for percent responsive to S+ per mouse\n'])


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
    set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

    hold on
    bar_offset=0;
    edges=[0:0.2:10];
    rand_offset=0.8;

    for grNo=groups
        these_pre=[];
        these_FV=[];
        these_odor=[];

        for mouseNo=unique(handles_out.mouseNo)
            switch handles.group_algo
                case 1
                    %Ming
                    switch grNo
                        case 1
                            files_included=(handles_out.perCorr>=45)&(handles_out.perCorr<=65)&(handles_out.mouseNo==mouseNo);
                        case 2
                            files_included=(handles_out.perCorr>65)&(handles_out.perCorr<80)&(handles_out.mouseNo==mouseNo);
                        case 3
                            files_included=(handles_out.perCorr>=80)&(handles_out.mouseNo==mouseNo);
                    end
                case 2
                    %Fabio
                    files_included=(handles.group==grNo)&(handles_out.mouseNo==mouseNo);
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



        bar_offset=bar_offset+1;
        bar(bar_offset,mean(these_odor),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
        if length(these_odor)>2
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(these_odor...
                ,edges,bar_offset,rand_offset,'k','k',5);
        end

        glm_sig.data(glm_sig_ii+1:glm_sig_ii+length(these_odor))=these_odor;
        glm_sig.grNo(glm_sig_ii+1:glm_sig_ii+length(these_odor))=grNo*ones(1,length(these_odor));
        glm_sig.epoch(glm_sig_ii+1:glm_sig_ii+length(these_odor))=1*ones(1,length(these_odor));
        glm_sig_ii=glm_sig_ii+length(these_odor);

        id_sig_ii=id_sig_ii+1;
        input_sig_data(id_sig_ii).data=these_odor;
        switch handles.group_algo
            case 1
                input_sig_data(id_sig_ii).description=['Odor ' fr_per_names{grNo}];
            case 2
                input_sig_data(id_sig_ii).description=['Odor ' handles.group_names{grNo}];
        end

        bar_offset=bar_offset+1;
    end

    switch handles.group_algo
        case 1
            %Ming

            xticks([1 3 5])
            labels='xticklabels({';
            for ii_label=2:4
                labels=[labels '''' per_names{ii_label} ''', '];
            end
            labels=[labels(1:end-2) '})'];
            eval(labels)
            xlim([-0.7 2.7])
        case 2
            %Fabio
            xticks(2*groups-1)
            labels='xticklabels({';
            for ii_label=1:length(groups)
                labels=[labels '''' handles.group_names{ii_label} ''', '];
            end
            labels=[labels(1:end-2) '})'];
            eval(labels)
    end

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



    %Now generate the summary figures for the divergent responses

    if handles_out.all_div_ii_dFF>5
        %Calculate the crosscorrelations
        %     croscorr_traces=abs(corrcoef(handles_out.all_div_dFFspm)); %please note that I am using the absolute value
        croscorr_traces=corrcoef(handles_out.all_div_dFFspm);

        %Set autocorrelations to zero
        for ii=1:size(croscorr_traces,1)
            croscorr_traces(ii,ii)=0;
        end
        Z = linkage(croscorr_traces,'complete','correlation');
        no_clusters=4;
        clusters = cluster(Z,'Maxclust',no_clusters);
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);
        %Do cutoff for 4 clusters
        cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
        [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
        set(H,'LineWidth',2)
        hFig=figure(figureNo);
        set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

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
        shading flat

        caxis([-1  1])
        title(['Cross correlations for all ROIs'])
        xlim([1 handles_out.all_div_ii_dFF])
        ylim([1 handles_out.all_div_ii_dFF])

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

        %Plot timecourses for all ROIs


        %S+
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
        hold on



        sorted_handles_out.all_div_dFFsplus=[];
        for ii=1:handles_out.all_div_ii_dFF
            sorted_handles_out.all_div_dFFsplus(ii,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
        end

        time_span_mat=repmat(time_span,ii,1);
        ROI_mat=repmat(1:ii,length(time_span),1)';

        pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
        colormap fire
        shading flat

        caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
        end

        xlim([-7 15])
        ylim([1 ii])
        title(['S+ for all ROIs'])
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


        sorted_handles_out.all_div_dFFsminus=[];
        for ii=1:handles_out.all_div_ii_dFF
            sorted_handles_out.all_div_dFFsminus(ii,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
        end

        time_span_mat=repmat(time_span,ii,1);
        ROI_mat=repmat(1:ii,length(time_span),1)';

        pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
        colormap fire
        shading flat

        caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
        end

        xlim([-7 15])
        ylim([1 ii])
        title(['S- for all ROIs'])
        xlabel('Time (sec)')
        ylabel('ROI number')

        %Plot the average timecourses per cluster
        for clus=1:max(clusters)
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);

            ax=gca;ax.LineWidth=3;
            set(hFig, 'units','normalized','position',[.3 .2 .3 .3])


            hold on

            %get the dF/F


            try

                %S-
                this_cluster_dFFsminus=[];
                ii_included=0;
                for ii=1:handles_out.all_div_ii_dFF
                    if clusters(ii)==clus
                        ii_included=ii_included+1;
                        this_cluster_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(ii,:);
                    end
                end


                CIpv = bootci(1000, @mean, this_cluster_dFFsminus);
                meanpv=mean(this_cluster_dFFsminus,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;


                [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_dFFsminus), CIpv','cmap',[158/255 31/255 99/255]);

                %S+
                this_cluster_dFFsplus=[];
                ii_included=0;
                for ii=1:handles_out.all_div_ii_dFF
                    if clusters(ii)==clus
                        ii_included=ii_included+1;
                        this_cluster_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(ii,:);
                    end
                end

                CIpv = bootci(1000, @mean, this_cluster_dFFsplus);
                meanpv=mean(this_cluster_dFFsplus,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;


                [hlpvl, hppvl] = boundedline(time_span, mean(this_cluster_dFFsplus), CIpv','cmap',[0 114/255 178/255]);


                plot(time_span',mean(this_cluster_dFFsminus)','Color',[158/255 31/255 99/255],'LineWidth',1.5);
                plot(time_span',mean(this_cluster_dFFsplus)','Color',[0 114/255 178/255],'LineWidth',1.5);



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


            title(['dFF divergent for cluster no ' num2str(clus)])

        end


        %Now plot the timecourses per group
        for grNo=1:3

            switch grNo
                case 1
                    ROIs_included=(handles_out.all_div_dFFspm_pcorr>=45)&(handles_out.all_div_dFFspm_pcorr<=65);
                case 2
                    ROIs_included=(handles_out.all_div_dFFspm_pcorr>65)&(handles_out.all_div_dFFspm_pcorr<80);
                case 3
                    ROIs_included=(handles_out.all_div_dFFspm_pcorr>=80);
            end

            if sum(ROIs_included)>0
                try
                    %S+
                    figureNo=figureNo+1;
                    try
                        close(figureNo)
                    catch
                    end

                    hFig = figure(figureNo);

                    set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
                    hold on



                    sorted_handles_out.all_div_dFFsplus=[];
                    ii_included=0;
                    for ii=1:handles_out.all_div_ii_dFF
                        if ROIs_included(ii)
                            ii_included=ii_included+1;
                            sorted_handles_out.all_div_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
                        end
                    end

                    time_span_mat=repmat(time_span,ii_included,1);
                    ROI_mat=repmat(1:ii_included,length(time_span),1)';

                    pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
                    colormap fire
                    shading flat

                    caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

                    for ii_te=1:length(timeEvents)
                        plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
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



                    sorted_handles_out.all_div_dFFsminus=[];
                    ii_included=0;
                    for ii=1:handles_out.all_div_ii_dFF
                        if ROIs_included(ii)
                            ii_included=ii_included+1;
                            sorted_handles_out.all_div_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
                        end
                    end

                    time_span_mat=repmat(time_span,ii_included,1);
                    ROI_mat=repmat(1:ii_included,length(time_span),1)';

                    pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
                    colormap fire
                    shading flat

                    caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])

                    for ii_te=1:length(timeEvents)
                        plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
                    end

                    xlim([-7 15])
                    ylim([1 ii_included])
                    title(['S- ' per_names{grNo+1}])
                    xlabel('Time (sec)')
                    ylabel('ROI number')

                catch
                end
            end
        end

        %Plot rainbow
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.49 .1 .05 .3])


        prain=[prctile(handles_out.all_div_dFFspm(:),1):(prctile(handles_out.all_div_dFFspm(:),99)-prctile(handles_out.all_div_dFFspm(:),1))/99:prctile(handles_out.all_div_dFFspm(:),99)];
        pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
        %             colormap jet
        colormap fire
        shading interp
        ax=gca;
        set(ax,'XTickLabel','')

        %Generate a histogram for the magnitude of the responses

        %     handles_out.all_div_t=[];
        %     handles_out.all_div_delta_dFFsplus=[];
        %     handles_out.all_div_delta_dFFsminus=[];

        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.2 .2 .4 .3])

        edges=-0.5:0.5:4;
        histogram(handles_out.all_div_delta_dFFsplus,edges,'FaceColor',[0 114/255 178/255])
        hold on
        histogram(handles_out.all_div_delta_dFFsminus,edges,'FaceColor',[158/255 31/255 99/255])

        ylabel('No of ROIs')
        xlabel('Delta zdFF')
        title('Histogram for Delta zdFF')

        this_ylim=ylim;

        text(2.5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
        text(2.5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)

        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

        plot(handles_out.all_div_delta_dFFsminus,handles_out.all_div_delta_dFFsplus,'ok')

        ylabel('S+')
        xlabel('S-')
        title('Delta zdFF')

        %Divergence time
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.2 .2 .4 .3])

        edges=[-1.75:0.5:5.25];
        histogram(handles_out.all_div_t(~isnan(handles_out.all_div_t)),edges)

        ylabel('No of ROIs')
        xlabel('Divegence time')
        title('Histogram for divergence time')

    end


    %Now generate the summary figures for the divergent responses

    if handles_out.all_spresp_ii_dFF+handles_out.all_smresp_ii_dFF>5
        %Calculate the crosscorrelations
        %     croscorr_traces=abs(corrcoef(handles_out.all_div_dFFspm)); %please note that I am using the absolute value
        all_spandmresp_dFFspm=[handles_out.all_spresp_dFFspm handles_out.all_smresp_dFFspm];
        all_spandmresp_iidFF=[1:handles_out.all_spresp_ii_dFF 1:handles_out.all_smresp_ii_dFF];
        all_spandmresp_sporsm=[ones(1,handles_out.all_spresp_ii_dFF) zeros(1,handles_out.all_smresp_ii_dFF)];
        croscorr_traces=corrcoef(all_spandmresp_dFFspm);

        %Set autocorrelations to zero
        for ii=1:size(croscorr_traces,1)
            croscorr_traces(ii,ii)=0;
        end
        Z = linkage(croscorr_traces,'complete','correlation');
        no_clusters=2; %This is set for Fabio's analysis
        clusters = cluster(Z,'Maxclust',no_clusters);
        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);
        %Do cutoff for 4 clusters
        cutoff = median([Z(end-(no_clusters-1),3) Z(end-(no_clusters-2),3)]);
        [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold',cutoff);
        set(H,'LineWidth',2)
        hFig=figure(figureNo);
        set(hFig, 'units','normalized','position',[.05 .1 .14 .8])

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
        shading flat

        caxis([-1    1])
        title(['Cross correlations for all ROIs'])
        xlim([1 handles_out.all_spresp_ii_dFF+handles_out.all_smresp_ii_dFF])
        ylim([1 handles_out.all_spresp_ii_dFF+handles_out.all_smresp_ii_dFF])

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

        %Plot timecourses for all ROIs for S+ and S-

        figureNo=figureNo+1;
        try
            close(figureNo)
        catch
        end

        hFig = figure(figureNo);

        set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
        hold on



        sorted_handles_out.all_spmresp_ii_dFF=0;
        sorted_handles_out.all_spmresp_sporsm=[];
        sorted_handles_out.all_sporsm_dFFspm=[];
        sorted_handles_out.clusters=[];
        for ii=1:length(all_spandmresp_iidFF)
            this_ii_dFF=all_spandmresp_iidFF(outperm(ii));
            this_sporsm=all_spandmresp_sporsm(outperm(ii));
            sorted_handles_out.all_spmresp_ii_dFF=sorted_handles_out.all_spmresp_ii_dFF+1;
            if this_sporsm==0
                %This is S-
                sorted_handles_out.all_sporsm_dFFspm(sorted_handles_out.all_spmresp_ii_dFF,:)=handles_out.all_smresp_dFFsminus(this_ii_dFF,:);
            else
                %This is S+
                sorted_handles_out.all_sporsm_dFFspm(sorted_handles_out.all_spmresp_ii_dFF,:)=handles_out.all_spresp_dFFsplus(this_ii_dFF,:);
            end
            sorted_handles_out.clusters(sorted_handles_out.all_spmresp_ii_dFF)=clusters(outperm(ii));
        end

        %pcolor does not show the first row
        pseudo_dFFspm=zeros(size(sorted_handles_out.all_sporsm_dFFspm,1)+1,size(sorted_handles_out.all_sporsm_dFFspm,2));
        pseudo_dFFspm(1:end-1,:)=sorted_handles_out.all_sporsm_dFFspm;

        time_span_mat=repmat(time_span,sorted_handles_out.all_spmresp_ii_dFF+1,1);
        ROI_mat=repmat(1:sorted_handles_out.all_spmresp_ii_dFF+1,length(time_span),1)';

        pcolor(time_span_mat,ROI_mat,pseudo_dFFspm)


%         time_span_mat=repmat(time_span,sorted_handles_out.all_spmresp_ii_dFF,1);
%         ROI_mat=repmat(1:sorted_handles_out.all_spmresp_ii_dFF,length(time_span),1)';
% 
%         pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_sporsm_dFFspm)

        colormap fire
        shading flat

        caxis([prctile(sorted_handles_out.all_sporsm_dFFspm(:),1) prctile(sorted_handles_out.all_sporsm_dFFspm(:),99)])

        for ii_te=1:length(timeEvents)
            plot([timeEvents(ii_te) timeEvents(ii_te)],[0 sorted_handles_out.all_spmresp_ii_dFF],'-r')
        end

        xlim([-7 15])
        ylim([1 ii])
        title(['Timecourses for S+ and S- with significant odor changes, all ROIs'])
        xlabel('Time (sec)')
        ylabel('ROI number')

       

        %Plot the average timecourses per cluster
        for clus=1:max(clusters)
            figureNo = figureNo + 1;
            try
                close(figureNo)
            catch
            end
            hFig=figure(figureNo);

            ax=gca;ax.LineWidth=3;
            set(hFig, 'units','normalized','position',[.3 .2 .3 .3])


            hold on

            %get the dF/F


            try

                
                this_cluster_dFF=[];
                ii_included=0;
                for ii=1:length(sorted_handles_out.clusters)
                    if sorted_handles_out.clusters(ii)==clus
                        ii_included=ii_included+1;
                        this_cluster_dFF(ii_included,:)=sorted_handles_out.all_sporsm_dFFspm(ii,:);
                    end
                end


                CIpv = bootci(1000, @mean, this_cluster_dFF);
                meanpv=mean(this_cluster_dFF,1);
                CIpv(1,:)=meanpv-CIpv(1,:);
                CIpv(2,:)=CIpv(2,:)-meanpv;


                [hlpvl, hppvl] = boundedline(time_span,mean(this_cluster_dFF), CIpv','cmap',[0 0 0]);

               

            catch
            end

            this_ylim=ylim;
            for ii_te=1:length(timeEvents)
                plot([timeEvents(ii_te) timeEvents(ii_te)],this_ylim,'-k')
            end

            xlim([-7 15])


            xlabel('Time(sec)')
            ylabel('dFF')


            title(['dFF for cluster no ' num2str(clus)])

        end


%         %Now plot the timecourses per group
%         for grNo=1:3
% 
%             switch grNo
%                 case 1
%                     ROIs_included=(handles_out.all_div_dFFspm_pcorr>=45)&(handles_out.all_div_dFFspm_pcorr<=65);
%                 case 2
%                     ROIs_included=(handles_out.all_div_dFFspm_pcorr>65)&(handles_out.all_div_dFFspm_pcorr<80);
%                 case 3
%                     ROIs_included=(handles_out.all_div_dFFspm_pcorr>=80);
%             end
% 
%             if sum(ROIs_included)>0
% 
%                 %S+
%                 figureNo=figureNo+1;
%                 try
%                     close(figureNo)
%                 catch
%                 end
% 
%                 hFig = figure(figureNo);
% 
%                 set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
%                 hold on
% 
% 
% 
%                 sorted_handles_out.all_div_dFFsplus=[];
%                 ii_included=0;
%                 for ii=1:handles_out.all_div_ii_dFF
%                     if ROIs_included(ii)
%                         ii_included=ii_included+1;
%                         sorted_handles_out.all_div_dFFsplus(ii_included,:)=handles_out.all_div_dFFsplus(outperm(ii),:);
%                     end
%                 end
% 
%                 time_span_mat=repmat(time_span,ii_included,1);
%                 ROI_mat=repmat(1:ii_included,length(time_span),1)';
% 
%                 pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsplus)
%                 colormap fire
%                 shading flat
% 
%                 caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])
% 
%                 for ii_te=1:length(timeEvents)
%                     plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
%                 end
% 
%                 xlim([-7 15])
%                 ylim([1 ii_included])
%                 title(['S+ ' per_names{grNo+1}])
%                 xlabel('Time (sec)')
%                 ylabel('ROI number')
% 
%                 %S-
%                 figureNo=figureNo+1;
%                 try
%                     close(figureNo)
%                 catch
%                 end
% 
%                 hFig = figure(figureNo);
% 
%                 set(hFig, 'units','normalized','position',[.05 .1 .3 .8])
%                 hold on
% 
% 
% 
%                 sorted_handles_out.all_div_dFFsminus=[];
%                 ii_included=0;
%                 for ii=1:handles_out.all_div_ii_dFF
%                     if ROIs_included(ii)
%                         ii_included=ii_included+1;
%                         sorted_handles_out.all_div_dFFsminus(ii_included,:)=handles_out.all_div_dFFsminus(outperm(ii),:);
%                     end
%                 end
% 
%                 time_span_mat=repmat(time_span,ii_included,1);
%                 ROI_mat=repmat(1:ii_included,length(time_span),1)';
% 
%                 pcolor(time_span_mat,ROI_mat,sorted_handles_out.all_div_dFFsminus)
%                 colormap fire
%                 shading flat
% 
%                 caxis([prctile(handles_out.all_div_dFFspm(:),1) prctile(handles_out.all_div_dFFspm(:),99)])
% 
%                 for ii_te=1:length(timeEvents)
%                     plot([timeEvents(ii_te) timeEvents(ii_te)],[0 handles_out.all_div_ii_dFF],'-r')
%                 end
% 
%                 xlim([-7 15])
%                 ylim([1 ii_included])
%                 title(['S- ' per_names{grNo+1}])
%                 xlabel('Time (sec)')
%                 ylabel('ROI number')
%             end
%         end
% 
%         %Plot rainbow
%         figureNo=figureNo+1;
%         try
%             close(figureNo)
%         catch
%         end
% 
%         hFig = figure(figureNo);
% 
%         set(hFig, 'units','normalized','position',[.49 .1 .05 .3])
% 
% 
%         prain=[prctile(handles_out.all_div_dFFspm(:),1):(prctile(handles_out.all_div_dFFspm(:),99)-prctile(handles_out.all_div_dFFspm(:),1))/99:prctile(handles_out.all_div_dFFspm(:),99)];
%         pcolor(repmat([1:10],100,1)',repmat(prain,10,1),repmat(prain,10,1))
%         %             colormap jet
%         colormap fire
%         shading interp
%         ax=gca;
%         set(ax,'XTickLabel','')
% 
%         %Generate a histogram for the magnitude of the responses
% 
%         %     handles_out.all_div_t=[];
%         %     handles_out.all_div_delta_dFFsplus=[];
%         %     handles_out.all_div_delta_dFFsminus=[];
% 
%         figureNo=figureNo+1;
%         try
%             close(figureNo)
%         catch
%         end
% 
%         hFig = figure(figureNo);
% 
%         set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
% 
%         edges=-0.5:0.5:4;
%         histogram(handles_out.all_div_delta_dFFsplus,edges,'FaceColor',[0 114/255 178/255])
%         hold on
%         histogram(handles_out.all_div_delta_dFFsminus,edges,'FaceColor',[158/255 31/255 99/255])
% 
%         ylabel('No of ROIs')
%         xlabel('Delta zdFF')
%         title('Histogram for Delta zdFF')
% 
%         this_ylim=ylim;
% 
%         text(2.5,this_ylim(1)+0.85*(this_ylim(2)-this_ylim(1)),'S+','Color',[0 114/255 178/255],'FontSize',16)
%         text(2.5,this_ylim(1)+0.75*(this_ylim(2)-this_ylim(1)),'S-','Color',[158/255 31/255 99/255],'FontSize',16)
% 
%         figureNo=figureNo+1;
%         try
%             close(figureNo)
%         catch
%         end
% 
%         hFig = figure(figureNo);
% 
%         set(hFig, 'units','normalized','position',[.2 .2 .3 .3])
% 
%         plot(handles_out.all_div_delta_dFFsminus,handles_out.all_div_delta_dFFsplus,'ok')
% 
%         ylabel('S+')
%         xlabel('S-')
%         title('Delta zdFF')
% 
%         %Divergence time
%         figureNo=figureNo+1;
%         try
%             close(figureNo)
%         catch
%         end
% 
%         hFig = figure(figureNo);
% 
%         set(hFig, 'units','normalized','position',[.2 .2 .4 .3])
% 
%         edges=[-1.75:0.5:5.25];
%         histogram(handles_out.all_div_t(~isnan(handles_out.all_div_t)),edges)
% 
%         ylabel('No of ROIs')
%         xlabel('Divegence time')
%         title('Histogram for divergence time')


    end
    %Uncomment this if you want to browse through the figures
%     tic
%     for figNo=1:figureNo
%         figure(figNo)
%         this_toc=toc;
%         while toc-this_toc<4
%         end
%     end
    pffft=1;

    save([choiceBatchPathName choiceFileName(1:end-4) '.mat'],'handles','handles_out','handles_choices','-v7.3')
end

pffft=1;