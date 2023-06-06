%summary_decoding_training_window
close all
clear all
figureNo=0;

read_window_label{1}='Pre';
read_window_label{2}='Odor';
read_window_label{3}='Reinf';

%Time windows
%     pre_t=[-1 0];
%     odor_t=[3.1 4.1];
%     reinf_t=[4.5 5.5];

pre_t=[-2.5 -1.5];
odor_t=[-1 0];
reinf_t=[2 4.1];

%This was used previously
%     time_periods_eu=[-5 -3;
%         -2.5 -1.5;
%         -1 0;
%         2 4.1];

time_periods_eu=[
    -1 0;
    3.1 4.1;
    4.4 5.4];

edges=[0:0.05:1];
rand_offset=0.8;


        per_corr_set_label{1}='naive';
        per_corr_set_label{2}='intermediate';
        per_corr_set_label{3}='proficient';


handles_pf.pre_per_outPathName='F:\SFTP\Ming Ma\Grin678and9_spm_noneg\';

fileID = fopen([handles_pf.pre_per_outPathName 'summary_decoding_training_windows.txt'],'w');

%They are ordered in training grant widow order, with the exception that
%the last one is the long 0.5 to 5.5 window
handles_pf.pre_per_FileName{6}='drgCaImAn_LDAfsdz_choices_spm_noneg_Grin678and9_05212023.mat';
handles_pf.pre_per_label{6}='0.5 to 5.5';

handles_pf.pre_per_FileName{1}='drgCaImAn_LDAfsdz_choices_spm_noneg_Grin678and9_05232023.mat';
handles_pf.pre_per_label{1}='-5 to -3';

handles_pf.pre_per_FileName{4}='drgCaImAn_LDAfsdz_choices_spm_noneg_Grin678and9_05252023.mat';
handles_pf.pre_per_label{4}='2 to 4';

handles_pf.pre_per_FileName{2}='drgCaImAn_LDAfsdz_choices_spm_noneg_Grin678and9_05282023.mat';
handles_pf.pre_per_label{2}='-2 to 0';

handles_pf.pre_per_FileName{5}='drgCaImAn_LDAfsdz_choices_spm_noneg_Grin678and9_05292023.mat';
handles_pf.pre_per_label{5}='4.5 to 6.5';

handles_pf.pre_per_FileName{3}='drgCaImAn_LDAfsdz_choices_spm_noneg_Grin678and9_05312023.mat';
handles_pf.pre_per_label{3}='0 to 2';

all_session_accs_pf=[];

for fileNo=1:length(handles_pf.pre_per_FileName)
    load([handles_pf.pre_per_outPathName handles_pf.pre_per_FileName{fileNo}])
    all_session_accs_pf(fileNo).all_accs=handles_out2.all_accs;
end

%Plot separate bar graphs for each read window
for read_window=2:3
    figureNo = figureNo + 1;
    try
        close(figureNo)
    catch
    end
    hFig=figure(figureNo);
    hold on

    ax=gca;ax.LineWidth=3;
    set(hFig, 'units','normalized','position',[.2 .2 .6 .3])

    glm_acc=[];
    glm_acc_ii=0;

    glm_acc_no_sh=[];
    glm_acc_no_sh_ii=0;


    id_acc_ii=0;
    input_acc_data=[];

    % for grNo=1:no_pcorr*length(these_groups)
    bar_offset=0;


    %
    %     window_label{1}='Base';
    %     window_label{2}='PreFV';
    %     window_label{3}='PreOdor';
    %     window_label{4}='Odor';





    %     per_corr_set_label{1}='naive';
    %     per_corr_set_label{2}='intermediate';
    %     per_corr_set_label{3}='proficient';

    
    for perCorr_no=[1 3]
        for train_window_no=1:length(handles_pf.pre_per_FileName)

            per_mouse_mean_accuracy=[];
            per_mouse_mean_accuracy_sh=[];
            all_session_accuracy=[];
            all_session_accuracy_sh=[];

             

            %Get per session and per mouse
            ii_m_included=0;
            for mouseNo=1:length(all_session_accs_pf(train_window_no).all_accs.perCorr(perCorr_no).decode_window(read_window).mouse)
               

                include_mouse=all_session_accs_pf(train_window_no).all_accs.perCorr(perCorr_no).decode_window(read_window).mouse(mouseNo).mouse_included;
%                 for grNo=group_sets(perCorr_no,:)
%                     if ~isempty(per_mouse_acc.group(grNo).mouse(mouseNo).mean_accuracy)
%                         include_mouse=1;
%                         for ii_sessions=1:size(per_mouse_acc.group(grNo).mouse(mouseNo).accuracy)
%                             t_from=time_periods_eu(train_window_no,1);
%                             t_to=time_periods_eu(train_window_no,2);
%                             this_ac=mean(per_mouse_acc.group(grNo).mouse(mouseNo).accuracy(ii_sessions,(time_span>=t_from)&(time_span<=t_to)) );
%                             this_ac_sh=mean(per_mouse_acc.group(grNo).mouse(mouseNo).accuracy_sh(ii_sessions,(time_span>=t_from)&(time_span<=t_to)) );
%                             these_accs=[these_accs this_ac];
%                             these_accs_sh=[these_accs_sh this_ac_sh];
%                         end
%                     end
%                 end
                if include_mouse==1
                    ii_m_included=ii_m_included+1;
                     these_accs=[];
                    these_accs_sh=[];
                    these_accs=all_session_accs_pf(train_window_no).all_accs.perCorr(perCorr_no).decode_window(read_window).mouse(mouseNo).per_mouse_accuracy;
                    these_accs_sh=all_session_accs_pf(train_window_no).all_accs.perCorr(perCorr_no).decode_window(read_window).mouse(mouseNo).per_mouse_accuracy_sh;
                    per_mouse_mean_accuracy(1,ii_m_included)=mean(these_accs);
                    per_mouse_mean_accuracy_sh(1,ii_m_included)=mean(these_accs_sh);
                    all_session_accuracy=[all_session_accuracy these_accs];
                    all_session_accuracy_sh=[all_session_accuracy_sh these_accs_sh];

                    glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=these_accs_sh;
                    glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=perCorr_no*ones(1,length(these_accs_sh));
                    glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=train_window_no*ones(1,length(these_accs_sh));
                    glm_acc.shuffled(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=ones(1,length(these_accs_sh));
                    glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs_sh))=mouseNo*ones(1,length(these_accs_sh));
                    glm_acc_ii=glm_acc_ii+length(these_accs_sh);

                    glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
                    glm_acc.pcorr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=perCorr_no*ones(1,length(these_accs));
                    glm_acc.window(glm_acc_ii+1:glm_acc_ii+length(these_accs))=train_window_no*ones(1,length(these_accs));
                    glm_acc.shuffled(glm_acc_ii+1:glm_acc_ii+length(these_accs))=zeros(1,length(these_accs));
                    glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=mouseNo*ones(1,length(these_accs));
                    glm_acc_ii=glm_acc_ii+length(these_accs);

                    glm_acc_no_sh.data(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=these_accs;
                    glm_acc_no_sh.pcorr(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=perCorr_no*ones(1,length(these_accs));
                    glm_acc_no_sh.window(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=train_window_no*ones(1,length(these_accs));
                    glm_acc_no_sh.shuffled(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=ones(1,length(these_accs));
                    glm_acc_no_sh.mouse_nos(glm_acc_no_sh_ii+1:glm_acc_no_sh_ii+length(these_accs))=mouseNo*ones(1,length(these_accs));
                    glm_acc_no_sh_ii=glm_acc_no_sh_ii+length(these_accs);
                end

            end

            id_acc_ii=id_acc_ii+1;
            input_acc_data(id_acc_ii).data=all_session_accuracy_sh;
            input_acc_data(id_acc_ii).description=['Shuffled ' per_corr_set_label{perCorr_no} ' ' handles_pf.pre_per_label{train_window_no}];


            id_acc_ii=id_acc_ii+1;
            input_acc_data(id_acc_ii).data=these_accs;
            input_acc_data(id_acc_ii).description=[per_corr_set_label{perCorr_no} ' ' handles_pf.pre_per_label{train_window_no}];

            handles_pf_out2.all_accs.perCorr(perCorr_no).decode_window(train_window_no).all_session_accuracy_sh=all_session_accuracy_sh;
            handles_pf_out2.all_accs.perCorr(perCorr_no).decode_window(train_window_no).all_session_accuracy=all_session_accuracy;

            %Accuracy shuffled
            bar(bar_offset,mean(all_session_accuracy_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[200/255 200/255 200/255])

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(all_session_accuracy_sh...
                ,edges,bar_offset,rand_offset,'k','k',2);



            bar_offset=bar_offset+1;

            %Accuracy
            bar(bar_offset,mean(all_session_accuracy),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])

            %Violin plot
            [mean_out, CIout]=drgViolinPoint(all_session_accuracy...
                ,edges,bar_offset,rand_offset,'k','k',2);



            for ii_mouse=1:ii_m_included
                plot([bar_offset-1 bar_offset],[per_mouse_mean_accuracy_sh(ii_mouse) per_mouse_mean_accuracy(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
            end

            bar_offset=bar_offset+1;


        end
        bar_offset=bar_offset+1;
    end

    %     xticks([1 3 5 7 9 11 13 15])
    %     xticklabels({'Base','PreFV','PreOdor','Odor','Base','PreFV','PreOdor','Odor'})
    xticks([0.5 2.5 4.5 6.5 8.5 10.5 14.5 16.5 18.5 20.5 22.5 24.5])
    xticklabels({'-5 to -3','-2 to 0','0 to 2','2 to 4','4.5 to 6.5','0.5 to 5.5','-5 to -3','-2 to 0','0 to 2','2 to 4','4.5 to 6.5','0.5 to 5.5'})


    title(['Accuracy for different training windows for ' read_window_label{read_window}])
    ylabel('Accuracy')
    ylim([0.4 1.1])
    xlim([-1 26.5])

    %     text(1,1,'45-65%')
    %     text(6,1,'65-80%')
    %     text(11,1,'>80%')

    %Perform the glm including shuffled
    fprintf(1, ['\nglm for decoding accuracy including shuffled\n'])
    fprintf(fileID, ['\nglm for decoding accuracy including shuffled\n']);

    fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
    fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

    tbl = table(glm_acc.data',glm_acc.pcorr',glm_acc.window',glm_acc.shuffled',...
        'VariableNames',{'accuracy','percent_correct','window','shuffled'});
    mdl = fitglm(tbl,'accuracy~percent_correct+window+shuffled+percent_correct*window*shuffled'...
        ,'CategoricalVars',[2,3,4])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);

    %Perform the glm not including shuffled
    fprintf(1, ['\nglm for decoding accuracy not including shuffled\n'])
    fprintf(fileID, ['\nglm for decoding accuracy not including shuffled\n']);

    fprintf(1, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n'])
    fprintf(fileID, ['\nwindow 0: Base, 1:PreFV, 2:PreOdor, 3:Odor\n']);

    tbl = table(glm_acc_no_sh.data',glm_acc_no_sh.pcorr',glm_acc_no_sh.window',...
        'VariableNames',{'accuracy','percent_correct','window'});
    mdl = fitglm(tbl,'accuracy~percent_correct+window+percent_correct*window'...
        ,'CategoricalVars',[2,3])


    txt = evalc('mdl');
    txt=regexp(txt,'<strong>','split');
    txt=cell2mat(txt);
    txt=regexp(txt,'</strong>','split');
    txt=cell2mat(txt);

    fprintf(fileID,'%s\n', txt);


    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for decoding accuracy\n'])
    fprintf(fileID, ['\n\nRanksum or t-test p values for decoding accuracy\n']);


    [output_data] = drgMutiRanksumorTtest(input_acc_data, fileID,0);

end

fclose(fileID);