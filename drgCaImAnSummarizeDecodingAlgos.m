%drgCaImAnSummarizeDecodingAlgos
%
%This program takes the output from
%drgCaImAn_analyze_batch_pre_per_to_decode_per_mouse_v2.m
%and generates a bar graph comaring decoding algorithms in the odor window


% time_periods_eu=[
%             -1 0;
%             3.1 4.1;
%             4.4 5.4];

clear all
close all

choiceBatchPathName='/Users/restrepd/Documents/Projects/SFTP/Fabio/algos/';
fileID = fopen([choiceBatchPathName 'drgCaImAnSummarizeDecodingAlgos.txt'],'w');

ii_pcorr=3; %proficient
ii_window=2; %odor window

%Load the files
these_groups_out=1:5;
group_legends{1}='GLM';
group_legends{2}='SVM';
group_legends{3}='BDT';
group_legends{4}='ANN';
group_legends{5}='LDA';

figureNo=0;


all_handles=[];

load('/Users/restrepd/Documents/Projects/SFTP/Fabio/algos/drgCaImAn_LDAfsdz_choices_spm_noneg_algo_Grin678and9_GLMDR.mat')
all_handles(1).handles_out2=handles_out2;

load('/Users/restrepd/Documents/Projects/SFTP/Fabio/algos/drgCaImAn_LDAfsdz_choices_spm_noneg_algo_Grin678and9_SVDRM.mat')
all_handles(2).handles_out2=handles_out2;

load('/Users/restrepd/Documents/Projects/SFTP/Fabio/algos/drgCaImAn_LDAfsdz_choices_spm_noneg_algo_Grin678and9_DTDR.mat')
all_handles(3).handles_out2=handles_out2;

load('/Users/restrepd/Documents/Projects/SFTP/Fabio/algos/drgCaImAn_LDAfsdz_choices_spm_noneg_algo_Grin678and9_ANRDR')
all_handles(4).handles_out2=handles_out2;

load('/Users/restrepd/Documents/Projects/SFTP/Fabio/algos/drgCaImAn_LDAfsdz_choices_spm_noneg_algo_Grin678and9_LDADR.mat')
all_handles(5).handles_out2=handles_out2;






%Plot a bar graph with of accuracy for each agorithm and shuffled
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

id_acc_ii=0;
input_acc_data=[];

% for grNo=1:no_pcorr*length(these_groups)
bar_offset=0;




edges=[0:0.05:1];
rand_offset=0.7;



for grNo=these_groups_out

        %For shuffled
        these_accs=[];
        these_mean_per_mouse_sh_accs=[];
        these_mice=[];
        for ii_mouse=1:4
           these_accs=[these_accs  all_handles(grNo).handles_out2.all_accs.perCorr(ii_pcorr).decode_window(ii_window).mouse(ii_mouse).per_mouse_accuracy_sh];
           these_mice=[these_mice ii_mouse*ones(1,length(all_handles(grNo).handles_out2.all_accs.perCorr(ii_pcorr).decode_window(ii_window).mouse(ii_mouse).per_mouse_accuracy_sh))];
           these_mean_per_mouse_sh_accs=[these_mean_per_mouse_sh_accs  mean(all_handles(grNo).handles_out2.all_accs.perCorr(ii_pcorr).decode_window(ii_window).mouse(ii_mouse).per_mouse_accuracy_sh)];
        end

        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[200/255 200/255 200/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_accs...
            ,edges,bar_offset,rand_offset,'k','k',3);

        

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        glm_acc.shift(glm_acc_ii+1:glm_acc_ii+length(these_accs))=ones(1,length(these_accs));
        glm_acc.algo(glm_acc_ii+1:glm_acc_ii+length(these_accs))=grNo*ones(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_mice;
        glm_acc_ii=glm_acc_ii+length(these_accs);

        id_acc_ii=id_acc_ii+1;
        input_acc_data(id_acc_ii).data=these_accs;
        input_acc_data(id_acc_ii).description=['Shuffled ' group_legends{grNo}];

        bar_offset=bar_offset+1;

        %Odor window
        
        these_accs=[];
        these_mean_per_mouse_accs=[];
        these_mice=[];
        for ii_mouse=1:4
           these_accs=[these_accs  all_handles(grNo).handles_out2.all_accs.perCorr(ii_pcorr).decode_window(ii_window).mouse(ii_mouse).per_mouse_accuracy];
           these_mice=[these_mice ii_mouse*ones(1,length(all_handles(grNo).handles_out2.all_accs.perCorr(ii_pcorr).decode_window(ii_window).mouse(ii_mouse).per_mouse_accuracy))];
           these_mean_per_mouse_accs=[these_mean_per_mouse_accs  mean(all_handles(grNo).handles_out2.all_accs.perCorr(ii_pcorr).decode_window(ii_window).mouse(ii_mouse).per_mouse_accuracy)];
        end

        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[120/255 120/255 120/255])

        %Violin plot
        [mean_out, CIout]=drgViolinPoint(these_accs...
            ,edges,bar_offset,rand_offset,'k','k',3);

          for ii_mouse=1:length( these_mean_per_mouse_accs)
                plot([bar_offset-1 bar_offset],[ these_mean_per_mouse_sh_accs(ii_mouse)  these_mean_per_mouse_accs(ii_mouse)],'-','Color',[150/255 150/255 150/255],'LineWidth',2 )
            end

        glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
        glm_acc.shift(glm_acc_ii+1:glm_acc_ii+length(these_accs))=zeros(1,length(these_accs));
        glm_acc.algo(glm_acc_ii+1:glm_acc_ii+length(these_accs))=grNo*ones(1,length(these_accs));
        glm_acc.mouse_nos(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_mice;
        glm_acc_ii=glm_acc_ii+length(these_accs);

        id_acc_ii=id_acc_ii+1;
        input_acc_data(id_acc_ii).data=these_accs;
        input_acc_data(id_acc_ii).description=[group_legends{grNo}];

        bar_offset=bar_offset+2;

  

end

xticks([0.5 3.5 6.5 9.5 12.5])
xticklabels({'GLM','SVM','BDT','ANN','LDA'})
% xtickangle(45)
title(['Decoding accuracy'])
ylabel('Accuracy')
ylim([0.4 1.1])
xlim([-1 14])


%Perform the glm
fprintf(1, ['\nglm for accuracy for different decoding algorithms\n'])
fprintf(fileID, ['\nglm for accuracy for different decoding algorithms\n']);

tbl = table(glm_acc.data',glm_acc.shift',glm_acc.algo',glm_acc.mouse_nos',...
    'VariableNames',{'accuracy','shift','algo','mouse'});
mdl = fitglm(tbl,'accuracy~shift+algo+mouse+shift*algo'...
    ,'CategoricalVars',[2,3,4])


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

fclose(fileID)


