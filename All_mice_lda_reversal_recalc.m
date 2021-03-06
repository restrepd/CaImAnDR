close all
clear all

no_mice=3;
 
%mmG7f09 processed
outFileName='20180608_mmG7f09_Cerebellum_lda.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180608_mmG7f09_Cerebellum new analysis/';
load([outPathName outFileName])
handles_par_per_mouse(1).handles_par=handles_par;
handles_par_per_mouse(1).handles_sig=handles_sig;

%mmPVG04 processed
outFileName='20180917and19_mmPVG04_Cerebellum_lda.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/';
load([outPathName outFileName])
handles_par_per_mouse(2).handles_par=handles_par;
handles_par_per_mouse(2).handles_sig=handles_sig;

%mmG06 processed
outFileName='20180419and23_mmG06_cerebellum_lda.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/';
load([outPathName outFileName])
handles_par_per_mouse(3).handles_par=handles_par;
handles_par_per_mouse(3).handles_sig=handles_sig;

edges=0:2.5:100;
rand_offset=0.8;

figNo=0;

winNo=1;
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.25 .2+0.1*(winNo-1) .5 .5])
hold on

x_val=0;

%Percent correct before reversal
sample_no=1;
for mouseNo=1:no_mice
    x_val=x_val+1;
    bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(1).discriminant_correct),'FaceColor',[1 0 0])
    plot(x_val*ones(1,length(handles_par_per_mouse(mouseNo).handles_sig.win(1).discriminant_correct)),handles_par_per_mouse(mouseNo).handles_sig.win(1).discriminant_correct,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    plot(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(1).discriminant_correct),'ok','MarkerFaceColor','k')
end

x_val=x_val+1;

%Percent correct before reversal (shufffled)
sample_no=2;
for mouseNo=1:no_mice
    x_val=x_val+1;
    bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(3).discriminant_correct),'FaceColor',[0 0 1])
    drgViolinPoint(handles_par_per_mouse(mouseNo).handles_sig.win(3).discriminant_correct,edges,x_val,rand_offset,'k','k',2);
end

x_val=x_val+3;

%Percent correct after reversal
sample_no=3;
for mouseNo=1:no_mice
    x_val=x_val+1;
    bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(4).discriminant_correct),'FaceColor',[1 0 0])
    plot(x_val*ones(1,length(handles_par_per_mouse(mouseNo).handles_sig.win(4).discriminant_correct)),handles_par_per_mouse(mouseNo).handles_sig.win(4).discriminant_correct,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
    plot(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(4).discriminant_correct),'ok','MarkerFaceColor','k')
end

x_val=x_val+1;

%Percent correct after reversal (shuffled)
sample_no=4;
for mouseNo=1:no_mice
    x_val=x_val+1;
    bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(6).discriminant_correct),'FaceColor',[0 0 1])
    drgViolinPoint(handles_par_per_mouse(mouseNo).handles_sig.win(6).discriminant_correct,edges,x_val,rand_offset,'k','k',2);
end

xlim([0 18])
ylim([30 110])
title(['Percent correct in decoding analysis for window No' num2str(winNo)])
text(3,107.5,'Forward','FontSize',18)
plot([0.5 7.5],[105 105],'-k','LineWidth',3)
plot([10.5 17.5],[105 105],'-k','LineWidth',3)
text(13,107.5,'Reversed','FontSize',18)
ylabel('Percent correct')
text(0.5,102,'Original','FontSize',18,'Color','r')
text(4.5,102,'Shuffled','FontSize',18,'Color','b')
text(10.5,102,'Original','FontSize',18,'Color','b')
text(14.5,102,'Shuffled','FontSize',18,'Color','r')

samples_to_compare=[1 2 3 4]
sample_description{1}='Original forward';
sample_description{2}='Shuffled forward';
sample_description{3}='Original reversed';
sample_description{4}='Shuffled reversed';

winNo=1
p_vals_LDA=0;
no_comps=0;
fprintf(1, ['\n\np values for dFF for window %d\n\n'],winNo);
glm_LDA=[];
glm_ii=0;
LDAs=[];
LDAs_ii=0;

these_LDAs=[];
for mouseNo=1:3
    these_LDAs(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(1).discriminant_correct);
end
glm_LDA.data(glm_ii+1:glm_ii+3)=these_LDAs;
glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=1;
glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=1;
glm_ii=glm_ii+3;

LDAs_ii=LDAs_ii+1;
LDAs(LDAs_ii).data=these_LDAs;
LDAs(LDAs_ii).description='Forward';


these_LDAs=[];
for mouseNo=1:3
    these_LDAs(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(3).discriminant_correct);
end
glm_LDA.data(glm_ii+1:glm_ii+3)=these_LDAs;
glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=2;
glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=1;
glm_ii=glm_ii+3;

LDAs_ii=LDAs_ii+1;
LDAs(LDAs_ii).data=these_LDAs;
LDAs(LDAs_ii).description='Forward shuffled';

these_LDAs=[];
for mouseNo=1:3
    these_LDAs(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(4).discriminant_correct);
end
glm_LDA.data(glm_ii+1:glm_ii+3)=these_LDAs;
glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=1;
glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=2;
glm_ii=glm_ii+3;

LDAs_ii=LDAs_ii+1;
LDAs(LDAs_ii).data=these_LDAs;
LDAs(LDAs_ii).description='Reverse';

these_LDAs=[];
for mouseNo=1:3
    these_LDAs(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(6).discriminant_correct);
end
glm_LDA.data(glm_ii+1:glm_ii+3)=these_LDAs;
glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=2;
glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=2;
glm_ii=glm_ii+1;

LDAs_ii=LDAs_ii+1;
LDAs(LDAs_ii).data=these_LDAs;
LDAs(LDAs_ii).description='Reverse shuffled';

fprintf(1, ['\n\nglm for LDA for window' num2str(winNo) '\n'])
tbl = table(glm_LDA.data',glm_LDA.orig_shuffled',glm_LDA.fwd_rev',...
    'VariableNames',{'per_corr_LDA','orig_shuffled','fwd_rev'});
mdl = fitglm(tbl,'per_corr_LDA~orig_shuffled+fwd_rev+orig_shuffled*fwd_rev'...
    ,'CategoricalVars',[2,3])

fprintf(1, ['\n\nRanksum or t-test p values for LDA for window ' num2str(winNo) '\n'])
try
    [output_data] = drgMutiRanksumorTtest(LDAs);
    fprintf(1, '\n\n')
catch
end


pffft=1;
