close all
clear all

no_mice=3;
 
%mmG7f09 processed
outFileName='20180608_mmG7f09_Cerebellum_lda.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180608_mmG7f09_Cerebellum/';
load([outPathName outFileName])
handles_par_per_mouse(1).handles_par=handles_par;
handles_par_per_mouse(1).handles_sig=handles_sig;

%mmPVG04 processed
outFileName='20180917and19_mmPVG04_Cerebellum_lda.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum/';
load([outPathName outFileName])
handles_par_per_mouse(2).handles_par=handles_par;
handles_par_per_mouse(2).handles_sig=handles_sig;

%mmG06 processed
outFileName='20180419and23_mmG06_cerebellum_lda.mat';
outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum/';
load([outPathName outFileName])
handles_par_per_mouse(3).handles_par=handles_par;
handles_par_per_mouse(3).handles_sig=handles_sig;

edges=0:2.5:100;
rand_offset=0.8;

figNo=0;

for winNo=3:-1:2
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
        bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct),'FaceColor',[1 0 0])
        plot(x_val*ones(1,length(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct)),handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
        plot(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct),'ok','MarkerFaceColor','k')
    end
    
    x_val=x_val+1;
    
     %Percent correct before reversal (shufffled)
     sample_no=2;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct),'FaceColor',[0 0 1])
        drgViolinPoint(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct,edges,x_val,rand_offset,'k','k',2);
    end
    
       x_val=x_val+3;
    
    %Percent correct after reversal
    sample_no=3;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct),'FaceColor',[1 0 0])
        plot(x_val*ones(1,length(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct)),handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2)
        plot(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct),'ok','MarkerFaceColor','k')
    end
    
    x_val=x_val+1;
     
    %Percent correct after reversal (shuffled)
     sample_no=4;
    for mouseNo=1:no_mice
        x_val=x_val+1;
        bar(x_val,mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct),'FaceColor',[0 0 1])
        drgViolinPoint(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct,edges,x_val,rand_offset,'k','k',2);
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
end
 
samples_to_compare=[1 2 3 4]
sample_description{1}='Original forward';
sample_description{2}='Shuffled forward';
sample_description{3}='Original reversed';
sample_description{4}='Shuffled reversed';

for winNo=3:-1:2
    p_vals_LDA=0;
    no_comps=0;
    fprintf(1, ['\n\np values for dFF for window %d\n\n'],winNo);
      glm_LDA=[];
    glm_ii=0;
    LDAs=[];
    LDAs_ii=0;
    for sampleNo1=1:4
        these_LDAs=[];
        for mouseNo=1:3
            sample_no=samples_to_compare(sampleNo1);
            these_LDAs(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct);
        end
        
        glm_LDA.data(glm_ii+1:glm_ii+3)=these_LDAs;
        LDAs_ii=LDAs_ii+1;
        LDAs(LDAs_ii).data=these_LDAs;
        switch sampleNo1
            case 1
                glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=1;
                glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=1;
                LDAs(LDAs_ii).description='Odor 1 Forward';
            case 2
                glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=2;
                glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=1;
                LDAs(LDAs_ii).description='Odor 2 Forward';
            case 3
                glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=1;
                glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=2;
                LDAs(LDAs_ii).description='Odor 1 Reverse';
            case 4
                glm_LDA.orig_shuffled(glm_ii+1:glm_ii+3)=2;
                glm_LDA.fwd_rev(glm_ii+1:glm_ii+3)=2;
                LDAs(LDAs_ii).description='Odor 2 Reverse';
        end
        glm_ii=glm_ii+3;
        for sampleNo2=sampleNo1+1:4
            no_comps=no_comps+1;
            LDA1=[];
            LDA2=[];
            for mouseNo=1:3
                sample_no=samples_to_compare(sampleNo1);
                LDA1(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct);
                sample_no=samples_to_compare(sampleNo2);
                LDA2(mouseNo)= mean(handles_par_per_mouse(mouseNo).handles_sig.win(winNo).disc_corr(sample_no).discriminant_correct);
            end
            [h p_vals_LDA(no_comps)]=ttest2(LDA1,LDA2);
            fprintf(1, ['p values ranksum for window %d ' sample_description{sampleNo1} ' vs. ' sample_description{sampleNo2} ' =%d\n'],winNo,p_vals_LDA(no_comps));
        end
    end
    pFDRLDA=drsFDRpval(p_vals_LDA);
    fprintf(1, ['pFDR for window %d  = %d\n\n'],winNo, pFDRLDA);
    
     fprintf(1, ['\n\nglm for LDA for window' num2str(winNo) '\n'])
    tbl = table(glm_LDA.data',glm_LDA.orig_shuffled',glm_LDA.fwd_rev',...
        'VariableNames',{'per_corr_LDA','orig_shuffled','fwd_rev'});
    mdl = fitglm(tbl,'per_corr_LDA~orig_shuffled+fwd_rev+orig_shuffled*fwd_rev'...
        ,'CategoricalVars',[2,3])
    
    fprintf(1, ['\n\nRanksum or t-test p values for LDA for window ' num2str(winNo) '\n'])
    try
        [output_data] = drgMutiRanksumorTtest(dFFs);
        fprintf(1, '\n\n')
    catch
    end
end

pffft=1;
