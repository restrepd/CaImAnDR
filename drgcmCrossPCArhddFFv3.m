%drgcmCrossPCArhddFFv2
clear all
close all

rng('shuffle')

use_raw=1; %If this is 1 the program uses the raw trace, otherwise it uses the inferred trace
 
seg_method=1; %If this variable is 0 the program uses random segmentation, otherwise it uses random circular shift

figNo=0;

k_fold=8;

%Theta
handles.lowF(1)=6;
handles.highF(1)=14;

%Beta
handles.lowF(2)=15;
handles.highF(2)=30;

%Low gamma
handles.lowF(3)=35;
handles.highF(3)=55;

%High gamma
handles.lowF(4)=65;
handles.highF(4)=95;

%swr
handles.lowF(4)=150;
handles.highF(4)=250;

%Names of bandwidths
handles.bw_names{1}='Theta';
handles.bw_names{2}='Beta';
handles.bw_names{3}='Low gamma';
handles.bw_names{4}='High gamma';
handles.bw_names{5}='swr';

%Ask user for the file with the traces
[choiceFileName,choiceBatchPathName] = uigetfile({'*batch_per_file.mat'},'Select the .mat file for analysis');
fprintf(1, ['\ndrgcmCrossPCArhddFF run for ' choiceFileName '\n\n']);



cd(choiceBatchPathName)
load(choiceFileName)

%rho
ii_rho=0;
rho_dFF=[];
pval_rho_dFF=[];
rho_ddFF=[];
pval_rho_ddFF=[];
rho_fileNo=[];
rho_dFFtrace_no=[];
rho_electrode_no=[];
rho_bwii=[];

%Shuffled rho
ii_rhos=0;
rhos_dFF=[];
pval_rhos_dFF=[];
rhos_ddFF=[];
pval_rhos_ddFF=[];
rhos_fileNo=[];
rhos_dFFtrace_no=[];
rhos_electrode_no=[];
rhos_bwii=[];

%rhos for glm
ii_rho_glm=0;
rho_glm=[];
rho_glm_bwii=[];
rho_glm_elctrode_no=[];

%Now let's read Connor's data
optical_data = readmatrix([choiceBatchPathName 'DFTemporal.csv']); %DF/F traces from CaImAn
highest_corr=zeros(4,16);
highest_corr_ROI=ones(4,16);
highest_rho_per_ROI=zeros(4,100);
elect_for_highest_rho_per_ROI=ones(4,100);
for fileNo=1:handles_per_file.no_files
    for bwii=1:4
        no_timepoints_dFF=size(handles_per_file.file(fileNo).dFFtraces,2);
        %         if use_raw==1
        %             dFFtraces=handles_per_file.file(fileNo).dFFtraces;
        %         else
        %             dFFtraces=handles_per_file.file(fileNo).dFFtraces_inferred;
        %         end
        dFFtraces=optical_data;
        no_traces_dFF=size(optical_data,1);
        
        LFPtime=handles_per_file.file(fileNo).log_P_time;
        dFFtime=handles_per_file.file(fileNo).dFF_time;
        dFFtime=dFFtime(1:size(dFFtraces,2));
        no_timepoints_dFF=size(dFFtraces,2);
        decimated_LFP_logPtraces=zeros(handles_per_file.file(fileNo).no_electrodes,no_timepoints_dFF);
        dt_dFF=dFFtime(2)-dFFtime(1);
        for elect_no=1:handles_per_file.file(fileNo).no_electrodes
            decimated_LFP_logPtraces(elect_no,1)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(elect_no,bwii,(LFPtime>dFFtime(1))&(LFPtime<=dFFtime(1)+(dt_dFF/2))),3);
            decimated_LFP_logPtraces(elect_no,end)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(elect_no,bwii,(LFPtime>dFFtime(end)-(dt_dFF/2))&(LFPtime<=dFFtime(end))),3);
            for ii_time=2:no_timepoints_dFF-1
                decimated_LFP_logPtraces(:,ii_time)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(:,bwii,(LFPtime>dFFtime(ii_time)-(dt_dFF/2))&(LFPtime<=dFFtime(ii_time)+(dt_dFF/2))),3);
            end
        end
        
        %Save the decimated data
        save([choiceFileName(1:end-19) '_dFF.txt'],'dFFtraces','-ascii')
        save([choiceFileName(1:end-19) '_LFP_' handles.bw_names{bwii} '.txt'],'decimated_LFP_logPtraces','-ascii')
        dFFtraces_t=dFFtraces';
        decimated_LFP_logPtraces_t=decimated_LFP_logPtraces';
        save([choiceFileName(1:end-19) '_dec_' handles.bw_names{bwii} '.mat'],'dFFtraces_t','decimated_LFP_logPtraces_t','dFFtime')
        
        %Calculate rho
        mask_dFF=~isnan(decimated_LFP_logPtraces(1,:));
        for elect_no=1:handles_per_file.file(fileNo).no_electrodes
            for trace_no=1:no_traces_dFF
                ii_rho=ii_rho+1;
                [rho_dFF(ii_rho),pval_rho_dFF(ii_rho)] = corr(dFFtraces(trace_no,mask_dFF)',decimated_LFP_logPtraces(elect_no,mask_dFF)');
                
                rho_fileNo(ii_rho)=fileNo;
                rho_dFFtrace_no(ii_rho)=trace_no;
                rho_electrode_no(ii_rho)=elect_no;
                rho_bwii(ii_rho)=bwii;
                
                if highest_corr(bwii,elect_no)<rho_dFF(ii_rho)
                    highest_corr(bwii,elect_no)=rho_dFF(ii_rho);
                    highest_corr_ROI(bwii,elect_no)=trace_no;
                end
                
                if highest_rho_per_ROI(bwii,trace_no)<rho_dFF(ii_rho)
                    highest_rho_per_ROI(bwii,trace_no)=rho_dFF(ii_rho);
                    elect_for_highest_rho_per_ROI(bwii,trace_no)=elect_no;
                end
                
            end
        end
        
        %Shuffle dFF
        no_shuffles=100;
        if seg_method==0
            %Now calculate rho for shuffled dFF by random segmentation
            no_segments=10;
            
            segments=[1:no_segments];
            per_seg=perms(segments);
            seg_ii=randi(size(per_seg,1),1,10000);
            ii_used=0;
            seg_length=floor(sum(mask_dFF)/no_segments);
            
            for shuf_ii=1:no_shuffles
                found_shuffle=0;
                while found_shuffle==0
                    ii_used=ii_used+1;
                    if sum(per_seg(seg_ii(ii_used),:)==segments)==0
                        found_shuffle=1;
                    end
                end
                
                for elect_no=1:handles_per_file.file(fileNo).no_electrodes
                    for trace_no=1:no_traces_dFF
                        shdFFtrace=zeros(1,no_segments*seg_length);
                        this_per_ii=per_seg(seg_ii(ii_used),:);
                        for ii=1:no_segments
                            shdFFtrace(1,(ii-1)*seg_length+1:ii*seg_length)=dFFtraces(trace_no,(this_per_ii(ii)-1)*seg_length+1:this_per_ii(ii)*seg_length);
                        end
                        ii_rhos=ii_rhos+1;
                        [rhos_dFF(ii_rhos),pval_rhos_dFF(ii_rhos)] = corr(shdFFtrace(1,:)',decimated_LFP_logPtraces(elect_no,1:length(shdFFtrace))');
                        
                        rhos_fileNo(ii_rhos)=fileNo;
                        rhos_dFFtrace_no(ii_rhos)=trace_no;
                        rhos_electrode_no(ii_rhos)=elect_no;
                        rhos_bwii(ii_rhos)=bwii;
                    end
                end
            end
            
        else
            %Now calculate rho for shuffled dFF using a circular shift
            
            
            shift_ii=ceil(0.1*size(dFFtraces,2)) + randi(ceil(0.9*size(dFFtraces,2)),1,1000);
            no_timepoints=size(dFFtraces,2);
            
            for shuf_ii=1:no_shuffles
                for elect_no=1:handles_per_file.file(fileNo).no_electrodes
                    for trace_no=1:no_traces_dFF
                        shdFFtrace=zeros(1,no_timepoints);
                        
                        %Shift forward
                        shdFFtrace(1,shift_ii(shuf_ii):end)=dFFtraces(trace_no,1:no_timepoints-shift_ii(shuf_ii)+1);
                        shdFFtrace(1,1:shift_ii(shuf_ii)-1)=dFFtraces(trace_no,no_timepoints-shift_ii(shuf_ii)+2:end);
                        
                        ii_rhos=ii_rhos+1;
                        [rhos_dFF(ii_rhos),pval_rhos_dFF(ii_rhos)] = corr(shdFFtrace(1,:)',decimated_LFP_logPtraces(elect_no,1:length(shdFFtrace))');
                        
                        rhos_fileNo(ii_rhos)=fileNo;
                        rhos_dFFtrace_no(ii_rhos)=trace_no;
                        rhos_electrode_no(ii_rhos)=elect_no;
                        rhos_bwii(ii_rhos)=bwii;
                    end
                end
            end
        end
        %          %Fit LFP with dFFtraces
        %         order=3;
        %         framelen=31;
        %         mask_dFF=~isnan(decimated_LFP_logPtraces(1,:));
        %         for elect_no=1:handles_per_file.file(fileNo).no_electrodes
        %             this_decimated_LFP_logPtraces=zeros(sum(mask_dFF),1);
        %             this_decimated_LFP_logPtraces(:,1)=decimated_LFP_logPtraces(elect_no,mask_dFF)';
        %             sgf_this_decimated_LFP_logPtraces = sgolayfilt(this_decimated_LFP_logPtraces,order,framelen);
        %             %Now do a k-fold glm fit where the data left out (1/8th) is predicted from
        %             %the data that was used for the glm (the other 7/8ths)
        %             k_fold_chunk=floor(length(sgf_this_decimated_LFP_logPtraces)/k_fold);
        %             LFPlog_pred=zeros(length(sgf_this_decimated_LFP_logPtraces),1);
        %             LFPlog_CI=zeros(length(sgf_this_decimated_LFP_logPtraces),2);
        %             for k_fold_ii=1:k_fold
        %                 mask_test_dFFtraces=zeros(1,length(sgf_this_decimated_LFP_logPtraces));
        %                 mask_test_dFFtraces((k_fold_ii-1)*k_fold_chunk+1:k_fold_ii*k_fold_chunk)=1;
        %                 mask_test_dFFtraces=logical(mask_test_dFFtraces);
        %                 mdl=fitglm(dFFtraces(:,~mask_test_dFFtraces)',sgf_this_decimated_LFP_logPtraces(~mask_test_dFFtraces),'linear');
        %                 [LFPlog_pred(mask_test_dFFtraces),LFPlog_CI(mask_test_dFFtraces,:)] = predict(mdl,dFFtraces(:,mask_test_dFFtraces)');
        %             end
        %             ii_rho_glm=ii_rho_glm+1;
        %             rho_glm_elctrode_no(ii_rho_glm)=elect_no;
        %             rho_glm_bwii(ii_rho_glm)=bwii;
        %             rho_glm(ii_rho_glm)=corr(LFPlog_pred,sgf_this_decimated_LFP_logPtraces);
        %
        %
        %             if figNo<5
        %                 figNo=figNo+1;
        %                 try
        %                     close(figNo)
        %                 catch
        %                 end
        %                 hFig=figure(figNo);
        %                 hold on
        %                 set(hFig, 'units','normalized','position',[.3 .3 .6 .3])
        %                 plot(dFFtime,sgf_this_decimated_LFP_logPtraces,'-k','LineWidth',3)
        %                 plot(dFFtime+3,LFPlog_pred,'-r')
        %                 title(['glm fit of logLFP of electrode No ' num2str(elect_no) ' for ' handles.bw_names{bwii}])
        %                 fprintf(1, ['\nProcessed glm for ' handles.bw_names{bwii} ' electrode No ' num2str(elect_no) '\n']);
        %             end
        %
        %             %Do shuffled glm
        %             LFPlog_pred=zeros(length(sgf_this_decimated_LFP_logPtraces),1);
        %             LFPlog_CI=zeros(length(sgf_this_decimated_LFP_logPtraces),2);
        %
        %             for k_fold_ii=1:k_fold
        %                 mask_test_dFFtraces=zeros(1,length(sgf_this_decimated_LFP_logPtraces));
        %                 mask_test_dFFtraces((k_fold_ii-1)*k_fold_chunk+1:k_fold_ii*k_fold_chunk)=1;
        %                 mask_test_dFFtraces=logical(mask_test_dFFtraces);
        %                 perm_k_fold_ii=randperm(k_fold);
        %                 sh_sgf_this_decimated_LFP_logPtraces=zeros(1,length(sgf_this_decimated_LFP_logPtraces));
        %                 for jj=1:k_fold
        %                     sh_sgf_this_decimated_LFP_logPtraces((jj-1)*k_fold_chunk+1:jj*k_fold_chunk)=...
        %                         sgf_this_decimated_LFP_logPtraces((perm_k_fold_ii(jj)-1)*k_fold_chunk+1:perm_k_fold_ii(jj)*k_fold_chunk);
        %                 end
        %                 mdl=fitglm(dFFtraces(:,~mask_test_dFFtraces)',sh_sgf_this_decimated_LFP_logPtraces(~mask_test_dFFtraces),'linear');
        %                 [LFPlog_pred(mask_test_dFFtraces),LFPlog_CI(mask_test_dFFtraces,:)] = predict(mdl,dFFtraces(:,mask_test_dFFtraces)');
        %             end
        %
        %             rho_glm_sh(ii_rho_glm)=corr(LFPlog_pred,sgf_this_decimated_LFP_logPtraces);
        %
        %         end
        
        %         %Calculate PCs per timepoint
        %         for ii_time=1:no_timepoints_dFF
        %             these_logPs=zeros(1,handles_per_file.file(fileNo).no_electrodes);
        %             these_logPs(1,:)=decimated_LFP_logPtraces(:,ii_time);
        %             this_PCs_logP=[];
        %             [coeff_logP,this_PCs_logP,latent_logP]=pca(these_logPs');
        %             if ii_time==1
        %                 PCs_logP=zeros(length(this_PCs_logP),no_timepoints_dFF);
        %             end
        %             PCs_logP(:,ii_time)=this_PCs_logP;
        %
        %             these_dFFs=zeros(1,no_traces_dFF);
        %             these_dFFs(1,:)=dFFtraces(:,ii_time);
        %             this_PCs_dFF=[];
        %             [coeff_logP,this_PCs_dFF,latent_logP]=pca(these_dFFs');
        %              if ii_time==1
        %                 PCs_dFF=zeros(length(this_PCs_dFF),no_timepoints_dFF);
        %             end
        %             PCs_dFF(:,ii_time)=this_PCs_dFF;
        %         end
        %
        %         figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        %
        %         hold on
        %
        %
        %
        %         PC1_dFF=PCs_dFF(1,:)';
        %         PC1_logP=PCs_logP(1,:)';
        %         PC_mask=(~isnan(PC1_dFF))&(~isnan(PC1_logP));
        %         PC1_dFF=PC1_dFF(PC_mask);
        %         PC1_logP=PC1_logP(PC_mask);
        %
        %         plot(PC1_dFF,PC1_logP,'.k')
        %         xlabel('PC1 for dFF')
        %         ylabel('PC1 for logP')
        %         title(['PC1 plot for ' handles.bw_names{bwii}])
        %
        %         [rho,pval] = corr(PC1_dFF,PC1_logP);
        %
        %         fprintf(1, ['\nFor file No %d ' handles.bw_names{bwii} 'PC1 dFF vs PC1 logP rho= %d, p value= %d\n'],fileNo,rho, pval);
        %
        %Now do dFF derivative
        %Convolve lick_freq using a window of 0.9 sec
        %         no_conv_points=5;
        %         conv_win=ones(1,no_conv_points);
        %         conv_dFFtraces=zeros(handles_per_file.file(fileNo).no_electrodes,no_timepoints_dFF);
        %
        %         for trace_no=1:no_traces_dFF
        %             conv_dFFtraces(trace_no,:)=conv(dFFtraces(trace_no,:),conv_win,'same')/no_conv_points;
        %         end
        %
        %         ddFFtraces=zeros(no_traces_dFF,size(dFFtraces,2));
        %         ddFFtraces(:,2:end)=(conv_dFFtraces(:,2:end)-conv_dFFtraces(:,1:end-1))/dt_dFF;
        %         ddFFtraces(:,1)=ddFFtraces(:,2);
        %         pffft=1;
        %               %Calculate PCs per timepoint
        %         for ii_time=1:no_timepoints_dFF
        %
        %
        %             these_ddFFtraces=zeros(1,no_traces_dFF);
        %             these_ddFFtraces(1,:)=ddFFtraces(:,ii_time);
        %             this_PCs_ddFF=[];
        %             [coeff_logP,this_PCs_ddFF,latent_logP]=pca(these_ddFFtraces');
        %              if ii_time==1
        %                 PCs_ddFF=zeros(length(this_PCs_ddFF),no_timepoints_dFF);
        %             end
        %             PCs_ddFF(:,ii_time)=this_PCs_ddFF;
        %         end
        %
        %         PCs_ddFF=PCs_ddFF(:,PC_mask);
        %         PC1_ddFF=PCs_ddFF(1,:);
        %
        %         figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        %
        %         hold on
        %
        %
        %         plot(PC1_ddFF,PC1_logP,'.k')
        %         xlabel('PC1 for derivative of dFF')
        %         ylabel('PC1 for logP')
        %         title(['PC1 plot for ' handles.bw_names{bwii}])
        %
        %         [rho,pval] = corr(PC1_dFF,PC1_logP);
        %
        %         fprintf(1, ['\nFor file No %d ' handles.bw_names{bwii} 'PC1 dFF derivative vs PC1 logP rho= %d, p value= %d\n'],fileNo,rho, pval);
        %
        %Now add all derivatives
        
        %         sumddFFtraces=sum(ddFFtraces,1);
        %          sumddFFtraces=sumddFFtraces(1,PC_mask);
        %          figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        %
        %         hold on
        %
        %
        %         plot(sumddFFtraces,PC1_logP,'.k')
        %         xlabel('Sum of derivatives fot dFF')
        %         ylabel('PC1 for logP')
        %         title(['PC1 plot for ' handles.bw_names{bwii}])
        %
        %         [rho,pval] = corr(sumddFFtraces',PC1_logP);
        %
        %         fprintf(1, ['\nFor file No %d ' handles.bw_names{bwii} 'sum of dFF derivatives vs PC1 logP rho= %d, p value= %d\n'],fileNo,rho, pval);
        %
    end
end
%
% %Plot correlation between glm prediction and LFP vs shuffled
% edges=[-1:0.05:1];
% rand_offset=0.8;
%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% set(hFig, 'units','normalized','position',[.3 .3 .6 .3])
% hold on
% bar_offset=0;
%
% glm_rho_glm=[];
% glm_ii=0;
%
% id_ii=0;
% input_data=[];
%
% for bwii=1:4
%
%     rho_glm_elctrode_no(ii_rho_glm)=elect_no;
%     rho_glm_bwii(ii_rho_glm)=bwii;
%     rho_glm(ii_rho_glm)=corr(LFPlog_pred,sgf_this_decimated_LFP_logPtraces);
%
%
%     %Shuffled rho
%     bar_offset=bar_offset+1
%     bar(bar_offset,mean(rho_glm_sh(rho_glm_bwii==bwii)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0.7 0.7 0.7])
%
%     %Violin plot
%     these_rhos=rho_glm_sh(rho_glm_bwii==bwii);
%     [mean_out, CIout]=drgViolinPoint(these_rhos...
%         ,edges,bar_offset,rand_offset,'k','k',3);
%
%     %Shuffled
%     glm_rho_glm.data(glm_ii+1:glm_ii+length(these_rhos))=these_rhos;
%     glm_rho_glm.bwii(glm_ii+1:glm_ii+length(these_rhos))=bwii*ones(1,length(these_rhos));
%     glm_rho_glm.shuffled(glm_ii+1:glm_ii+length(these_rhos))=ones(1,length(these_rhos));
%
%     glm_ii=glm_ii+length(these_rhos);
%
%     id_ii=id_ii+1;
%     input_data(id_ii).data=these_rhos;
%     input_data(id_ii).description=['rho for shuffled ' handles.bw_names{bwii}];
%
%
%     %rho
%     bar_offset=bar_offset+1
%     bar(bar_offset,mean(rho_glm(rho_glm_bwii==bwii)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 0.7 0.7])
%
%     %Violin plot
%     these_rhos=rho_glm(rho_glm_bwii==bwii);
%     [mean_out, CIout]=drgViolinPoint(these_rhos...
%         ,edges,bar_offset,rand_offset,'k','k',3);
%
%
%     bar_offset=bar_offset+1
%
%     %Original data
%     glm_rho_glm.data(glm_ii+1:glm_ii+length(these_rhos))=these_rhos;
%     glm_rho_glm.bwii(glm_ii+1:glm_ii+length(these_rhos))=bwii*ones(1,length(these_rhos));
%     glm_rho_glm.shuffled(glm_ii+1:glm_ii+length(these_rhos))=zeros(1,length(these_rhos));
%
%     glm_ii=glm_ii+length(these_rhos);
%
%     id_ii=id_ii+1;
%     input_data(id_ii).data=these_rhos;
%     input_data(id_ii).description=['rho for original ' handles.bw_names{bwii}];
%
%
% end
%
% title('Correlation between glm predict and LFP power')
% ylabel('rho')
%
% xticks([1.5 4.5 7.5 10.5])
% xticklabels({handles.bw_names{1}, handles.bw_names{2}, handles.bw_names{3}, handles.bw_names{4}})
%
% %Perform the glm
% fprintf(1, ['glm for rho between LFP log power and predicted using dFF\n'])
%
%
% tbl = table(glm_rho_glm.data',glm_rho_glm.bwii',glm_rho_glm.shuffled',...
%     'VariableNames',{'rho','bwii','shuffled'});
% mdl = fitglm(tbl,'rho~bwii+shuffled+bwii*shuffled'...
%     ,'CategoricalVars',[2,3])
%
%
% %Do the ranksum/t-test
% fprintf(1, ['\n\nRanksum or t-test p values for rho between LFP log power and predicted using dFF\n'])
%
% try
%     [output_data] = drgMutiRanksumorTtest(input_data);
% catch
% end

%Plot correlation histograms
edges=[-0.2:0.01:0.2];

pFDR=drsFDRpval(pval_rho_dFF);
fprintf(1,'\n\npFDR for rho p value  = %d\n',pFDR);

pFDRs=drsFDRpval(pval_rhos_dFF);
fprintf(1, ['\n\npFDR for shuffled rho p value = %d\n'],pFDR);

for bwii=1:4
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .6])
    
    subplot(2,1,1)
    hold on
    histogram(rho_dFF(rho_bwii==bwii),edges);

    title('Original')
    xlabel('rho')
    
    
    
    %Plot shuffled correlation histograms
    
    subplot(2,1,2)
    

    hold on
    histogram(rhos_dFF(rhos_bwii==bwii),edges,'FaceColor',[0.3010 0.7450 0.9330]);

    title('Shuffled')
    xlabel('rho')
    
    sgtitle(['rho for dFF x LFP log P for ' handles.bw_names{bwii} ])
    
    [h,p] = kstest2(rhos_dFF(rhos_bwii==bwii),rho_dFF(rho_bwii==bwii),'Alpha',0.01);
    fprintf(1, ['\np value for thr KS test difference between original and shuffled ' handles.bw_names{bwii}  ' = %d\n\n'],p);
    
    
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    hold on
    
    [f_rho,x_rho] = drg_ecdf(rho_dFF);
    plot(x_rho,f_rho,'b','LineWidth',3)
    
    [f_rhos,x_rhos] = drg_ecdf(rhos_dFF);
    plot(x_rhos,f_rhos,'Color',[0.3010 0.7450 0.9330],'LineWidth',3)
    
    text(0.2,0.6,'Original','Color','b','FontSize',12)
    text(0.2,0.5,'Shuffled','Color',[0.3010 0.7450 0.9330],'FontSize',12)

    title(['Cumulative probability for rho ' handles.bw_names{bwii}])
    xlabel('rho')
    ylabel('Cumulative probability')
    
end


figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.3 .3 .4 .4])

hold on

ii=0;
for bwii=1:4
    ii=ii+1;
    bar(ii,100*sum(pval_rho_dFF(rho_bwii==bwii)<pFDR)/sum(rho_bwii==bwii),'r')
    ii=ii+1;
    bar(ii,100*sum(pval_rhos_dFF(rhos_bwii==bwii)<pFDRs)/sum(rhos_bwii==bwii),'b')
    ii=ii+1;
end

xticks([1.5 4.5 7.5 10.5])
xticklabels({'Theta','Beta','Low gamma','High gamma'})
ylim([0 50])
ylabel('Percent significant rho')
title('Percent significant correlations')
text(10,38,'Original','Color','red','FontSize',12)
text(10,36,'Shuffled','Color','blue','FontSize',12)


figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.3 .3 .3 .3])

hold on

ii=0;
for bwii=1:4
    bar(ii,100*sum(pval_rho_dFF(rho_bwii==bwii)<pFDR)/sum(rho_bwii==bwii)-100*sum(pval_rhos_dFF(rhos_bwii==bwii)<pFDRs)/sum(rhos_bwii==bwii),'b')
    ii=ii+2;
end

xticks([0 1 2 3])
xticklabels({'Theta','Beta','Low gamma','High gamma'})
ylim([0 25])
ylabel('Percent significant rho')
title('Percent significant correlations')


writematrix(highest_corr,[choiceBatchPathName 'highest_corr.csv'])
writematrix(highest_corr_ROI,[choiceBatchPathName 'highest_corr_ROI.csv'])
writematrix(highest_rho_per_ROI',[choiceBatchPathName 'highest_rho_per_ROI.csv'])
writematrix(elect_for_highest_rho_per_ROI',[choiceBatchPathName 'elect_for_highest_rho_per_ROI.csv'])


pffft=1;