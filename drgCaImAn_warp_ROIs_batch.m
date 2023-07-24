function drgCaImAn_warp_ROIs_batch
%% drgCaImAn_warp_ROIs_batch.m
%
% This code aligns the ROIs for different runs within a single behavioral
% session to allow studies in the go-no go reversal
%
close all
clear all

overlap_threshold=0.3;

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_rev_choices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_LDA_fsdz run for ' choiceFileName '\n\n']);

tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

figNo=0;

%Reference file
this_pre_per=handles.FileName_pre_per{handles.refImage};
ref_outFileName=[this_pre_per(1:end-20) '.mat'];

thisoutPathName=handles.PathName_pre_per{handles.refImage};

this_outFullFileName=[thisoutPathName ref_outFileName];

load(this_outFullFileName)

ref_output=output;

ref_stacked_weights=zeros(size(output.spatial_weights,1),size(output.spatial_weights,2));

for ii=1:size(ref_output.spatial_weights,3)
    for xx=1:size(ref_output.spatial_weights,1)
        for yy=1:size(ref_output.spatial_weights,2)
            if output.spatial_weights(xx,yy,ii)>0
                ref_stacked_weights(xx,yy)=ref_stacked_weights(xx,yy)+ref_output.spatial_weights(xx,yy,ii);
            end
        end
    end
end

ref_stacked_weights=ref_stacked_weights/max(max(ref_stacked_weights));

figNo=figNo+1;
try
    close(figNo)
catch
end

%Plot the timecourse
hFig = figure(figNo);

hold on

imshow(ref_stacked_weights)
title(['Reference ' ref_outFileName])
set(hFig, 'units','normalized','position',[.01 .1 .4 .4])

working_ROIs_included_in_reference_table=zeros(handles.no_files,1000); %Number of working ROI for each reference ROI
no_ref_ROIs=size(ref_output.spatial_weights,3);
working_ROIs_included_in_reference_table(handles.refImage,1:no_ref_ROIs)=1:no_ref_ROIs;
all_ROI_spatial_weights=ref_output.spatial_weights;

%Enter the reference image ROIs
for fileNo=handles.first_file:handles.no_files
    if fileNo~=handles.refImage
        %File for alignment
        this_pre_per=handles.FileName_pre_per{fileNo};
        working_outFileName=[this_pre_per(1:end-20) '.mat'];


        this_outFullFileName=[thisoutPathName working_outFileName];

        load(this_outFullFileName)

        working_output=output;

        working_stacked_weights=zeros(size(output.spatial_weights,1),size(output.spatial_weights,2));

        for ii=1:size(working_output.spatial_weights,3)
            for xx=1:size(working_output.spatial_weights,1)
                for yy=1:size(working_output.spatial_weights,2)
                    if output.spatial_weights(xx,yy,ii)>0
                        working_stacked_weights(xx,yy)=working_stacked_weights(xx,yy)+working_output.spatial_weights(xx,yy,ii);
                    end
                end
            end
        end

        working_stacked_weights=working_stacked_weights/max(max(working_stacked_weights));

        figNo=figNo+1;
        try
            close(figNo)
        catch
        end

        hFig = figure(figNo);

        hold on

        imshow(working_stacked_weights)
        title(['Working ' working_outFileName])
        set(hFig, 'units','normalized','position',[.01 .1 .4 .4])

%         %Show the overlap in an RGB image
%         overlap_image=zeros(size(working_stacked_weights,1),size(working_stacked_weights,2),3);
%         overlap_image(:,:,1)=ref_stacked_weights;
%         overlap_image(:,:,2)=working_stacked_weights;
% 
%         figNo=figNo+1;
%         try
%             close(figNo)
%         catch
%         end
% 
%         hFig = figure(figNo);
% 
%         hold on
% 
%         imshow(overlap_image)
%         title(['Combined '])
%         set(hFig, 'units','normalized','position',[.01 .1 .4 .4])
%         figure(figNo)
%         title([working_outFileName(1:15) ' vs. ' ref_outFileName(1:15)])

        %Now assign each ROI in the working file with an existing ROI
        for ii=1:size(working_output.spatial_weights,3)
            these_spatial_weights=zeros(size(output.spatial_weights,1),size(output.spatial_weights,2));
            these_spatial_weights(:,:)=working_output.spatial_weights(:,:,ii);
            positive_working_pixels=sum(these_spatial_weights(:)>0);
            %Is this ROI in the reference list?
            found_ROI=0;
            for ii_all=1:size(all_ROI_spatial_weights,3)
                these_all_spatial_weights=zeros(size(output.spatial_weights,1),size(output.spatial_weights,2));
                these_all_spatial_weights(:,:)=all_ROI_spatial_weights(:,:,ii_all);
                overlapping_pixels=sum((these_all_spatial_weights(:).*these_spatial_weights(:))>0);
                positive_all_pixels=sum(these_all_spatial_weights(:)>0);
                mean_positive_pixels=mean([positive_all_pixels positive_working_pixels]);
                if (overlapping_pixels/mean_positive_pixels)>overlap_threshold
                    found_ROI=1;
                    these_working_ROIs_included_in_reference_table=zeros(1,1000);
                    these_working_ROIs_included_in_reference_table(1,:)=working_ROIs_included_in_reference_table(fileNo,:);
                    if sum(these_working_ROIs_included_in_reference_table==ii)>0
                        ii_allb=find(these_working_ROIs_included_in_reference_table==ii);
                        these_allb_spatial_weights=zeros(size(output.spatial_weights,1),size(output.spatial_weights,2));
                        these_allb_spatial_weights(:,:)=all_ROI_spatial_weights(:,:,ii_allb);
                        overlapping_pixelsb=sum((these_allb_spatial_weights(:).*these_spatial_weights(:))>0);
                        %Assign to the ROI with the largest number of overallping pixels
                        if overlapping_pixels>overlapping_pixelsb
                            working_ROIs_included_in_reference_table(fileNo,ii_allb)=0;
                            working_ROIs_included_in_reference_table(fileNo,ii_all)=ii;
                        end
                    else
                        working_ROIs_included_in_reference_table(fileNo,ii_all)=ii;
                    end
                end
            end
            if found_ROI==0
                %The new ROI did not overlap with existing ROIs. Add ROI to all ROIs
                no_ref_ROIs=no_ref_ROIs+1;
                all_ROI_spatial_weights(:,:,no_ref_ROIs)=these_spatial_weights;
                working_ROIs_included_in_reference_table(fileNo,no_ref_ROIs)=ii;
            end
        end

    end
end

%Now do the accounting
ROI_included_in_no_files=zeros(1,handles.no_files);
for ii_all=1:no_ref_ROIs
    ROI_included_in_no_files(sum(working_ROIs_included_in_reference_table(:,ii_all)~=0))=ROI_included_in_no_files(sum(working_ROIs_included_in_reference_table(:,ii_all)~=0))+1;
end


%Show the ROIs occuring in different numbers of recordings
for ii_overlap_figs=1:2
    overlap_image=zeros(size(working_stacked_weights,1),size(working_stacked_weights,2),3);
    for ii_all=1:no_ref_ROIs
        if sum(working_ROIs_included_in_reference_table(:,ii_all)~=0)==1+3*(ii_overlap_figs-1)
            overlap_image(:,:,1)=overlap_image(:,:,1)+all_ROI_spatial_weights(:,:,ii_all);
        end
        if sum(working_ROIs_included_in_reference_table(:,ii_all)~=0)==2+3*(ii_overlap_figs-1)
            overlap_image(:,:,2)=overlap_image(:,:,2)+all_ROI_spatial_weights(:,:,ii_all);
        end
        if sum(working_ROIs_included_in_reference_table(:,ii_all)~=0)==3+3*(ii_overlap_figs-1)
            overlap_image(:,:,3)=overlap_image(:,:,3)+all_ROI_spatial_weights(:,:,ii_all);
        end
    end


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    hold on

    imshow(overlap_image)
    
    set(hFig, 'units','normalized','position',[.01 .1 .4 .4])
    figure(figNo)
    title(['Combined ' num2str(ii_overlap_figs) ' r=' num2str(1+3*(ii_overlap_figs-1)) ' g=' num2str(2+3*(ii_overlap_figs-1)) ' b=' num2str(3+3*(ii_overlap_figs-1))])
end

duplicate_assignemnts=0;
for fileNo=handles.first_file:handles.no_files
    these_aligned_ROIIs=zeros(1,size(working_ROIs_included_in_reference_table,2));
    these_aligned_ROIIs(1,:)=working_ROIs_included_in_reference_table(fileNo,:);
    for ii=1:no_ref_ROIs
        if sum(these_aligned_ROIIs==ii)>1
            duplicate_assignemnts=1;
            fprintf(1,['Multiple aligned ROIs assigned for file number ' num2str(fileNo) ' ROI number ' num2str(ii) ' number of assignemts ' num2str(sum(these_aligned_ROIIs==ii)), '\n'])
        end
    end
end

if duplicate_assignemnts==0
            fprintf(1,'No duplicate assignments\n')
end
save([choiceBatchPathName choiceFileName(1:end-2) handles.suffix '.mat'],'working_ROIs_included_in_reference_table','ROI_included_in_no_files')






