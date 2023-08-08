function [mean_points CI_points] = drgViolinPoint(points,edges,x_val,rand_offset,which_color_mean,which_color_point,point_size)
%This function plots a violin data 
% 
% points has the data to be plotted
% edges are the edges of a histogram encompassing all data
% x_val is the x location
% rand_offset is the width of the violin plot
% which_color is the color for the points e.g. 'k'
% point_size is the size of the points i.e. 1

h=histogram(points,edges,'Visible','off');
normval=h.Values/max(h.Values);
random_offsets=rand(1,length(points))-0.5;

%Plot violin points
for ii=1:length(points)
%     [min_hist ii_histo]=min(abs(edges-points(ii)));
%     if ii_histo>length(normval)
%         ii_histo=length(normval);
%     end
    violin_x_val(ii)=random_offsets(ii)*rand_offset+x_val;
end 
plot(violin_x_val,points,'o','MarkerSize',point_size,'MarkerFaceColor',which_color_mean,'MarkerEdgeColor',which_color_point)
 
%Plot the mean
mean_points=mean(points);
% plot(x_val,mean_points,['o' which_color_mean],'MarkerFaceColor',which_color_mean,'MarkerSize',10)

%Plot CI
CI_points=bootci(1000, @mean, points);
plot([x_val x_val],CI_points,['-' which_color_mean],'LineWidth',3)

pffft=1;
