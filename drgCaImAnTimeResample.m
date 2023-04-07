function out_timecourse = drgCaImAnTimeResample(time_span,in_timecourse,resample_time_span)
%Resamples the in_variable to resample_time_span
%   Detailed explanation goes here
if length(time_span)>length(in_timecourse)
    time_span=time_span(1:length(in_timecourse));
end
out_timecourse=zeros(1,length(resample_time_span));
for ii_t=1:length(resample_time_span)
    this_jj_t_after=find(time_span>=resample_time_span(ii_t),1,'first');
    this_jj_t_before=find(time_span<resample_time_span(ii_t),1,'last');
    if isempty(this_jj_t_before)
        out_timecourse(ii_t)=in_timecourse(this_jj_t_after);
    else
       if isempty(this_jj_t_after)
           out_timecourse(ii_t)=in_timecourse(this_jj_t_before);
       else
        %Interpolate
            out_timecourse(ii_t)=in_timecourse(this_jj_t_before)+(resample_time_span(ii_t)-time_span(this_jj_t_before))*...
                (in_timecourse(this_jj_t_after)-in_timecourse(this_jj_t_before))/(time_span(this_jj_t_after)-time_span(this_jj_t_before));
       end
    end
end