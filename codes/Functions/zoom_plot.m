%% FUNCTION FOR ZOOMING IN A FIGURE
function zoom_plot(error_RENLMS,error_NLMS,error_LMS,error_ENLMS,start,range)
% This function produces a zoomed in plot of the original plot from the 
% given index of the original plot (defined by 'start'), and the range of 
% indices to be zoomed.
% Inputs  -> plot_params (structure containing quantities to be plotted)
%         -> start (Starting index for zooming)
%         -> range (Range of indices to be zoomed)
% Outputs -> Figure with a zoomed in portion.


figure;
a1 = axes(); % Define an axis handle for figure properties
plot(a1,error_RENLMS);hold on;plot(a1,error_NLMS);hold on;...
    plot(a1,error_LMS);hold on;plot(a1,error_ENLMS)

a2 = axes();  % Define new axis handle

% Define position in the figure window to be zoomed
a2.Position = [0.45 0.45 0.3 0.3]; % [xlocation, ylocation, xsize, ysize]

% Plot the part of the original figure to be zoomed 
endZoom = start + range; % Ending index for zooming
plot(a2,error_RENLMS(start:endZoom));hold on;...
    plot(a2,error_NLMS(start:endZoom)); plot(a2,error_LMS(start:endZoom));...
    plot(a2,error_enlms(start:endZoom));

axis tight; % Hold the current axis handle

% Annotate the zoomed portion with the corresponding shape and position
% co-ordinates
annotation('ellipse',[.75 0.1 .1 .1]) % Shape of the zoomed part
annotation('arrow',[.8 .7],[.1 .4])   % Arrow pointing to the zoomed part
legend(a1,'u')

end
