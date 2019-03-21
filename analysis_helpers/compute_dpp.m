function dpp = compute_dpp(distance_to_monitor, monitor_width, monitor_resolution)
%% 
% compute degree per pixel
%
% EXAMPLE: 
% dpp = compute_dpp(ex.setup.viewingDistance, ex.setup.monitorWidth,
% ex.setup.screenRect(3:4));
%

dpp = rad2deg(atan2(0.5*monitor_width*(monitor_resolution(2)/monitor_resolution(1)), distance_to_monitor)/(0.5*monitor_resolution(2)));


