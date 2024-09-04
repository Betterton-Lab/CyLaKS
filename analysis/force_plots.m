clear variables;

sim_name = 'Test_extreme_wca_2'; %'Difp1_N30_sym3_eps3_asym50_off500_dif04_10uM_2';
start_time = 0;          % in seconds 
end_time = -1;           % in seconds; -1 defaults to full data 
smoothing_window  = 1;   % in seconds
didnt_end=0;
% Load parameter structure
file_dir = '../runs';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));
        data_steps_file=fopen(sprintf('%s/%s_%s', file_dir, sim_name,"steps.txt"), 'r');
    if data_steps_file == -1
        %error('File could not be opened. Check the file path and permissions.');
        didnt_end=didnt_end+1;
    else
       number_of_steps = fgetl(data_steps_file);
       %end_frame=str2num(line);
       %params.n_datapoints=str2num(number_of_steps);
       params.n_datapoints=floor(str2num(number_of_steps) / 100) * 100;
    end
% Open data files
xlink_filename = sprintf('%s/%s_xlink_force.file', file_dir, sim_name);
xlink_data = zeros(params.n_dims, params.n_datapoints);
xlink_data = load_data(xlink_data, xlink_filename, 'double');

i_start = int32(start_time / params.time_per_datapoint) + 1;
if(end_time == -1)
    end_time = params.n_datapoints * params.time_per_datapoint;
end
i_end = int32(end_time / params.time_per_datapoint);
n_datapoints_active = i_end - i_start + 1;
smoothing_window_size = smoothing_window / params.time_per_datapoint;

smoothed_x = smooth(xlink_data(1, :), smoothing_window_size);
smoothed_y = smooth(xlink_data(2, :), smoothing_window_size);
smoothed_tot = smooth(sqrt(xlink_data(1, :).*xlink_data(2, :)), smoothing_window_size);

fig1 = figure();
set(fig1, 'Position', [100, 100, 1000, 500]);
hold on
plot(linspace(start_time, end_time, n_datapoints_active), smoothed_x(i_start:i_end), 'LineWidth', 2);
plot(linspace(start_time, end_time, n_datapoints_active), smoothed_y(i_start:i_end), 'LineWidth', 2);
plot([start_time end_time], [0 0], '--', 'LineWidth', 2);
%plot(linspace(start_time, end_time, n_datapoints), smoothed_tot, 'LineWidth', 2);


ylabel('Applied force (pN)');
xlabel('Time (s)');
% cip off boundaries b/c smoothing makes them artificially large 
%xlim([start_time+smoothing_window end_time-smoothing_window]);
%ylim([-2.5 2.5]);
grid off
legend({'Horizontal component', 'Vertical component'}, 'location', 'best', 'FontSize', 12);
legend('boxoff');