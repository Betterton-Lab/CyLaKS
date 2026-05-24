clear variables;

sim_name = 'out_final_johann/shep_0.504nM_33.52nM_1_500_0.0kT_0.01xCeff_0';
sim_name = 'out_final_johann/shep_johann_simpleMotors';

size_x = 500;
size_y = 400;
size_dataWindow = 18000;

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

n_datapoints = int32(params.n_datapoints);

% Open occupancy data file
occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
occupancy = load_data(occupancy, occupancy_filename, '*int');

xlink_speciesID = 1;
motor_speciesID = 2;
altMAP_speciesID = 3;

xlink_raw_data = occupancy;
motor_raw_data = occupancy;

xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;
motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;

xlink_avg_occupancy = zeros([params.max_sites params.n_mts]);
motor_avg_occupancy = zeros([params.max_sites params.n_mts]);

motor_avg_occupancy_tot = zeros([params.max_sites 1]);
xlink_avg_occupancy_tot = zeros([params.max_sites 1]);

fig1 = figure('Position', [50, 250, size_x, size_y]);

for i = n_datapoints-size_dataWindow:1:params.n_datapoints
    for i_pf = 1 : 1 : params.n_mts
        motor_avg_occupancy(:, i_pf) = motor_avg_occupancy(:, i_pf) + double(motor_raw_data(:, i_pf, i)) ./ size_dataWindow;
        xlink_avg_occupancy(:, i_pf) = xlink_avg_occupancy(:, i_pf) + double(xlink_raw_data(:, i_pf, i)) ./ size_dataWindow;
        motor_avg_occupancy_tot(:, 1) = motor_avg_occupancy_tot(:, 1) + double(motor_raw_data(:, i_pf, i)) ./ (size_dataWindow * params.n_mts);
        xlink_avg_occupancy_tot(:, 1) = xlink_avg_occupancy_tot(:, 1) + double(xlink_raw_data(:, i_pf, i)) ./ (size_dataWindow * params.n_mts);
    end
end

%{
smooth_window = 32; % should be equivalent to diffraction limit
motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
motor_occupancy_tot = smoothdata(motor_avg_occupancy_tot, 'movmean', smooth_window);
xlink_occupancy_tot = smoothdata(xlink_avg_occupancy_tot, 'movmean', smooth_window);
%}

xlink_avg_occupancy_tot = flip(xlink_avg_occupancy_tot);
motor_avg_occupancy_tot = flip(motor_avg_occupancy_tot);

yyaxis left
plot(linspace(0, params.max_sites, params.max_sites), ...
    xlink_avg_occupancy_tot, 'LineWidth', 2.5);
ylabel('MAP average occupancy', 'FontSize', 14);       
ylim([0 1]);
yticks([0.25 0.5 0.75]);

yyaxis right
plot(linspace(0, params.max_sites, params.max_sites), ...
    motor_avg_occupancy_tot, 'LineWidth', 2.5); 
ylabel('Motor average occupancy', 'FontSize', 14);  
ylim([0 0.1]);
yticks([0.025 0.05 0.075]);
xlabel('Site index', 'FontSize', 14);


%five_percent = params.max_sites * params.site_size / 20.0;
%xlim([-five_percent params.max_sites * params.site_size + five_percent]);
set(gca, 'FontSize', 14);
axis = gca;
%axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
%legendLabel = {'Xlinks (avg)', 'Motors(avg)'};
%legend(legendLabel, 'Location', 'northeastoutside');
%pbaspect([1 1 1]);
