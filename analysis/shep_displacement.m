
clear variables;

seeds = [0, 1, 2, 3, 4, 5];
%seeds = [0];

site_SID = 0;
xlink_SID = 1;
motor_SID = 2;
endtag_boundary = 375;
flux_win_size = 500;
vel_limit_upper = 3000;
vel_limit_lower = 15;

xlink_analysis_only = true;

output_folder = 'plots_geoVel';
file_dir = '../out_final/';

%sim_name_base = 'out_final_xlinkOnly/shep_0.1nM_0.0nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final/shep_1nM_100nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final_xlinkOnly/shep_0.1nM_0.0nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final/shep_0.1nM_100nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final/shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0'
%sim_name_base = '/out_final_xlinkLifetime/shep_0.01nM_100nM_8_1000_0.6kT_3x_5x_0_xlink_0.1x'
sim_name_base = 'shep_1nM_100nM_8_1000_0.6kT_3x_5x';

% Open log file and parse it into param labels & their values
log_file = sprintf('%s/%s', file_dir, sprintf('%s_0', sim_name_base));
params = load_parameters(log_file);
n_seeds = length(seeds);

% motor ID is unique; make following arrays serial w/ ID as index
runlengths = zeros([2 100 * n_seeds * params.n_mts * params.max_sites]);
lifetimes = zeros([2 100 * n_seeds  * params.n_mts * params.max_sites]);
velocities = zeros([2 100 * n_seeds  * params.n_mts * params.max_sites]);
n_runs = zeros([2 1]);
n_runs_vel = zeros([2 1]);

for i_seed = 1:n_seeds
    sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed))
    %sim_name = sim_name_base;

    proteinFileStruct = '%s_protein_id.file';
    proteinFileName = sprintf("%s/%s", file_dir, sprintf(proteinFileStruct, sim_name));
    protein_data_file = fopen(proteinFileName);
    raw_motor_data = fread(protein_data_file, params.n_mts * params.max_sites * params.n_datapoints, '*int');
    protein_data = reshape(raw_motor_data, params.max_sites, params.n_mts, params.n_datapoints);
    fclose(protein_data_file);

    occuFileStruct = '%s_occupancy.file';
    occuFileName = sprintf("%s/%s", file_dir, sprintf(occuFileStruct, sim_name));
    occu_data_file = fopen(occuFileName);
    raw_occu_data = fread(occu_data_file, params.n_mts * params.max_sites * params.n_datapoints, '*int');
    occupancy_data = reshape(raw_occu_data, params.max_sites, params.n_mts, params.n_datapoints);
    fclose(occu_data_file);

    % have an active list for each MT
    active_proteins = zeros([2 params.n_mts params.n_mts * params.max_sites * 10]);
    n_active = zeros([2 params.n_mts]);
    starting_site = zeros([2 100 * params.n_mts * params.max_sites]) - 1;
    starting_mt = zeros([2 100 * params.n_mts * params.max_sites]) - 1;
    starting_datapoint = zeros([2 100 * params.n_mts * params.max_sites]) - 1;
    for i_data = 1 : params.n_datapoints - 1
        for i_mt = 1: 1 : params.n_mts
            protein_IDs = protein_data(:, i_mt, i_data);
            % Scan through IDs of bound motors (-1 means no motor on that site)
            for i_site = 1 : 1 : params.max_sites
                i_species = occupancy_data(i_site, i_mt, i_data);
                if i_species == site_SID || (xlink_analysis_only && i_species == motor_SID)
                    continue;
                end
                protein_ID = protein_IDs(i_site);
                if protein_ID <= 0
                    disp("err");
                    return
                end
                % Always count protein on first site
                if i_site == params.max_sites
                    % Record the protein's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(i_species, protein_ID) == -1
                        starting_site(i_species, protein_ID) = i_site;
                        starting_datapoint(i_species, protein_ID) = i_data;
                        starting_mt(i_species, protein_ID) =  i_mt;
                        n_active(i_species, i_mt) = n_active(i_species, i_mt) + 1;
                        active_proteins(i_species, i_mt, n_active(i_species, i_mt)) = protein_ID;
                    end
                    % Otherwise if a motor is found, only count first
                elseif protein_IDs(i_site + 1) ~= protein_ID
                    % Record the motor's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(i_species, protein_ID) == -1
                        starting_site(i_species, protein_ID) = i_site;
                        starting_datapoint(i_species, protein_ID) = i_data;
                        starting_mt(i_species, protein_ID) =  i_mt;
                        n_active(i_species, i_mt) = n_active(i_species, i_mt) + 1;
                        active_proteins(i_species, i_mt, n_active(i_species, i_mt)) = protein_ID;
                    end
                end
            end
            % Check one datapoint into the future to see if any motors unbound
            for i_species = 1 : 1 : 2
                n_deleted = 0;
                for i_protein = 1:1:n_active(i_species, i_mt)
                    i_adj = i_protein - n_deleted;
                    protein_ID = active_proteins(i_species, i_mt, i_adj);
                    unbound = true;
                    for j_mt = i_mt : 1 : params.n_mts + i_mt - 1
                        if j_mt > params.n_mts
                            j_mt_adj = j_mt - params.n_mts;
                        else
                            j_mt_adj = j_mt;
                        end
                        future_site = find(protein_data(:, j_mt_adj, i_data + 1 ) == protein_ID, 1);
                        if ~isempty(future_site)
                            % motors have high turnover rate; make sure they do not unbind+rebind in single step
                            if i_species == motor_SID
                                if j_mt == starting_mt(i_species, protein_ID) && future_site(1) < starting_site(i_species, protein_ID)
                                    unbound = false;
                                end
                            else
                                unbound = false;
                            end
                        end
                    end
                    if unbound
                        % Calculate distance traveled
                        end_site = -1;
                        end_mt = -1;
                        for j_mt = i_mt : 1 : params.n_mts + i_mt - 1
                            if j_mt > params.n_mts
                                j_mt_adj = j_mt - params.n_mts;
                            else
                                j_mt_adj = j_mt;
                            end
                            end_site = find(protein_data(:, j_mt_adj, i_data) == protein_ID, 1);
                            if ~isempty(end_site)
                                end_mt = j_mt_adj;
                                break;
                            end
                        end
                        if isempty(end_site)
                            disp('Error in finding end site');
                            return;
                        end
                        start_site = starting_site(i_species, protein_ID);
                        delta = end_site(1) - start_site; % we want actual displacement, not MSD
                        run_length = delta * params.site_size;
                        % Calculate time bound
                        start_datapoint = starting_datapoint(i_species, protein_ID);
                        delta_time = i_data - start_datapoint;
                        run_time = delta_time * params.time_per_datapoint;
                        velocity = (run_length / run_time) * 1000; % convert to nm/s
                        % If time bound is above time cutoff, add to data
                        if run_time > params.time_per_datapoint && abs(velocity) > vel_limit_lower && abs(velocity) < vel_limit_upper
                            n_runs(i_species) = n_runs(i_species) + 1;
                            runlengths(i_species, n_runs(i_species)) = run_length;
                            lifetimes(i_species, n_runs(i_species)) = run_time;
                            if end_site(1) > endtag_boundary
                                n_runs_vel(i_species) = n_runs_vel(i_species) + 1;
                                velocities(i_species, n_runs_vel(i_species)) = velocity;
                            end
                        end
                        starting_site(i_species, protein_ID) = -1;
                        starting_datapoint(i_species, protein_ID) = -1;
                        starting_mt(i_species, protein_ID) = -1;
                        % Switch now-deleted entry with last entry in active_motors
                        active_proteins(i_species, i_mt, i_adj) = active_proteins(i_species, i_mt, n_active(i_species, i_mt));
                        active_proteins(i_species, i_mt, n_active(i_species, i_mt)) = -1;
                        n_active(i_species, i_mt) = n_active(i_species, i_mt) - 1;
                        n_deleted = n_deleted + 1;
                    end
                end
            end
        end
    end
end



% trim arrays to get rid of un-used containers
%runlengths = runlengths(1:n_runs);
%lifetimes = lifetimes(1:n_runs);
%velocities = velocities(1:n_runs);
%}

%disp(velocities(find(velocities <= 0)))
%[avg_run, err_run] = get_stats(xlink_SID, runlengths, n_runs, 'normal', "Run length (um)", output_folder, sim_name, "run");
%[avg_time, err_time] = get_stats(xlink_SID, lifetimes, n_runs, 'exponential', "Lifetime (s)", output_folder, sim_name, "time");
%[avg_vel, err_vel] = get_stats(xlink_SID, velocities, n_runs_vel, 'normal', "Velocity (nm/s)", output_folder, sim_name, "vel");
[avg_vel_geo, err_vel_geo] = get_stats(xlink_SID, abs(velocities), n_runs_vel, 'lognormal', "Geometric velocity (nm/s)", output_folder, sim_name, "velGeo");

%fprintf("Run: %.2f ± %.2f | Time: %.2f ± %.2f | Vel: %.2f ± %.2f (%.2f ± %.2f)\n", avg_run, err_run, avg_time, err_time, avg_vel, err_run, avg_vel_geo, err_vel_geo)
%}

fprintf("geo vel: %.2f ± %.2f \n", avg_vel_geo, err_vel_geo)



function [avg, sem] = get_stats(SID, data_arr, n_runs_arr, distribution, dataLabel, output_folder, sim_name, plotname)
n_runs = n_runs_arr(SID);
data = data_arr(SID, 1:n_runs);
if strcmp(distribution, 'lognormal') == 1
    avg = exp(sum(log(data)) / n_runs);
    sd = exp(std(log(data)));
    sem = sd / sqrt(n_runs);
else
    avg = sum(data) / n_runs;
    sd = std(data);
    sem = sd / sqrt(n_runs);
end
if n_runs < 10
    disp("Not enough statistics for histogram fitting")
    return
end
if strcmp(distribution, 'exponential') == 1
    min_val = min(data);
    data = data - min_val;
end
pd = fitdist(data', distribution);
conf_inv = paramci(pd);
if strcmp(distribution, 'normal') == 1
    mu = pd.mu;
    sigma = pd.sigma;
    stderr = (conf_inv(2) - conf_inv(1)) / (2*1.96); % get std_err
elseif strcmp(distribution, 'exponential') == 1
    mu = pd.mu + min_val;
    sigma = nan;
    stderr = (conf_inv(2) - conf_inv(1)) / (2*1.96); % get std_err
elseif strcmp(distribution, 'lognormal') == 1
    mu = exp(pd.mu);
    sigma = exp(pd.sigma);
    stderr = (exp(conf_inv(2)) - exp(conf_inv(1))) / (2*1.96); % get std_err
else
    disp("unrecognized probability distribution")
    return
end
if abs(avg - mu) > 1e-6
    fprintf("Error calculating %s for species %i\n", dataLabel, SID)
    avg = nan;
    sem = nan;
    return;
end
switch plotname
    case 'run'
        bin_lower_lim = 0.1;
    case 'time'
        bin_lower_lim = 0.2;
    case 'vel'
        bin_lower_lim = 50;
    case 'velGeo'
        bin_lower_lim = 50;
    otherwise
        disp('error')
        return
end
if SID == 1
    speciesLabel = "xlink";
elseif SID == 2
    speciesLabel = "motor";
end
fig = figure('Position', [50 50 720 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
if strcmp(plotname, "velGeo") == 1
    bin_width = 25;
else
    bin_width = max(bin_lower_lim, 2*iqr(data)/(n_runs^(1/3)));
end
n_bins = ceil((max(data) - min(data))/bin_width);
h = histfit(data, n_bins, distribution);
h(1).FaceAlpha = 0.5;
h(1).EdgeAlpha = 0.5;
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(1).BarWidth = 0.75;
h(2).LineWidth = 3;
if strcmp(plotname, "velGeo") == 1
    xlim([0 400])
    %ylim([0 200])
else
    xlim([min(0, min(data)) max(data)])
end
%xlabel(sprintf("%s of %ss", dataLabel, speciesLabel));
xlabel( "Velocity (nm/s)");
ylabel("Frequency");
dim1 = [0.625 0.7 0.2 0.2];
str1 = sprintf('Avg = %#.3g ± %#.2g', avg, sem);
dim2 = [0.625 0.65 0.2 0.2];
str2 = sprintf('Std = %#.3g', sd);
dim3 = [0.625 0.6 0.2 0.2];
str3 = sprintf('Mu = %#.3g ± %#.2g', mu, stderr);
dim4 = [0.625 0.55 0.2 0.2];
str4 = sprintf('Sigma = %#.3g', sigma);
dim5 = [0.625 0.5 0.2 0.2];
str5 = sprintf('N = %i', n_runs);
dim6 = [0.625 0.45 0.2 0.2];
str6 = sprintf('Bin width = %#.3g', bin_width);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim4, 'String', str4, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim5, 'String', str5, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim6, 'String', str6, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
set(gca,'box','off')
set(gca, 'FontSize', 28);
set(gca, 'FontName', 'Arial');
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title(sprintf("%s (%ss)", sim_name, speciesLabel),'Interpreter', 'none', "FontName", "Arial", "FontSize", 10);
exportgraphics(fig, sprintf('%s/%s_%s_%s.png', output_folder, sim_name, speciesLabel, plotname),'ContentType','image');
exportgraphics(fig, sprintf('%s/%s_%s_%s.pdf', output_folder, sim_name, speciesLabel, plotname),'ContentType','vector');
saveas(fig, sprintf('%s/%s_%s_%s.svg', output_folder, sim_name, speciesLabel, plotname),'svg');
%saveas(fig, sprintf('%s/%s_%s_%s.png', output_folder, sim_name, speciesLabel, plotname),'png');
%close(fig);
%fprintf("%.3g +/- %.1g\n", mu, stderr);
%fprintf("%.3g +/- %.1g\n", avg, sem);
end
