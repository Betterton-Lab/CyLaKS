
vel = [430, 210, 120, 69, 54;
       400, 76, 30, 21, 17];

vel = [430, 215, 122, 69, 53.5;
       400, 76.3, 29.7, 21.1, 17.5];

run = [14.7, 10.3, 7.50, 5.19, 4.29;
       14.7, 5.78, 2.70, 1.86, 1.35];

n_proto = [1, 2, 3, 5, 8];

vN = @(N) 5473 / (9.49*N + 1.205); 

%set(0, 'DefaultAxesFontName', 'Arial');
%set(0, 'DefaultTextFontName', 'Arial');

% velocity plot
fig1 = figure('Position', [50 50 600 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
hold on
plot(n_proto, vel, '.', 'MarkerSize', 50)
fplot(vN, [1 8], 'LineWidth', 3,'Color', 'k');

ylabel("Average MAP velocity (nm/s)");
ylim([0 550])
yticks([0 100 200 300 400 500])
xlabel("Protofilament Number")
xlim([0 10]);
xticks([1 2 3 5 8]);

legendLabel = ["50 nM motor", "10 nM motor", "Analytic (50 nM)"];
legend(legendLabel, "Location", "northeast")

set(gca,'box','off')
set(gca, 'FontSize', 20);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);

% run plot
fig2 = figure('Position', [50 50 600 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
hold on
plot(n_proto, run, '.', 'MarkerSize', 50)

ylabel("Average MAP displacement (um)");
ylim([0 16.5])
yticks([0 3 6 9 12 15])
xlabel("Protofilament Number")
xlim([0 10]);
xticks([1 2 3 5 8]);

legendLabel = ["50 nM motor", "10 nM motor"]; %, "Analytic (50 nM)"];
legend(legendLabel, "Location", "northeast")

set(gca,'box','off')
set(gca, 'FontSize', 20);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);