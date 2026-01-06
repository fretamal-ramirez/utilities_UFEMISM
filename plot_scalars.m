% plot ice volumes
path_init='/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20_FCMP/scalar_output_ANT_00001.nc';
path_PD='/Users/frre9931/Desktop/tetralith_results/results_ant_PD_control5000_FCMP/scalar_output_ANT_00001.nc';
%path_PD='/Users/frre9931/Desktop/UFEMISM2.0_porting/results_ant_PD_control5000_PMP/scalar_output_ANT_00001.nc'; 

ice_volume_init=ncread(path_init,'ice_volume');
ice_volume_af_init=ncread(path_init,'ice_volume_af');
time_init=ncread(path_init,'time');

ice_volume_PD=ncread(path_PD,'ice_volume');
ice_volume_af_PD=ncread(path_PD,'ice_volume_af');
time_PD=ncread(path_PD,'time');

% convert PD time starting from 1980 to 50,000 to plot both together
starting_yr=time_init(end);
yrs_of_PD_simulation=time_PD(end)-time_PD(1);
yrs_continued=linspace(starting_yr,starting_yr+yrs_of_PD_simulation,length(time_PD))';

% plot
figure()
plot(time_init,ice_volume_af_init,yrs_continued,ice_volume_af_PD,'LineWidth',2);
grid on
xlabel('time (yrs)');
ylabel('ice volume af');
title('Simulation using FCMP')
% 
% figure()
% plot(time_init,ice_volume_init,yrs_continued,ice_volume_PD,'LineWidth',2);
% grid on
% xlabel('time (yrs)');
% ylabel('ice volume');

%% something similar but now to calculate the volume_af difference between control and retreat
folder_ctrl='results_ant_PD_maxphi_20_HR_ctrl2500_SMB_and_phi_50percent';
folder_run='results_ant_PD_maxphi_20_HR_retreat_SMB_and_phi_50percent';
%path_ctrl = '/Users/frre9931/Desktop/tetralith_results/results_ant_PD_maxphi_20_ctrl2500_SMB_and_phi_50percent/scalar_output_ANT_00001.nc';
path_ctrl = ['/Users/frre9931/Desktop/tetralith_results/',folder_ctrl,'/scalar_output_ANT_00001.nc'];
path_run  = ['/Users/frre9931/Desktop/tetralith_results/',folder_run,'/scalar_output_ANT_00001.nc'];
PD_SL = 55.7; % extracted from UFEMISM output ice_volume_af_PD

ice_volume_ctrl=ncread(path_ctrl,'ice_volume');
ice_volume_af_ctrl=ncread(path_ctrl,'ice_volume_af');

ice_volume_run=ncread(path_run,'ice_volume');
ice_volume_af_run=ncread(path_run,'ice_volume_af');

final_sea_level_contribution = ice_volume_af_ctrl(end) - ice_volume_af_run(end);

%% plot using time series of every simulation
% the time-step is not the same in all simulations, this means that we need
% to plot each of them with their relative time.

% First idea: easiest one, read ice volume af and time, plot them with a
% for

% ==== DEFINE OUTPUTS ====
outputs = { ...
    'results_ant_PD_maxphi_30_SHR_ctrl2500',...
    'results_ant_PD_maxphi_30_SHR_ctrl2500_SMB_and_phi_10percent',...
    'results_ant_PD_maxphi_30_SHR_ctrl2500_SMB_and_phi_20percent',...
    'results_ant_PD_maxphi_30_SHR_ctrl2500_SMB_and_phi_30percent',...
    'results_ant_PD_maxphi_30_SHR_ctrl2500_SMB_and_phi_40percent',...
    'results_ant_PD_maxphi_30_SHR_ctrl2500_SMB_and_phi_50percent',...
    'results_ant_PD_maxphi_30_SHR_retreat',...
    'results_ant_PD_maxphi_30_SHR_retreat_SMB_and_phi_10percent',...
    'results_ant_PD_maxphi_30_SHR_retreat_SMB_and_phi_20percent',...
    'results_ant_PD_maxphi_30_SHR_retreat_SMB_and_phi_30percent',...
    'results_ant_PD_maxphi_30_SHR_retreat_SMB_and_phi_40percent',...
    'results_ant_PD_maxphi_30_SHR_retreat_SMB_and_phi_50percent',...
};
legend_name = { ...
    'control',...
    'SMB&ϕ 10%',...
    'SMB&ϕ 20%',...
    'SMB&ϕ 30%',...
    'SMB&ϕ 40%',...
    'SMB&ϕ 50%', ...
    'retreat',...
    'retreat SMB&ϕ 10%', ...
    'retreat SMB&ϕ 20%', ...
    'retreat SMB&ϕ 30%', ...
    'retreat SMB&ϕ 40%', ...
    'retreat SMB&ϕ 50%', ...
};

basepath = '/Users/frre9931/Desktop/tetralith_results/';

% define SL contribtion for PD AIS, extracted from UFEMISM
PD_SL = 55.7; 

% figure config
fig_width  = 500;
fig_height = 300;

% limits for plot
xmin=1980; xmax=2500;
ymin=-25; ymax=350;

%start figure
H.fig = figure('Units','pixels','Position',[100 100 fig_width fig_height],'Visible','on');

%set colors inspired in ColorBrewer
colors = [27 158 119;
    217 95 2;
    117 112 179;
    231 41 138;
    102 166 30;
    230 171 2;
    166 118 29;
    102 102 102] / 255;

% try this colors. https://jfly.uni-koeln.de/color/#cudo
PALETTE4 = [
    230 159   0;
    128   0 128;
    240 228  66;
      0 114 178;
    213  94   0;
     86 180 233;
    143  72   0;
      0 158 115;
    204 121 167
] / 255;

palette_jfly = [
    102 102 102;
    230 159   0;
    86  180 233;
    0   158 115;
    240 226 66;
    0   114 178;
    213  94   0;
    204 121 167;
] / 255;

plots_with_same_style=6;

% Use tiledlayout with minimal padding and spacing
tl = tiledlayout('flow'); % 'flow' creates a single tile for one plot
tl.Padding = 'compact';
tl.TileSpacing = 'compact';
nexttile
% loop for each UFEMISM scalar output
for i=1:length(outputs)
    output_folder = outputs{i};
    filepath = fullfile(basepath, output_folder, 'scalar_output_ANT_00001.nc');
    ice_volume_af=ncread(filepath,'ice_volume_af');
    time = ncread(filepath,"time");
    
    % sea level contribution respect to PD
    SLC_respect_PD = PD_SL - ice_volume_af ;
    % sea level contribution respect to initial state
    % so it will start from 0
    SLC_respect_init = ice_volume_af(1) - ice_volume_af;
    % plot time series, *1000 to show it in mm instead of meters
    if i <= plots_with_same_style % it was set to 8 before
        plot(time, SLC_respect_init*1000,'LineWidth',2,'Color',palette_jfly(i,:))
    else
        plot(time, SLC_respect_init*1000,'LineWidth',2,'LineStyle','--','Color',palette_jfly(i-plots_with_same_style,:))
    end
    xlabel('Time (yr)');
    ylabel('SLR contribution (mm)');
    grid on
    hold on
end
axis([xmin xmax ymin ymax]);
%legend(legend_name,'Location','northwest');

% ==== LEGEND FOR EXPERIMENTS (COLORS ONLY) ====
hexp = gobjects(plots_with_same_style,1);
for i = 1:plots_with_same_style % set to 8 before
    hexp(i) = plot(NaN, NaN, ...
        'Color', palette_jfly(i,:), ...
        'LineWidth', 2);
end

lgd_exp = legend(hexp, legend_name(1:plots_with_same_style), ...
    'Location', 'northwest');
%lgd_exp.Box = 'off';
lgd_exp.AutoUpdate = 'off';

ax = gca;
pos = ax.Position;

ax2 = axes('Position',[pos(1)+0.29 pos(2)+0.73 0.07 0.05]);
plot([0 1],[0.7 0.7],'k-','LineWidth',2); hold on
plot([0 1],[0.3 0.3],'k--','LineWidth',2)

% text(1.1,0.7,'max\phi = 20°','FontSize',10)
% text(1.1,0.3,'max\phi = 30°','FontSize',10)
text(1.1,0.7,'control','FontSize',10)
text(1.1,0.3,'retreat','FontSize',10)

axis off

print(H.fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/scalarSLE.png','-dpng','-r300');