%create a transient dT_ocean forcing to use as input in UFEMISM2.0
clear all
% let's start from 2000 to 2500
time = 2000:1:2500;
time_rel = time - time(1);
% dT to perturb the ocean field
dT = 2; 
%warming_rate = dT/ (time(end)-time(1)); % C yr-1
warming_rate = 4e-3; % C yr-1
warming_rate2 = 8e-3;
% create a vector with the values
dT_constant = ones(size(time)).*dT;
%dT_ramp = linspace(0,dT,length(time));

dT_ramp = zeros(size(time));
dT_ramp2 = dT_ramp;
for i= 1:length(time)-1
    % create ramp forcing that increases by a set rate
    dT_ramp(i+1) = dT_ramp(i) + warming_rate;
    dT_ramp2(i+1) = dT_ramp2(i) + warming_rate2;
end
limit_dT = 2.0; % maximum temperature allowed in forcing

% limit maximum value reached to limit_dT
dT_ramp_limited = min( dT_ramp, limit_dT);
dT_ramp_limited2 = min( dT_ramp2, limit_dT);

%find the time when the cap will be reached
time_cap = limit_dT / warming_rate;
time_cap2 = limit_dT / warming_rate2;
% normalized time (0 to 1)
%t_norm = (time - time(1)) / (time(end) - time(1));
% normalized time (0–1 until cap)
t_norm_cap = min(time_rel / time_cap, 1);
t_norm_cap2 = min(time_rel / time_cap2, 1);
% power ramp
p = 2;  %
%dT_power = dT * t_norm.^p;
%cap_dT_power = min(dT_power,limit_dT);

cap_dT_power = limit_dT * (t_norm_cap).^p;
cap_dT_power2 = limit_dT * (t_norm_cap2).^p;
%% plot data
% figure config
fig_width  = 17.9; % 179 mm, two columns JOG
fig_height = 10.0; % 254 mm maximum height JOG

% limits for plot
xmin=2000; xmax=2500;
ymin=0; ymax=2.2; %

H.fig = figure('Visible','off');
set(H.fig,'PaperUnits','centimeters');
set(H.fig,'PaperSize',[fig_width fig_height]);
set(H.fig,'PaperPosition',[0 0 fig_width fig_height]);

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

% Use tiledlayout with minimal padding and spacing
tl = tiledlayout('flow'); % 'flow' creates a single tile for one plot
tl.Padding = 'compact';
tl.TileSpacing = 'compact';
nexttile

plot(time,dT_constant,'LineWidth',2,'Color',palette_jfly(1,:))
hold on
plot(time,dT_ramp_limited,'LineWidth',2,'Color',palette_jfly(2,:),'LineStyle','--')
plot(time,cap_dT_power,'LineWidth',2,'Color',palette_jfly(3,:),'LineStyle','-.')
plot(time,dT_ramp_limited2,'LineWidth',2,'Color',palette_jfly(4,:),'LineStyle','--')
plot(time,cap_dT_power2,'LineWidth',2,'Color',palette_jfly(5,:),'LineStyle','-.')
grid on
xlabel('Year');
ylabel('Ocean temperature offset (°C)');
axis([xmin xmax ymin ymax]);
legend('OC 2.0','LR4e-3','PW4e-3','LR8e-3','PW8e-3','Location','southeast')
print(H.fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/scalardTocns.pdf','-dpdf','-vector');

%% save file as netCDF
% possible names for variable in file 'dT||dT_ocean||dTo'
% filename with path
%filename = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/ocean_forcings/dT_power_6e1_2p_max_1e1.nc';
%filename = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/ocean_forcings/dT_rt_5e-3_max_1e1.nc';
filename = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/ocean_forcings/dT_pw2_rt_4e-3_max_2e1.nc';

% save as netcdf file
nccreate(filename, 'time', 'Dimension', {'time', length(time)});
nccreate(filename, 'dT_ocean', 'Dimensions',{'time', length(dT_constant)});

% Write data
ncwrite(filename, 'dT_ocean', cap_dT_power);
ncwrite(filename, 'time', time);

% Add metadata
ncwriteatt(filename,'time','units','years');
ncwriteatt(filename,'dT_ocean','units','degrees');