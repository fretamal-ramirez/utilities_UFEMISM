% code to plot plot_ice_volumes from scalar_ROI...
clear all; close all; clc;

path_to_files = '/Users/frre9931/Desktop/arrhenius_results/';

folders = {...
    'results_direct_R-LIS_init_BCAP_H2_O1_SMB1_GHF1_SL1_2km',...
    'results_direct_R-LIS_init_BCAP_H1_O1_SMB1_GHF1_SL2_2km',...
    'results_direct_R-LIS_init_STD_H1_O1_SMB1_GHF1_SL1_2km',...
    'results_direct_R-LIS_init_STDdH_H1_O1_SMB1_GHF1_SL1_2km',...
    'results_direct_R-LIS_init_BCAP_H1_O1_SMB2_GHF1_SL1_2km',...
     };
plot_names={...
    'BCAP H2',...
    'BCAP Budd',...
    'STD',...
    'STDdH',...
    'BCAP SMB2',...
};

figure()
tiledlayout(2,1)
nexttile
for i= 1:length(folders)
    file_to_read = fullfile(path_to_files,folders{i},'/scalar_output_ANT_ROI_RiiserLarsen_00001.nc');
    ice_vol = ncread(file_to_read,'ice_volume');
    time = ncread(file_to_read,'time');

    plot(time/1000,ice_vol,'LineWidth',2)
    hold on
end
legend(plot_names);
grid on
title('Ice volume')
ylabel('Ice volume (m)')
%hold off

nexttile
ax = gca;
ax.XAxis.Exponent = 0; % Disables scientific notation on the X-axis
for i= 1:length(folders)
    file_to_read = fullfile(path_to_files,folders{i},'/scalar_output_ANT_ROI_RiiserLarsen_00001.nc');
    ice_vol_af = ncread(file_to_read,'ice_volume_af');
    time = ncread(file_to_read,'time');

    plot(time/1000,ice_vol_af,'LineWidth',2)
    hold on
end
legend(plot_names);
grid on
title('Ice volume above floatation')
ylabel('Ice volume af (m)')
xlabel('Time (kyr)')
%hold off