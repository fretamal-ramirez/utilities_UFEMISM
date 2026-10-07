% code to plot plot_ice_volumes
clear all; close all; clc;

path_to_files = '/Users/frre9931/Desktop/arrhenius_results/';

folders = {...
    'results_R-LIS_init_BCAPdH_H1_O1_SMB1_GHF1_2km_ext',...
    'results_R-LIS_init_BCAP_H1_O1_SMB1_GHF1_2km_ext',...
    'results_R-LIS_init_BCAP_H1_O1_SMB1_GHF2_2km',...
    'results_R-LIS_init_STDdH_H1_O1_SMB1_GHF1_2km_extended',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_M11_Hb-2000to-250m_SHR',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_max30_SHR_new',...
    };
plot_names={...
    'BCAP dH/dt','BCAP','BCAP GHF2','STD dH/dt',...
    'BCAP MS1','STD MS1'
};

figure()
tiledlayout(2,1)
nexttile
for i= 1:length(folders)
    file_to_read = fullfile(path_to_files,folders{i},'/scalar_output_ANT_00001.nc');
    ice_vol = ncread(file_to_read,'ice_volume');
    time = ncread(file_to_read,'time');

    plot(time,ice_vol,'LineWidth',2)
    hold on
end
legend(plot_names);
grid on
title('Ice volume using grounded ice mask protection')
%hold off

nexttile
for i= 1:length(folders)
    file_to_read = fullfile(path_to_files,folders{i},'/scalar_output_ANT_00001.nc');
    ice_vol_af = ncread(file_to_read,'ice_volume_af');
    time = ncread(file_to_read,'time');

    plot(time,ice_vol_af,'LineWidth',2)
    hold on
end
legend(plot_names);
grid on
title('Ice volume above floatation using grounded ice mask protection')
%hold off
