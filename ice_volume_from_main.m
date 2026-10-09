% ice volumes from main_output_UFE. only grounded for now

clear all; close all; clc;

path_to_files = '/Users/frre9931/Desktop/arrhenius_results/';

folders = {...
    'results_R-LIS_init_BCAP_H1_O1_SMB1_GHF1_2km_ext',...
    'results_R-LIS_init_STDdH_H1_O1_SMB1_GHF1_2km_extended',...
    'results_R-LIS_init_STD_H1_O1_SMB1_GHF1_SL2_2km_ext',...
    'results_direct_R-LIS_init_BCAP_H2_O1_SMB1_GHF1_SL1_2km',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_M11_Hb-2000to-250m_SHR',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_max30_SHR_new',...
    'results_R-LIS_init_BCAP_H1_O1_SMB1_GHF1_2km_MS1mod',...;
     };
plot_names={...
    'BCAP',...
    'STD dH/dt',...
    'STD Budd',...
    'BCAP H2',...
    'BCAP MS1',...
    'STD MS1',...
    'BCAP MS1 40kyr',...
};

% Ocean area
A_ocean = 3.618e14; % m^2, global ocean area
area_pixel = 5000*5000; % 5 km res

figure()

for i= 1:length(folders)
file_to_read = fullfile(path_to_files,folders{i},'/main_output_ANT_grid_ROI_RiiserLarsen.nc');
time = ncread(file_to_read,'time');
Hi = ncread(file_to_read,'Hi');

ice_vol = zeros(size(time));
for k = 1:size(Hi,3) % iterate over time
    Hi_pixel = Hi(:,:,k)*area_pixel;
    ice_vol(k) = sum(sum(Hi_pixel));
end
ice_vol_sle = 910/1028 * ice_vol/A_ocean;

plot(time/1000,ice_vol_sle,'LineWidth',2)
hold on
end
legend(plot_names);
grid on
ylabel('Ice volume (m)')
xlabel('Time (kyr)')
title('Ice volume calculated from gridded output for ROI km resolution')