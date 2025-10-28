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
path_ctrl = '/Users/frre9931/Desktop/tetralith_results/results_ant_PD_maxphi_20_ctrl2500_SMB_and_phi_50percent/scalar_output_ANT_00001.nc';
path_ctrl = '/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20_PMP/scalar_output_ANT_00001.nc';
path_run = '/Users/frre9931/Desktop/tetralith_results/results_ant_PD_maxphi_20_retreat_mask_code_SMB_and_phi_50percent/scalar_output_ANT_00001.nc';
PD_SL = 55.7; % extracted from UFEMISM output ice_volume_af_PD

ice_volume_ctrl=ncread(path_ctrl,'ice_volume');
ice_volume_af_ctrl=ncread(path_ctrl,'ice_volume_af');

ice_volume_run=ncread(path_run,'ice_volume');
ice_volume_af_run=ncread(path_run,'ice_volume_af');

sea_level_contribution = ice_volume_af_ctrl(end) - ice_volume_af_run(end);
