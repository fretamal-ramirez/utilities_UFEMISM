% plot ice volumes
ice_volume_init=ncread('/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20/scalar_output_ANT_00001.nc','ice_volume');
ice_volume_af_init=ncread('/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20/scalar_output_ANT_00001.nc','ice_volume_af');
time_init=ncread('/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20/scalar_output_ANT_00001.nc','time');

ice_volume_PD=ncread('/Users/frre9931/Desktop/UFEMISM2.0_porting/results_ant_PD_control5000/scalar_output_ANT_00001.nc','ice_volume');
ice_volume_af_PD=ncread('/Users/frre9931/Desktop/UFEMISM2.0_porting/results_ant_PD_control5000/scalar_output_ANT_00001.nc','ice_volume_af');
time_PD=ncread('/Users/frre9931/Desktop/UFEMISM2.0_porting/results_ant_PD_control5000/scalar_output_ANT_00001.nc','time');

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

figure()
plot(time_init,ice_volume_init,yrs_continued,ice_volume_PD,'LineWidth',2);
grid on
xlabel('time (yrs)');
ylabel('ice volume');