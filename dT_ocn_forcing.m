%create a transient dT_ocean forcing to use as input in UFEMISM2.0
clear all
% let's start from 1990 to 2020
time = 2000:50:2500;
% dT to perturb the ocean field
dT = 1.0; 
% create a vector with the values
dT_constant = ones(size(time)).*dT;
dT_ramp = linspace(0,dT,length(time));

% possible names for variable in file 'dT||dT_ocean||dTo'
% filename with path
filename = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/ocean_forcings/dT_const_1e1.nc';

% save as netcdf file
nccreate(filename, 'time', 'Dimension', {'time', length(time)});
nccreate(filename, 'dT_ocean', 'Dimensions',{'time', length(dT_constant)});

% Write data
ncwrite(filename, 'dT_ocean', dT_constant);
ncwrite(filename, 'time', time);

% Add metadata
ncwriteatt(filename,'time','units','years');
ncwriteatt(filename,'dT_ocean','units','degrees');