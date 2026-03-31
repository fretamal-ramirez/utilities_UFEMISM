% Code to implement gradual retreat masks. The code starts with the mask of
% PD created from Bedmachine v2, where 1 are assigned to open ocean and 0
% where is ice. In the code in UFEMISM all the values with 1 are deleted

% The idea of this code is to start from there and use the provided
% polygons from QGIS to change these values from 0 to 1

% Is it better to also keep the ocean as 1? It should not be needed as the
% code does not allow ice to grow ice further than PD in the forcing
% simulations. So if I do a file just with 1s in the mask should also work

% read mask from present-day
PDmask = ncread('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving_bedmachine_2km.nc','mask')';
%PDmask = permute(PDmask,[2,1]);
x_PDmask = ncread('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving_bedmachine_2km.nc','x');
y_PDmask = ncread('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving_bedmachine_2km.nc','y');
[x_PDgrid,y_PDgrid]=meshgrid(x_PDmask,y_PDmask);

buffers_name = { ...
    'CF_buffer25km', ...
    'CF_buffer50km', ...
    'CF_buffer75km', ...
    'CF_buffer100km', ...
    'CF_buffer125km', ...
    'CF_buffer150km', ...
    'CF_buffer175km', ...
    'CF_buffer200km', ...
};

time_init = 2000;

% in which years each buffer will be applied
%years_for_retreat = [2100,2200,2300,2400] ;
years_for_retreat = [2050,2100,2150,2200,2250,2300,2350,2400] ;

% create a vector with the time init and years of retreat masks
time_slides = [time_init, years_for_retreat];

datapath = '/Users/frre9931/Documents/PhD/RiiserLarsen/';

masks_stacked = PDmask;

% loop trough every buffer
for i=1:length(buffers_name)

filepath = [datapath, buffers_name{i},'.shp'];

% read buffered data
buffer = shaperead(filepath);

% find points inside the mask
[in,on]=inpolygon(x_PDgrid,y_PDgrid,buffer.X,buffer.Y);

% copy PDmask to new variable and assign 1 value
mask_retreat = PDmask;
mask_retreat(in)=1;

% add each mask to the 3D file
masks_stacked(:,:,i+1) = mask_retreat;

end

% % check if the masks are well implemented
% figure()
% contourf(x_PDgrid,y_PDgrid,masks_stacked(:,:,2))

%% save the masks

% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_gradual_retreat.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x_PDmask,1));
dim_y = netcdf.defDim(ncid,'y',size(y_PDmask,1));
dim_time = netcdf.defDim(ncid,'time',size(time_slides,2));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
id_time = netcdf.defVar(ncid,'time','double',dim_time);

id_mask  = netcdf.defVar(ncid,'mask','double',[dim_x, dim_y, dim_time]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 
netcdf.putAtt(ncid,id_time,'standard_name','time');
netcdf.putAtt(ncid,id_time,'units','years');

netcdf.putAtt(ncid,id_mask,'standard_name','mask'); 
netcdf.putAtt(ncid,id_mask,'units','unitless');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x_PDmask);
netcdf.putVar(ncid,id_y,y_PDmask);
netcdf.putVar(ncid,id_time,time_slides);

netcdf.putVar(ncid,id_mask,permute(masks_stacked,[2,1,3]));

% Close file
% ==========

netcdf.close(ncid);