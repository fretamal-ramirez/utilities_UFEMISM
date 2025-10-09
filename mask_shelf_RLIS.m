% code to create a mask of 1s to whatever is outside the continental ice
% inside the ROI

ice_shelf=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/shape_outside_icesheet.shp');
% use the bounding box to create a square grid of x and y
grid_res=1000; % 5 km
x=ice_shelf.BoundingBox(1,1):grid_res:ice_shelf.BoundingBox(2,1);
y=ice_shelf.BoundingBox(1,2):grid_res:ice_shelf.BoundingBox(2,2);
[x_grid,y_grid]=meshgrid(x,y);

[in,on]=inpolygon(x_grid,y_grid,ice_shelf.X,ice_shelf.Y);

figure()
contourf(x_grid,y_grid,in)

% is looking good, create an actual data with 0s to fill with ones
mask_shelf=zeros(size(x_grid));
mask_shelf(in)=1;

% create a time vector
faketime=[1980:10:2500];

mask_shelf_time=zeros([size(mask_shelf),length(faketime)]);
for k=1:length(faketime)
    mask_shelf_time(:,:,k)=mask_shelf(:,:);
end

%% save the dataset as a file with x,y

% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_iceshelf_calving.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,2));
dim_y = netcdf.defDim(ncid,'y',size(y,2));
dim_time = netcdf.defDim(ncid,'time',size(faketime,2));

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

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);
netcdf.putVar(ncid,id_time,faketime);

netcdf.putVar(ncid,id_mask,permute(mask_shelf_time,[2,1,3]));

% Close file
% ==========

netcdf.close(ncid);