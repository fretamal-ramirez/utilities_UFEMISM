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

%% code to create a mask where is PD ocean according to MEaSUREs
ice_shape=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Coastline_Antarctica_v02_simplified2km.shp');

% use the bounding box to create a square grid of x and y
grid_res=5000; % km

%boundingBox=[-2800225,-2800025;
%    2800025,2800225];
boundingBox=[-3040000,-3040000;
    3040000,3040000];

x=boundingBox(1,1):grid_res:boundingBox(2,1);
y=boundingBox(1,2):grid_res:boundingBox(2,2);
[x_grid,y_grid]=meshgrid(x,y);

[in,on]=inpolygon(x_grid,y_grid,ice_shape.X,ice_shape.Y);

% is looking good, create an actual data with 0s to fill with ones
mask_shelf=ones(size(x_grid));
mask_shelf(in)=0;

figure()
contourf(x_grid,y_grid,mask_shelf)

%% save the dataset as a file with x,y

% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving1.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,2));
dim_y = netcdf.defDim(ncid,'y',size(y,2));
%dim_time = netcdf.defDim(ncid,'time',size(faketime,2));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
%id_time = netcdf.defVar(ncid,'time','double',dim_time);

id_mask  = netcdf.defVar(ncid,'mask','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 
%netcdf.putAtt(ncid,id_time,'standard_name','time');
%netcdf.putAtt(ncid,id_time,'units','years');

netcdf.putAtt(ncid,id_mask,'standard_name','mask'); 
netcdf.putAtt(ncid,id_mask,'units','unitless');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);
%netcdf.putVar(ncid,id_time,faketime);

netcdf.putVar(ncid,id_mask,permute(mask_shelf,[2,1]));

% Close file
% ==========

netcdf.close(ncid);

%% something similar but using the mask from Bedmachine
bedmachine_mask=ncread('/Users/frre9931/Documents/PhD/BedMachine/BedMachineAntarctica-v3.nc','mask');
bedmachine_x=ncread('/Users/frre9931/Documents/PhD/BedMachine/BedMachineAntarctica-v3.nc','x');
bedmachine_y=ncread('/Users/frre9931/Documents/PhD/BedMachine/BedMachineAntarctica-v3.nc','y');

% use the bounding box to create a square grid of x and y
grid_res=500; % m

boundingBox=[-3040000,-3040000;
    3040000,3040000];

x=boundingBox(1,1):grid_res:boundingBox(2,1);
y=boundingBox(1,2):grid_res:boundingBox(2,2);
[x_grid,y_grid]=meshgrid(x,y);

mask_outgrid=interp2(double(bedmachine_x),double(bedmachine_y),double(bedmachine_mask)',x_grid,y_grid);
figure()
contourf(x_grid,y_grid,mask_outgrid)

mask_filled = imfill(mask_outgrid > 0, 'holes'); % fill holes in the binary mask
mask_final = ones(size(mask_filled)).*~mask_filled;

figure()
contourf(x_grid,y_grid,mask_final)

%% save the dataset as a file with x,y

% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving_bedmachine.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,2));
dim_y = netcdf.defDim(ncid,'y',size(y,2));
%dim_time = netcdf.defDim(ncid,'time',size(faketime,2));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
%id_time = netcdf.defVar(ncid,'time','double',dim_time);

id_mask  = netcdf.defVar(ncid,'mask','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 
%netcdf.putAtt(ncid,id_time,'standard_name','time');
%netcdf.putAtt(ncid,id_time,'units','years');

netcdf.putAtt(ncid,id_mask,'standard_name','mask'); 
netcdf.putAtt(ncid,id_mask,'units','unitless');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);
%netcdf.putVar(ncid,id_time,faketime);

netcdf.putVar(ncid,id_mask,permute(mask_final,[2,1]));

% Close file
% ==========

netcdf.close(ncid);

%% do something similar but using data from Bedmachine in 2km from Ufe2.0 to
% create a mask at 2km.
path_bedmachine='/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc';
%path_bedmachine='/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/external/data/topography/Antarctica/BedMachineAntarctica_v3_2km.nc';

bedmachine_hi=ncread(path_bedmachine,'Hi');
bedmachine_x=ncread(path_bedmachine,'x');
bedmachine_y=ncread(path_bedmachine,'y');

% use the bounding box to create a square grid of x and y
grid_res=5000; % m

boundingBox=[-3040000,-3040000;
    3040000,3040000];

x=boundingBox(1,1):grid_res:boundingBox(2,1);
y=boundingBox(1,2):grid_res:boundingBox(2,2);
[x_grid,y_grid]=meshgrid(x,y);

bedmachine_mask=zeros(size(bedmachine_hi));
for i= 1:size(bedmachine_hi,1)
    for j=1:size(bedmachine_hi,2)
        if bedmachine_hi(i,j)>0
            bedmachine_mask(i,j)=1;
        end
    end
end

mask_outgrid=interp2(double(bedmachine_x),double(bedmachine_y),double(bedmachine_mask)',x_grid,y_grid);
figure()
contourf(x_grid,y_grid,mask_outgrid)

mask_filled = imfill(mask_outgrid > 0, 'holes'); % fill holes in the binary mask
mask_final = ones(size(mask_filled)).*~mask_filled;

%mask_final(905,471:473)=0; % make 0 one of the points giving problems while running
%mask_final(906:907,470:471)=0;
figure()
contourf(x_grid,y_grid,mask_final)

%%
% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving_bedmachine_5km.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,2));
dim_y = netcdf.defDim(ncid,'y',size(y,2));
%dim_time = netcdf.defDim(ncid,'time',size(faketime,2));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
%id_time = netcdf.defVar(ncid,'time','double',dim_time);

id_mask  = netcdf.defVar(ncid,'mask','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 
%netcdf.putAtt(ncid,id_time,'standard_name','time');
%netcdf.putAtt(ncid,id_time,'units','years');

netcdf.putAtt(ncid,id_mask,'standard_name','mask'); 
netcdf.putAtt(ncid,id_mask,'units','unitless');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);
%netcdf.putVar(ncid,id_time,faketime);

netcdf.putVar(ncid,id_mask,permute(mask_final,[2,1]));

% Close file
% ==========

netcdf.close(ncid);

%% do something similar but using data from Bedmachine in 5km from Ufe2.0 to
% create a mask at 5km.
bedmachine_hi=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/external/data/topography/Antarctica/BedMachineAntarctica_v3_2km.nc','Hi');
bedmachine_x=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/external/data/topography/Antarctica/BedMachineAntarctica_v3_2km.nc','x');
bedmachine_y=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/external/data/topography/Antarctica/BedMachineAntarctica_v3_2km.nc','y');

% use the bounding box to create a square grid of x and y
grid_res=2000; % km

boundingBox=[-3040000,-3040000;
    3040000,3040000];

x=boundingBox(1,1):grid_res:boundingBox(2,1);
y=boundingBox(1,2):grid_res:boundingBox(2,2);
[x_grid,y_grid]=meshgrid(x,y);

bedmachine_mask=zeros(size(bedmachine_hi));
for i= 1:size(bedmachine_hi,1)
    for j=1:size(bedmachine_hi,2)
        if bedmachine_hi(i,j)>0
            bedmachine_mask(i,j)=1;
        end
    end
end

mask_outgrid=interp2(double(bedmachine_x),double(bedmachine_y),double(bedmachine_mask)',x_grid,y_grid);
figure()
contourf(x_grid,y_grid,mask_outgrid)

mask_filled = imfill(mask_outgrid > 0, 'holes'); % fill holes in the binary mask
mask_final = ones(size(mask_filled)).*~mask_filled;

figure()
contourf(x_grid,y_grid,mask_final)

%%
% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/mask_ocean_calving_bedmachine_2km.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,2));
dim_y = netcdf.defDim(ncid,'y',size(y,2));
%dim_time = netcdf.defDim(ncid,'time',size(faketime,2));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
%id_time = netcdf.defVar(ncid,'time','double',dim_time);

id_mask  = netcdf.defVar(ncid,'mask','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 
%netcdf.putAtt(ncid,id_time,'standard_name','time');
%netcdf.putAtt(ncid,id_time,'units','years');

netcdf.putAtt(ncid,id_mask,'standard_name','mask'); 
netcdf.putAtt(ncid,id_mask,'units','unitless');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);
%netcdf.putVar(ncid,id_time,faketime);

netcdf.putVar(ncid,id_mask,permute(mask_final,[2,1]));

% Close file
% ==========

netcdf.close(ncid);