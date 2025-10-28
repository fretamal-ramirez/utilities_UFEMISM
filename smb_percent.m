% code to decrease homogeneusly the SMB in a certain percentage

% first read input data
smb_input=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/RACMO/Antarctica/RACMO_SMB_1979_2016.nc','SMB');
x=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/RACMO/Antarctica/RACMO_SMB_1979_2016.nc','x');
y=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/RACMO/Antarctica/RACMO_SMB_1979_2016.nc','y');

percent=25;
smb_decreased=ones(size(smb_input)); %allocate variable
smb_decreased=smb_input*(percent/100);

%% save new data set to force UFEMISM

% create netcdf file
ncid = netcdf.create(['/Users/frre9931/Documents/PhD/RACMO/RACMO_decreased/RACMO_SMB_1979_2016_',char(string(percent)),'percent.nc'],'CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,1));
dim_y = netcdf.defDim(ncid,'y',size(y,1));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);

id_SMB  = netcdf.defVar(ncid,'SMB','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 

netcdf.putAtt(ncid,id_SMB,'standard_name','Surface mass balance'); 
netcdf.putAtt(ncid,id_SMB,'units','meters per year');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);

netcdf.putVar(ncid,id_SMB,smb_decreased);

% Close file
% ==========

netcdf.close(ncid);

%% apply the same to the basal roughness after inversion,

% this will have to be implemented in the mesh file rather than grid
% version
clear all
pathtofile='/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20_PMP_roughness_max20';
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));

mesh1 = read_mesh_from_file([pathtofile,'/main_output_ANT_00001.nc']);
till_angle=ncread([pathtofile,'/main_output_ANT_00001.nc'],'till_friction_angle');

percent=25;
till_decreased=ones(size(till_angle)); %allocate variable
till_decreased=till_angle*(percent/100);

%% create a copy of the file original file and modify it.
% this is done to keep all the "mesh" parameters.
path_output='/Users/frre9931/Documents/PhD/bed_roughness/modified/till_friction_angle_maxphi20_';
copyfile([pathtofile,'/main_output_ANT_00001.nc'],[path_output,char(string(percent)),'percent.nc'])
ncid = netcdf.open([path_output,char(string(percent)),'percent.nc'],"NC_WRITE");
varid = netcdf.inqVarID(ncid,"till_friction_angle");
netcdf.putVar(ncid,varid, till_decreased)
netcdf.close(ncid)

