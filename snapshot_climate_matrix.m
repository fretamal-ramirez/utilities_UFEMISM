%% script to load information from PMIP experiments to force the climate matrix

% the input data from snapshots for the climate matrix is T2m, Precip,
% orography (Hs), uas and vas. We need to load the data, and save it in one file
% for LGM and other one for piControl, UFEMISM will read the variables 
% Precipitation field needs to be in meter per month and T2m in Kelvin
% Precipitation from PMIP need to be changed as it comes in mm/month
clear all
model='MIROC-ES2L' ; % could be AWIESM1, MPI-ESM1-2, INM-CM4-8 or MIROC-ES2L  

switch model
    case 'AWIESM1'
        folder_path='/Users/frre9931/Desktop/PMIP-data/AWI-ESM-1-1-LR/';
    case 'MPI-ESM1-2'
        folder_path='/Users/frre9931/Desktop/PMIP-data/MPI-ESM1-2/';
    case 'INM-CM4-8'
        folder_path='/Users/frre9931/Desktop/PMIP-data/INM-CM4-8/';        
    case 'MIROC-ES2L'
        folder_path='/Users/frre9931/Desktop/PMIP-data/MIROC-ES2L/'; 
    otherwise
        warning('unexpected PMIP model')
end

T2m_lgm=ncread([folder_path,'PMIP4_tas_Amon_',model,'_lgm_monClim.nc'],'tas');
T2m_piControl=ncread([folder_path,'PMIP4_tas_Amon_',model,'_piControl_monClim.nc'],'tas');
Precip_lgm=ncread([folder_path,'PMIP4_pr_Amon_',model,'_lgm_monClim.nc'],'pr')/1000;
Precip_piControl=ncread([folder_path,'PMIP4_pr_Amon_',model,'_piControl_monClim.nc'],'pr')/1000;
Hs_lgm=ncread([folder_path,'PMIP4_orog_fx_',model,'_lgm.nc'],'orog');
Hs_piControl=ncread([folder_path,'PMIP4_orog_fx_',model,'_piControl.nc'],'orog');
uas_lgm=ncread([folder_path,'PMIP4_uas_Amon_',model,'_lgm_monClim.nc'],'uas');
uas_piControl=ncread([folder_path,'PMIP4_uas_Amon_',model,'_piControl_monClim.nc'],'uas');
vas_lgm=ncread([folder_path,'PMIP4_vas_Amon_',model,'_lgm_monClim.nc'],'vas');
vas_piControl=ncread([folder_path,'PMIP4_vas_Amon_',model,'_piControl_monClim.nc'],'vas');
lon=ncread([folder_path,'PMIP4_tas_Amon_',model,'_lgm_monClim.nc'],'lon');
lat=ncread([folder_path,'PMIP4_tas_Amon_',model,'_lgm_monClim.nc'],'lat');
time=ncread([folder_path,'PMIP4_tas_Amon_',model,'_lgm_monClim.nc'],'time');
month=[1:12]';

% Check if the data has a regular o irregular grid in lat
diff_lat=lat( 2) - lat( 1); % if is regular this value should be similar
irregular=false;
for i = 3: length(lat)
        delta = lat( i) - lat( i-1);
        if (abs( 1 - delta / diff_lat) > 1e-5) 
            irregular=true;
            disp('Latitude coordinate is irregular')
            break
        end
end

% If is irregular interpolate to a regular grid and replace previous values of lat
% create the mesh first
[lon_grid,lat_grid]=meshgrid(lon,lat);

if irregular

% Interpolate to a regular grid
% ===================

lat_interp_regular=linspace(lat(1),lat(end),length(lat))';

[out_lon_grid,out_lat_grid]=meshgrid(lon,lat_interp_regular);
Hs_lgm=interp2(lon_grid,lat_grid,Hs_lgm',out_lon_grid,out_lat_grid,'linear')';
Hs_piControl=interp2(lon_grid,lat_grid,Hs_piControl',out_lon_grid,out_lat_grid,'linear')';
% variables with time
    for m=1:length(month)
        T2m_lgm(:,:,m)          = interp2(lon_grid,lat_grid,T2m_lgm(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        T2m_piControl(:,:,m)    = interp2(lon_grid,lat_grid,T2m_piControl(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        Precip_lgm(:,:,m)       = interp2(lon_grid,lat_grid,Precip_lgm(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        Precip_piControl(:,:,m) = interp2(lon_grid,lat_grid,Precip_piControl(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        uas_lgm(:,:,m)          = interp2(lon_grid,lat_grid,uas_lgm(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        uas_piControl(:,:,m)    = interp2(lon_grid,lat_grid,uas_piControl(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        vas_lgm(:,:,m)          = interp2(lon_grid,lat_grid,vas_lgm(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
        vas_piControl(:,:,m)    = interp2(lon_grid,lat_grid,vas_piControl(:,:,m)',out_lon_grid,out_lat_grid,'linear')';
    end

lat=lat_interp_regular;
disp('Latitude coordinate has been interpolated to a regular grid')
else % regular, just to make it work in the following code..
    out_lon_grid=lon_grid;
    out_lat_grid=lat_grid;
end

% check if Precip is lower than 0
for m=1:length(month)
    for i=1:length(Precip_lgm(:,1,1))
        for j=1:length(Precip_lgm(1,:,1))
            if Precip_lgm(i,j,m) < 0
                disp('LGM Negative precipitation, making value equal 0')
                Precip_lgm(i,j,m)=0;
            elseif Precip_piControl(i,j,m) < 0
                disp('piControl Negative precipitation, making value equal 0')
                Precip_piControl(i,j,m)=0;
            end
        end
    end
end
%% reproject data to UFEMISM mesh (x,y) instead of (lon,lat)
% add functions that I copied from Ufemism to reproject
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0-main_old_Larsen_modifications/tools/matlab'));
x=zeros(size(lon_grid)); y=zeros(size(lat_grid));
for i=1:length(lat)
    for j=1:length(lon)
        if out_lon_grid(i,j)<0
            aux_out=out_lon_grid(i,j)+360;
        else
            aux_out=out_lon_grid(i,j);
        end
        [x_out,y_out]=oblique_sg_projection(aux_out,out_lat_grid(i,j));
        x(i,j)=x_out;
        y(i,j)=y_out;
    end
end
lon_ufe=ncread('/Users/frre9931/Documents/PhD/RACMO/RACMO_clim_1979_2021.nc','lon');
lat_ufe=ncread('/Users/frre9931/Documents/PhD/RACMO/RACMO_clim_1979_2021.nc','lat');
x_ufe=ncread('/Users/frre9931/Documents/PhD/RACMO/RACMO_clim_1979_2021.nc','x');
y_ufe=ncread('/Users/frre9931/Documents/PhD/RACMO/RACMO_clim_1979_2021.nc','y');
[x_ufegrid,y_ufegrid]=meshgrid(x_ufe,y_ufe);
Hs_lgm_ufegrid       = griddata(x,y,squeeze(Hs_lgm(:,:))',x_ufegrid,y_ufegrid,'linear')';
Hs_piControl_ufegrid = griddata(x,y,squeeze(Hs_piControl(:,:))',x_ufegrid,y_ufegrid,'linear')';
for m=1:length(T2m_lgm(1,1,:))
    T2m_lgm_ufegrid(:,:,m)          = griddata(x,y,squeeze(T2m_lgm(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    T2m_piControl_ufegrid(:,:,m)    = griddata(x,y,squeeze(T2m_piControl(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    Precip_lgm_ufegrid(:,:,m)       = griddata(x,y,squeeze(Precip_lgm(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    Precip_piControl_ufegrid(:,:,m) = griddata(x,y,squeeze(Precip_piControl(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    uas_lgm_ufegrid(:,:,m)          = griddata(x,y,squeeze(uas_lgm(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    uas_piControl_ufegrid(:,:,m)    = griddata(x,y,squeeze(uas_piControl(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    vas_lgm_ufegrid(:,:,m)          = griddata(x,y,squeeze(vas_lgm(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
    vas_piControl_ufegrid(:,:,m)    = griddata(x,y,squeeze(vas_piControl(:,:,m))',x_ufegrid,y_ufegrid,'linear')';
end
%% plot to check if its fine
path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));

figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon_ufe,lat_ufe,Precip_piControl_ufegrid(:,:,1),20,'LineColor','none')
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('T2m');
t.Units='normalized';
t.Position(2)=1.05;
%% create the netCDF output for each snapshot
% Create empty file
% =================

% create netcdf file following UFEMISM
% improve the naming according to the input!!
ncid = netcdf.create(['/Users/frre9931/Documents/PhD/ANT_UFEMISM/preprocessing_input/snapshot_lgm_',model,'.nc'],'CLOBBER');

% Define dimensions
% =================

dim_lon = netcdf.defDim(ncid,'lon',size(lon,1));
dim_lat = netcdf.defDim(ncid,'lat',size(lat,1));
dim_month = netcdf.defDim(ncid,'month',size(month,1));

% Define variables
% ================

id_lon = netcdf.defVar(ncid,'lon','double',dim_lon);
id_lat = netcdf.defVar(ncid,'lat','double',dim_lat);
id_month = netcdf.defVar(ncid,'month','double',dim_month);
id_Hs  = netcdf.defVar(ncid,'Hs','double',[dim_lon, dim_lat]);
id_T2m = netcdf.defVar(ncid,'T2m','double',[dim_lon,dim_lat,dim_month]);
id_Precip = netcdf.defVar(ncid,'Precip','double',[dim_lon,dim_lat,dim_month]);
id_uas = netcdf.defVar(ncid,'uas','double',[dim_lon,dim_lat,dim_month]);
id_vas = netcdf.defVar(ncid,'vas','double',[dim_lon,dim_lat,dim_month]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_lon,'standard_name','longitude'); 
netcdf.putAtt(ncid,id_lat,'standard_name','latitude');
netcdf.putAtt(ncid,id_month,'standard_name','month');

netcdf.putAtt(ncid,id_T2m,'standard_name','Temperature 2 m above surface'); 
netcdf.putAtt(ncid,id_T2m,'units','K');

netcdf.putAtt(ncid,id_Precip,'standard_name','Precipitation'); 
netcdf.putAtt(ncid,id_Precip,'units','mm per month'); % do I need to change it? check!

netcdf.putAtt(ncid,id_Hs,'standard_name','Surface topography'); 
netcdf.putAtt(ncid,id_Hs,'units','meters');

netcdf.putAtt(ncid,id_uas,'standard_name','Eastward near-surface wind'); 
netcdf.putAtt(ncid,id_uas,'units','meters per second'); % should be m per month?

netcdf.putAtt(ncid,id_vas,'standard_name','Northward near-surface wind'); 
netcdf.putAtt(ncid,id_vas,'units','meters per second');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_lon,lon);
netcdf.putVar(ncid,id_lat,lat);
netcdf.putVar(ncid,id_month,month);

netcdf.putVar(ncid,id_T2m,T2m_lgm);
netcdf.putVar(ncid,id_Precip,Precip_lgm);
netcdf.putVar(ncid,id_Hs,Hs_lgm);
netcdf.putVar(ncid,id_uas,uas_lgm);
netcdf.putVar(ncid,id_vas,vas_lgm);

% Close file
% ==========

netcdf.close(ncid);


%% now repeat the same to create the files for piControl

% create the netcdf file
ncid = netcdf.create(['/Users/frre9931/Documents/PhD/ANT_UFEMISM/preprocessing_input/snapshot_piControl_',model,'.nc'],'CLOBBER');

% Define dimensions
% =================

dim_lon = netcdf.defDim(ncid,'lon',size(lon,1));
dim_lat = netcdf.defDim(ncid,'lat',size(lat,1));
dim_month = netcdf.defDim(ncid,'month',size(month,1));

% Define variables
% ================

id_lon = netcdf.defVar(ncid,'lon','double',dim_lon);
id_lat = netcdf.defVar(ncid,'lat','double',dim_lat);
id_month = netcdf.defVar(ncid,'month','double',dim_month);
id_Hs  = netcdf.defVar(ncid,'Hs','double',[dim_lon, dim_lat]);
id_T2m = netcdf.defVar(ncid,'T2m','double',[dim_lon,dim_lat,dim_month]);
id_Precip = netcdf.defVar(ncid,'Precip','double',[dim_lon,dim_lat,dim_month]);
id_uas = netcdf.defVar(ncid,'uas','double',[dim_lon,dim_lat,dim_month]);
id_vas = netcdf.defVar(ncid,'vas','double',[dim_lon,dim_lat,dim_month]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_lon,'standard_name','longitude'); 
netcdf.putAtt(ncid,id_lat,'standard_name','latitude');
netcdf.putAtt(ncid,id_month,'standard_name','month');

netcdf.putAtt(ncid,id_T2m,'standard_name','Temperature 2 m above surface'); 
netcdf.putAtt(ncid,id_T2m,'units','K');

netcdf.putAtt(ncid,id_Precip,'standard_name','Precipitation'); 
netcdf.putAtt(ncid,id_Precip,'units','mm per month'); % do I need to change it? check!

netcdf.putAtt(ncid,id_Hs,'standard_name','Surface topography'); 
netcdf.putAtt(ncid,id_Hs,'units','meters');

netcdf.putAtt(ncid,id_uas,'standard_name','Eastward near-surface wind'); 
netcdf.putAtt(ncid,id_uas,'units','meters per second'); % should be m per month?

netcdf.putAtt(ncid,id_vas,'standard_name','Northward near-surface wind'); 
netcdf.putAtt(ncid,id_vas,'units','meters per second');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_lon,lon);
netcdf.putVar(ncid,id_lat,lat);
netcdf.putVar(ncid,id_month,month);

netcdf.putVar(ncid,id_T2m,T2m_piControl);
netcdf.putVar(ncid,id_Precip,Precip_piControl);
netcdf.putVar(ncid,id_Hs,Hs_piControl);
netcdf.putVar(ncid,id_uas,uas_piControl);
netcdf.putVar(ncid,id_vas,vas_piControl);

% Close file
% ==========

netcdf.close(ncid);

%% save netcdf file with ANT projection (x,y)

% Example variables:
filename = ['/Users/frre9931/Documents/PhD/ANT_UFEMISM/preprocessing_input/snapshot_ANT_lgm_',model,'.nc'];

% Create NetCDF file
nccreate(filename, 'month', 'Dimension', {'month', length(month)});
nccreate(filename, 'x', 'Dimensions',{'x', length(x_ufe)});
nccreate(filename, 'y', 'Dimensions',{'y', length(y_ufe)});
nccreate(filename, 'uas', 'Dimensions', {'x', size(uas_lgm_ufegrid, 1), 'y', size(uas_lgm_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'vas', 'Dimensions', {'x', size(vas_lgm_ufegrid, 1), 'y', size(vas_lgm_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'T2m', 'Dimensions', {'x', size(T2m_lgm_ufegrid, 1), 'y', size(T2m_lgm_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'Precip', 'Dimensions', {'x', size(Precip_lgm_ufegrid, 1), 'y', size(Precip_lgm_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'Hs', 'Dimensions', {'x', size(Hs_lgm_ufegrid, 1), 'y', size(Hs_lgm_ufegrid, 2)});
nccreate(filename, 'lat', 'Dimensions', {'x', size(lat_ufe, 1), 'y', size(lat_ufe, 2)});
nccreate(filename, 'lon', 'Dimensions', {'x', size(lon_ufe, 1), 'y', size(lon_ufe, 2)});

% Write data
ncwrite(filename, 'uas', uas_lgm_ufegrid);
ncwrite(filename, 'vas', vas_lgm_ufegrid);
ncwrite(filename, 'T2m', T2m_lgm_ufegrid);
ncwrite(filename, 'Precip', Precip_lgm_ufegrid);
ncwrite(filename, 'Hs', Hs_lgm_ufegrid);
ncwrite(filename, 'lat', lat_ufe);
ncwrite(filename, 'lon', lon_ufe);
ncwrite(filename,'month', month);
ncwrite(filename,'x', x_ufe);
ncwrite(filename,'y', y_ufe);

% Add metadata
ncwriteatt(filename, 'lat', 'long_name', 'latitude');
ncwriteatt(filename, 'lon', 'long_name', 'longitude');
ncwriteatt(filename, 'uas', 'units', 'm/s');
ncwriteatt(filename, 'vas', 'units', 'm/s');
ncwriteatt(filename, 'T2m', 'units', 'K');
ncwriteatt(filename, 'Precip', 'units', 'meters (total accumulation over that month :) )');
ncwriteatt(filename, 'Hs', 'units', 'm');
ncwriteatt(filename,'month', 'units', 'months')
ncwriteatt(filename,'x','units','meters');
ncwriteatt(filename,'y','units','meters');

% Example variables:
filename = ['/Users/frre9931/Documents/PhD/ANT_UFEMISM/preprocessing_input/snapshot_ANT_piControl_',model,'.nc'];

% Create NetCDF file
nccreate(filename, 'month', 'Dimension', {'month', length(month)});
nccreate(filename, 'x', 'Dimensions',{'x', length(x_ufe)});
nccreate(filename, 'y', 'Dimensions',{'y', length(y_ufe)});
nccreate(filename, 'uas', 'Dimensions', {'x', size(uas_piControl_ufegrid, 1), 'y', size(uas_piControl_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'vas', 'Dimensions', {'x', size(vas_piControl_ufegrid, 1), 'y', size(vas_piControl_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'T2m', 'Dimensions', {'x', size(T2m_piControl_ufegrid, 1), 'y', size(T2m_piControl_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'Precip', 'Dimensions', {'x', size(Precip_piControl_ufegrid, 1), 'y', size(Precip_piControl_ufegrid, 2), 'month', length(month)});
nccreate(filename, 'Hs', 'Dimensions', {'x', size(Hs_piControl_ufegrid, 1), 'y', size(Hs_piControl_ufegrid, 2)});
nccreate(filename, 'lat', 'Dimensions', {'x', size(lat_ufe, 1), 'y', size(lat_ufe, 2)});
nccreate(filename, 'lon', 'Dimensions', {'x', size(lon_ufe, 1), 'y', size(lon_ufe, 2)});

% Write data
ncwrite(filename, 'uas', uas_piControl_ufegrid);
ncwrite(filename, 'vas', vas_piControl_ufegrid);
ncwrite(filename, 'T2m', T2m_piControl_ufegrid);
ncwrite(filename, 'Precip', Precip_piControl_ufegrid);
ncwrite(filename, 'Hs', Hs_piControl_ufegrid);
ncwrite(filename, 'lat', lat_ufe);
ncwrite(filename, 'lon', lon_ufe);
ncwrite(filename,'month', month);
ncwrite(filename,'x', x_ufe);
ncwrite(filename,'y', y_ufe);

% Add metadata
ncwriteatt(filename, 'lat', 'long_name', 'latitude');
ncwriteatt(filename, 'lon', 'long_name', 'longitude');
ncwriteatt(filename, 'uas', 'units', 'm/s');
ncwriteatt(filename, 'vas', 'units', 'm/s');
ncwriteatt(filename, 'T2m', 'units', 'K');
ncwriteatt(filename, 'Precip', 'units', 'meters (total accumulation over that month :) )');
ncwriteatt(filename, 'Hs', 'units', 'm');
ncwriteatt(filename,'month', 'units', 'months')
ncwriteatt(filename,'x','units','meters');
ncwriteatt(filename,'y','units','meters');