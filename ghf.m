% ghf comparison in wDML
path = '/Users/frre9931/Desktop/UFEMISM2.0_porting/external/data/geothermal_heat_flux/';

% Shapiro 2004
shapiro.ghf = ncread([path,'ShapiroRitzwoller2004_global.nc'],'hflux');
shapiro.lon = ncread([path,'ShapiroRitzwoller2004_global.nc'],'Longitude');
shapiro.lat = ncread([path,'ShapiroRitzwoller2004_global.nc'],'Latitude');

% Staal topocorr
stal.ghf = ncread([path,'Antarctica/Staal_2020_topocorr_4km.nc'],'hflux');
stal.x = ncread([path,'Antarctica/Staal_2020_topocorr_4km.nc'],'x');
stal.y = ncread([path,'Antarctica/Staal_2020_topocorr_4km.nc'],'y');

%Loesing topocorr
loesing.ghf = ncread([path,'Antarctica/LoesingEbbing_2021_topocorr_4km.nc'],'hflux');
loesing.x = ncread([path,'Antarctica/LoesingEbbing_2021_topocorr_4km.nc'],'x');
loesing.y = ncread([path,'Antarctica/LoesingEbbing_2021_topocorr_4km.nc'],'y');

%% plot
figure()
tiledlayout(1,3)
nexttile
contourf(shapiro.lon,shapiro.lat,shapiro.ghf')
clim([0 0.7]);
title('Shapiro')

nexttile
contourf(stal.x,stal.y,stal.ghf')
clim([0 0.7]);
title('Stal')

nexttile
contourf(loesing.x,loesing.y,loesing.ghf')
colorbar
clim([0 0.7]);
title('Loesing')

% next step zoom in in R-LIS :)
