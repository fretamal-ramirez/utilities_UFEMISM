% code to load and read variables from UFEMISM simulations
clear all
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));

% load the coastline from MEaSUREs
coast_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Coastline_Antarctica_v02.shp');
%function from Antarctic Maping Tools to project from x,y to lon,lat
[coast_lat,coast_lon]=ps2ll(coast_MEaSUREs.X,coast_MEaSUREs.Y);
%u_MEaSUREs=ncread('/Users/frre9931/Documents/PhD/MEaSUREs/antarctica_ice_velocity_450m_v2.nc','VX');
%v_MEaSUREs=ncread('/Users/frre9931/Documents/PhD/MEaSUREs/antarctica_ice_velocity_450m_v2.nc','VY');
uabs_MEaSUREs=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/MEaSUREs/Antarctica/surface_velocity_measures_2km.nc','uabs_surf');
x_MEaSUREs=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/MEaSUREs/Antarctica/surface_velocity_measures_2km.nc','x');
y_MEaSUREs=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/MEaSUREs/Antarctica/surface_velocity_measures_2km.nc','y');
% change -10000 to NaN for plots
for i=1:length(uabs_MEaSUREs(:,1))
    for j=1:length(uabs_MEaSUREs(1,:))
        if uabs_MEaSUREs(i,j)==-10000
            uabs_MEaSUREs(i,j)=NaN;
        end
    end
end

% UFEMISM output file to read
% ==========================

% path to the folder that has the outputs from UFEMISM simulation 
ufe_folder_path='/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/results_ant_PD_inversion_dHdt/';
% filename to load the main output in grid format
filename= [ufe_folder_path, 'main_output_ANT_grid.nc'];

Hi=ncread(filename,'Hi');
u_surf=ncread(filename,'u_surf');
v_surf=ncread(filename,'v_surf');
uabs_surf=ncread(filename,'uabs_surf');
Hb=ncread(filename,'Hb');
SMB=ncread(filename,'SMB');
BMB=ncread(filename,'BMB');
lon=ncread(filename,'lon');
lat=ncread(filename,'lat');
x=ncread(filename,'x');
y=ncread(filename,'y');
Hb_diff_first5kyr=Hb(:,:,end)-Hb(:,:,1);
Hi_diff_first5kyr=Hi(:,:,end)-Hi(:,:,1);
%% figure of Hb difference
figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon,lat,Hb_diff_first5kyr(:,:),'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
cbar=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar.Position(1) = cbar.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Bedrock elevation change after 10 kyr (m)');
t.Units='normalized';
t.Position(2)=1.05;
%print('Hb_diff10kyr','-dpng','-r300')

%% plot ice thickness
% after 10kyrs
figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon,lat,Hi_diff_first5kyr(:,:),20,'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','k');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Ice thickness change after 10 kyr (m)');
t.Units='normalized';
t.Position(2)=1.05;
colormap('jet');
clim([-1500 1500])
%print('Hi_diff10kyr','-dpng','-r300')

% initial simulation state
figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon,lat,Hi(:,:,1),20,'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Initial ice thickness (m)');
t.Units='normalized';
t.Position(2)=1.05;
clim([0 4000]);
%print('Hi_t0','-dpng','-r300')

% final simulation state
figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon,lat,Hi(:,:,end),20,'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Final ice thickness (m)');
t.Units='normalized';
t.Position(2)=1.05;
clim([0 4000]);
%print('Hi_tf','-dpng','-r300')
%% plot of surface ice velocity
% initial simulation state
figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon,lat,uabs_surf(:,:,1),20,'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Initial surface ice velocity (m/yr)');
t.Units='normalized';
t.Position(2)=1.05;
clim([0 1500]);
%print('uabs_t0','-dpng','-r300')

% final simulation state
figure('position',[100 100 500 500])
m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
m_contourf(lon,lat,uabs_surf(:,:,end),20,'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Final surface ice velocity (m/yr)');
t.Units='normalized';
t.Position(2)=1.05;
clim([0 1500]);
%print('uabs_tf','-dpng','-r300')

% observed from MEaSUREs
figure('position',[100 100 500 500])
%m_proj('stereographic','lat',-90,'long',0,'radius',37,'rectbox','on');
hold on
contourf(x_MEaSUREs,y_MEaSUREs,uabs_MEaSUREs(:,:)',20,'LineColor','none');
%m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
%m_line(coast_lon,coast_lat,'color','r');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Observed surface ice velocity MEaSUREs (m/yr)');
t.Units='normalized';
t.Position(2)=1.05;
clim([0 1500]);
%print('uabs_MEaSUREs','-dpng','-r300')
%% load the mesh

mesh_path_first= [ufe_folder_path, 'main_output_ANT_00001.nc']; %initial state
allow_mesh_update = false; % if remeshing is allowed in simulation
time_slice=1 ; % what time is going to be loaded, if load all of them is just NaNs idk why
allow_plot_mesh = true; % if we want to plot the mesh

% read_mesh_from_file + CL, GL and CF
mesh_first=read_mesh_from_file(mesh_path_first);
CL=ncread(mesh_path_first,'coastline',[1,1,1], [11788,2,time_slice]);
GL=ncread(mesh_path_first,'grounding_line',[1,1,1], [11788,2,time_slice]);
CF=ncread(mesh_path_first,'calving_front',[1,1,1], [11788,2,time_slice]);

if allow_mesh_update
    % path to the new mesh, for now just changing the number here
    mesh_path_update = [ufe_folder_path, 'main_output_ANT_00002.nc']; 
    mesh_updated=read_mesh_from_file(mesh_path_update);
    CL_mesh2=ncread(mesh_path_update,'coastline',[1,1,1], [11788,2,time_slice]);
    GL_mesh2=ncread(mesh_path_update,'grounding_line',[1,1,1], [11788,2,time_slice]);
    CF_mesh2=ncread(mesh_path_update,'calving_front',[1,1,1], [11788,2,time_slice]);
end
if allow_plot_mesh
    plot_mesh(mesh_first);
    hold on
    plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
    hold off
elseif allow_plot_mesh & allow_mesh_update
    plot_mesh(mesh)
    hold on
    plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
    plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
    plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
    hold off
end
%% try to add contour lines from the mesh to my plots

[Hi_fix, maskHi0]= Hi0_to_NaN(Hi);

figure('position',[100 100 500 500])
hold on
contourf(x,y,Hi_fix(:,:,end)',20,'LineColor','none');

cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Ice thickness (m)');
t.Units='normalized';
t.Position(2)=1.05;
%clim([0.001 1500]);
%colormap('jet');

plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');

%% now do same plot with velocities adding the maskHi0
figure('position',[100 100 500 500])
hold on
contourf(x,y,(uabs_surf(:,:,end).*maskHi0(:,:,end))',20,'LineColor','none');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Final surface ice velocity (m/yr)');
t.Units='normalized';
t.Position(2)=1.05;
clim([0 1500]);
%print('uabs_tf','-dpng','-r300')

plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');

%% plot for Antarctic Peninsula - Weddell Sea

figure()
m_proj('stereographic','lat',-69,'long',-60,'radius',12,'rectbox','on');
%m_coast('patch',[.7 .7 .7],'edgecolor','k');
hold on
m_contourf(lon,lat,Hi(:,:,1),'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-75 -70 -65 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
%m_grid;

%% plot AIS

% load data
filename_ant='/Users/frre9931/Desktop/UFEMISM2.0-main_old_Larsen_modifications/results_ant_Larsen10km/main_output_ANT_grid.nc';

Hi_ant=ncread(filename_ant,'Hi');
u_surf_ant=ncread(filename_ant,'u_surf');
v_surf_ant=ncread(filename_ant,'v_surf');
SMB_ant=ncread(filename_ant,'SMB');
BMB_ant=ncread(filename_ant,'BMB');
lon_ant=ncread(filename_ant,'lon');
lat_ant=ncread(filename_ant,'lat');

figure()
m_proj('stereographic','lat',-90,'long',0,'radius',30,'rectbox','on');
%m_coast('patch',[.7 .7 .7],'edgecolor','k');
hold on
m_contourf(lon_ant,lat_ant,Hi_ant(:,:,1),'LineColor','none');
m_grid('xtick',12,'XaxisLocation','top','ytick',[-80 -70 -60],'linest','--','box','off');
m_line(coast_lon,coast_lat,'color','r');
%m_grid;
%%
antmap('grid','lats',[-60 -80])
contourf(x_MEaSUREs,y_MEaSUREs,uabs_MEaSUREs(:,:),20,'LineColor','none');

%contourf(lat,lon,Hi(:,:,1))

figure('position',[100 100 1000 400])

subplot(1,3,1)
antmap('grid')              % initializes left map and plots grid
patchm(lat,long,'y')        % plots yellow continent atop  grid

subplot(1,3,2)  
antmap                      % initializes center map
patchm(lat,long,'b')        % plots blue continent
antmap('grid','color','m')  % overlays magenta grid

subplot(1,3,3) 
antmap                      % initializes right map
patchm(lat,long,[.5 .5 .5]) % plot gray continent
antmap('lats',-80:10:-50,...% plots lines at 50, 60, 70, & 80 S
    'lons',0:45:180,...     % and 0, 45, 95, etc longitude
    'color','red',...       % in red
    'linewidth',2,...       % kinda thick
    'linestyle',':')        % and dotted. 