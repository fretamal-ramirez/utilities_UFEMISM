% code to load and read variables from UFEMISM simulations
clear all; clc;

%========= PATH TO OUTPUT UFEMISM DIRECTORY ==========
output_folder = 'results_ant_PD_inversion_dHdt_init_R-LIS_gamma10_HR';
%ufe_folder_path=['/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/', output_folder];
%ufe_folder_path=['/Users/frre9931/Desktop/UFEMISM2.0_porting/', output_folder];
ufe_folder_path=['/Users/frre9931/Desktop/tetralith_results/', output_folder];
allow_plot_mesh = true; % if we want to plot the mesh
allow_save_plots = false;
path_save = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/';
allow_mesh_update = false; % if remeshing is allowed in simulation
% number pointing the file with updated mesh
number_mesh ='2'; % i.e. 2 for main_output_ANT_00002.nc
% in one simulation more than one ROI can be set, think about it
ROI_set = false; % if exist at least one ROI
%ROI='LarsenC';
ROI='Antarctic_Peninsula'; 
plot_Voronoi=true;
%========= END OF CONFIGURATION ======================

% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

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

% load the ice sheet reconstruction from RAISED
RAISED_20ka_small=shaperead('/Users/frre9931/Documents/PhD/RAISED_Shapefiles/RAISED_20ka_AntarcticPeninsula_small.shp');
RAISED_20ka_medium=shaperead('/Users/frre9931/Documents/PhD/RAISED_Shapefiles/RAISED_20ka_AntarcticPeninsula_medium.shp');
RAISED_20ka_big=shaperead('/Users/frre9931/Documents/PhD/RAISED_Shapefiles/RAISED_20ka_AntarcticPeninsula_big.shp');

% filename to load the main output in grid format
filename= [ufe_folder_path, '/main_output_ANT_grid.nc'];

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
Hi_diff=Hi(:,:,end)-Hi(:,:,1);

% calculate mask where Hi>0, outside is NaN
[Hi_fix, maskHi0]= Hi0_to_NaN(Hi);

% do the same for ROI
if ROI_set
filename_ROI= [ufe_folder_path, '/main_output_ANT_grid_ROI_',ROI,'.nc'];
Hi_ROI=ncread(filename_ROI,'Hi');
uabs_surf_ROI=ncread(filename_ROI,'uabs_surf');
x_ROI=ncread(filename_ROI,'x');
y_ROI=ncread(filename_ROI,'y');
[~, maskHi0_ROI]= Hi0_to_NaN(Hi_ROI);
end
%% load the mesh

mesh_path_first= [ufe_folder_path, '/main_output_ANT_00001.nc']; %initial state

% read_mesh_from_file + CL, GL and CF
mesh_first=read_mesh_from_file(mesh_path_first);

time_slice_init=1 ; % first time
time_slice1=ncread(mesh_path_first,'time'); % load the time of the mesh

CL=ncread(mesh_path_first,'coastline',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
GL=ncread(mesh_path_first,'grounding_line',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
CF=ncread(mesh_path_first,'calving_front',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
IM=ncread(mesh_path_first,'ice_margin',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
%GIC=ncread(mesh_path_first,'grounded_ice_contour',[1,1,1],[11788,2,time_slice_init]); % not really useful

%add variables of the mesh during the "last" time,
CL2=ncread(mesh_path_first,'coastline',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
GL2=ncread(mesh_path_first,'grounding_line',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
CF2=ncread(mesh_path_first,'calving_front',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
IM2=ncread(mesh_path_first,'ice_margin',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);

CL3=ncread(mesh_path_first,'coastline',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
GL3=ncread(mesh_path_first,'grounding_line',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
CF3=ncread(mesh_path_first,'calving_front',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
IM3=ncread(mesh_path_first,'ice_margin',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);


if allow_mesh_update
    % path to the new mesh, for now just changing the number here
    mesh_path_update = [ufe_folder_path, '/main_output_ANT_0000',number_mesh,'.nc']; 
    mesh_updated=read_mesh_from_file(mesh_path_update);
    time_slice2=ncread(mesh_path_update,'time');
    CL_mesh2=ncread(mesh_path_update,'coastline',[1,1,1], [11788,2,length(time_slice2)]);
    GL_mesh2=ncread(mesh_path_update,'grounding_line',[1,1,1], [11788,2,length(time_slice2)]);
    CF_mesh2=ncread(mesh_path_update,'calving_front',[1,1,1], [11788,2,length(time_slice2)]);
end
if allow_plot_mesh
    plot_mesh(mesh_first);
    hold on
    plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
    hold off
    if allow_save_plots
        print([path_save,output_folder,'_mesh_1'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh(mesh_updated);
        hold on
        plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
        plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
        hold off
        if allow_save_plots
            print([path_save,output_folder,'_mesh_',number_mesh],'-dpng','-r300')
        end
    else
        plot_mesh(mesh_first);
        hold on
        plot(CF2(:,1),CF2(:,2),'LineWidth',2,'Color','red');
        plot(CL2(:,1),CL2(:,2),'LineWidth',2,'Color','blue');
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','green');
        hold off
    end    
end
%% add plots with Voronois cells
mesh_first.SMB            = ncread( mesh_path_first,'SMB');
mesh_first.Hi            = ncread( mesh_path_first,'Hi');
mesh_first.T2m            = ncread( mesh_path_first,'T2m');
mesh_first.Precip            = ncread( mesh_path_first,'Precip');
mesh_first.uabs              = ncread( mesh_path_first,'uabs_surf');
mesh_first.BMB               = ncread( mesh_path_first,'BMB');
[Hi_fix_mesh, maskHi0_mesh]= Hi0_to_NaN_mesh(mesh_first.Hi);
mesh_first.mask = ncread( mesh_path_first, 'mask');

if allow_mesh_update
    % path to the new mesh, for now just changing the number here
    mesh_path_update = [ufe_folder_path, '/main_output_ANT_0000',number_mesh,'.nc']; 
    mesh_updated=read_mesh_from_file(mesh_path_update);
    time_slice2=ncread(mesh_path_update,'time');
    CL_mesh2=ncread(mesh_path_update,'coastline',[1,1,1], [11788,2,length(time_slice2)]);
    GL_mesh2=ncread(mesh_path_update,'grounding_line',[1,1,1], [11788,2,length(time_slice2)]);
    CF_mesh2=ncread(mesh_path_update,'calving_front',[1,1,1], [11788,2,length(time_slice2)]);
    IM_mesh2=ncread(mesh_path_update,'ice_margin',[1,1,1], [11788,2,length(time_slice2)]);
    mesh_updated.SMB    = ncread(mesh_path_update, 'SMB');
    mesh_updated.Hi     = ncread(mesh_path_update, 'Hi');
    mesh_updated.T2m    = ncread(mesh_path_update, 'T2m');
    mesh_updated.Precip = ncread(mesh_path_update, 'Precip');
    mesh_updated.uabs   = ncread(mesh_path_update,'uabs_surf');
    mesh_updated.BMB    = ncread(mesh_path_update, 'BMB');
    [Hi_fix_mesh_updated, maskHi0_mesh_updated]= Hi0_to_NaN_mesh(mesh_updated.Hi);
end
% polygon for R-LIS
V = [
-0.6469e6, 1.6448e6
-0.6507e6, 1.7370e6
  -0.6411e6, 1.8005e6
  -0.5989e6, 1.8370e6
  -0.5508e6, 1.8639e6
  -0.5104e6, 1.9081e6
  -0.4758e6, 1.9331e6
  -0.4451e6, 1.9542e6
  -0.4393e6, 1.9946e6
  -0.3336e6, 1.9720e6
  -0.3048e6, 1.9292e6
  -0.2644e6, 1.9081e6
  -0.2029e6, 1.8927e6
  -0.1741e6, 1.8716e6
  -0.1558e6, 1.8351e6
  -0.1414e6, 1.8043e6
  -0.1222e6, 1.7659e6
  -0.1057e6, 1.7269e6
  -0.1318e6, 1.6928e6
  -0.1644e6, 1.6640e6
  -0.2125e6, 1.6275e6
  -0.2394e6, 1.5948e6
  -0.2663e6, 1.5833e6
  -0.3259e6, 1.5813e6
  -0.3778e6, 1.5717e6
  -0.4201e6, 1.5640e6
  -0.4528e6, 1.5640e6
  -0.4931e6, 1.5660e6
  -0.5354e6, 1.5698e6
  -0.5758e6, 1.5871e6
  -0.6142e6, 1.6102e6    ];
Tri = 1:length(V);
%fprintf('  poly(%3d,:) = [%.4e_dp,%.4e_dp]\n',[(Tri)',V(:,1),V(:,2)].')
ice_boundaries=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/IceBoundaries_Antarctica_v02.shp');
% plots
% ice thickness
if allow_plot_mesh
    plot_mesh_data(mesh_first,mesh_first.Hi(:,1).*maskHi0_mesh(:,1));
    clim([0 4000]);
    hold on
    %patch('vertices',V,'faces',Tri,'facecolor','none','edgecolor','r','linewidth',3,'marker','o');
    %plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    %plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','red');
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
    t=title(output_folder);
    t.Units='normalized';
    %plot(ice_boundaries(5).X,ice_boundaries(5).Y,'LineWidth',2)
    if allow_save_plots
        print([path_save,output_folder,'_mesh_1_Hi_t0'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh_data(mesh_updated,mesh_updated.Hi(:,length(time_slice2)).*maskHi0_mesh_updated(:,length(time_slice2)));
        clim([0 4000]);
        hold on
        plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');
        t=title(output_folder);
        t.Units='normalized';
    else
        plot_mesh_data(mesh_first,mesh_first.Hi(:,end).*maskHi0_mesh(:,end));
        clim([0 4000]);
        hold on
        %plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
        %plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','red');
        plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
        t=title(output_folder);
        t.Units='normalized';
        if allow_save_plots
            print([path_save,output_folder,'_mesh_1_Hi_tf'],'-dpng','-r300')
        end
    end
end
% ====================
%  plot for velocities
% ====================
if allow_plot_mesh
    plot_mesh_data(mesh_first,log(mesh_first.uabs(:,1)));
    hold on
    plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
    set(gca,'ColorScale','log')
    clim([10^-1 10^1]);
    if allow_save_plots
        %print([path_save,output_folder,'_mesh_1_uabs_t0'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh_data(mesh_updated,log(mesh_updated.uabs(:,length(time_slice2))));
        hold on
        plot(IM_mesh2(:,1),IM_mesh2(:,2),'LineWidth',2,'Color','black');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');
        set(gca,'ColorScale','log')
        clim([10^-1 10^1]);
    else
        plot_mesh_data(mesh_first,log(mesh_first.uabs(:,end)));
        hold on
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','red');
        plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
        set(gca,'ColorScale','log')
        clim([10^-1 10^1]);
        if allow_save_plots
            print([path_save,output_folder,'_mesh_1_uabs_tf'],'-dpng','-r300')
        end
    end
end
%========================
% plot for BMB
% =======================
if allow_plot_mesh
    plot_mesh_data(mesh_first,mesh_first.BMB(:,1));
    hold on
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','black');
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-2 2]);
    t=title(output_folder);
    t.Units='normalized';
    if allow_save_plots
        print([path_save,output_folder,'_mesh_1_BMB_t0'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh_data(mesh_updated,mesh_updated.BMB(:,length(time_slice2)));
        hold on
        plot(IM_mesh2(:,1),IM_mesh2(:,2),'LineWidth',2,'Color','black');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');
        cptcmap('GMT_polar','flip',true,'ncol',100);
        clim([-2 2]);
        if allow_save_plots
          print([path_save,output_folder,'_mesh_2_BMB_tf'],'-dpng','-r300')
        end
    else
        plot_mesh_data(mesh_first,mesh_first.BMB(:,end-10));
        hold on
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','black');
        plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
        plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green','linestyle', ...
            '--');
        plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','red','LineStyle','--');
        cptcmap('GMT_polar','flip',true,'ncol',100);
        clim([-2 2]);
        t=title(output_folder);
        t.Units='normalized';
        if allow_save_plots
            print([path_save,output_folder,'_mesh_1_BMB_tf'],'-dpng','-r300')
        end
    end
end
%% add contour lines from the mesh to my plots
% Set up GUI
wa = 750;
ha = 750;

margins_hor = [125,125];
margins_ver = [125,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);

H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[min(x)+3e4 max(x)-3e4],'ylim',[min(y)+3e4 max(y)-3e4],'fontsize',24,'xgrid','on','ygrid','on');
%figure('position',[100 100 750 750])
hold on
contourf(x,y,Hi_fix(:,:,1)',20,'LineColor','none');
cbar2=colorbar;
%set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
%cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Ice thickness (m)');
t.Units='normalized';
%t.Position(2)=1.05;
clim([0 4000]);
%colormap('jet');
%plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','black');
%plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','black');
if allow_save_plots
    print([path_save,output_folder,'_Hi_t0'],'-dpng','-r300')
end

% now final state
figure('position',[100 100 750 750])
hold on
contourf(x,y,Hi_fix(:,:,end)',20,'LineColor','none');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Ice thickness final state (m)');
t.Units='normalized';
t.Position(2)=1.05;
%colormap('jet');
clim([0 4000]);
if allow_mesh_update
    plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
    plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
    plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black','LineStyle','-.');
    %plot(GIC(:,1),GIC(:,2),'LineWidth',2,'Color','yellow','LineStyle','-.');
else
    plot(CF2(:,1),CF2(:,2),'LineWidth',2,'Color','red');
    plot(CL2(:,1),CL2(:,2),'LineWidth',2,'Color','blue');
    plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','green');
    plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','yellow','LineStyle','-.');
end
if allow_save_plots
    print([path_save,output_folder,'_Hi_tf'],'-dpng','-r300')
end

% difference
figure('position',[100 100 750 750])
hold on
contourf(x,y,Hi_diff(:,:)',20,'LineColor','none');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Ice thickness initial state (m)');
t.Units='normalized';
t.Position(2)=1.05;
clim([min(min(Hi_diff)) max(max(Hi_diff))]);
%colormap('jet');
plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
if allow_save_plots
    print([path_save,output_folder,'_Hi_diff'],'-dpng','-r300')
end

%% now do same plot with velocities adding the maskHi0
figure('position',[100 100 500 500])
hold on
contourf(x,y,log((uabs_surf(:,:,1).*maskHi0(:,:,1))'),20,'LineColor','none');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90]); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Initial surface ice velocity (m/yr)');
t.Units='normalized';
t.Position(2)=1.05;
set(gca,'ColorScale','log')
clim([10^-1 10^1]);
colormap('hsv');
plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
if allow_save_plots
    print([path_save,output_folder,'_uabs_t0'],'-dpng','-r300')
end

figure('position',[100 100 500 500])
hold on
contourf(x,y,log((uabs_surf(:,:,end).*maskHi0(:,:,end))'),50,'LineColor','none');
cbar2=colorbar;
set(gca, 'Position', [0.035, 0.03, 0.83, 0.90],'ColorScale','log'); 
cbar2.Position(1) = cbar2.Position(1) + 0.03;  % Shift it 0.06 units to the right
t=title('Final surface ice velocity (m/yr)');
t.Units='normalized';
t.Position(2)=1.05;
%set(gca,'ColorScale','log')
clim([10^-2 10^1]);
if allow_mesh_update
    plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
    plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
    plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
else
    plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
end
if allow_save_plots
    print([path_save,output_folder,'_uabs_tf'],'-dpng','-r300')
end
%% NEXT TASK IS TO WORK ON THE SPECIFIC REGIONS, FOR EXAMPLE A PLOT FOR AP
if ROI_set

% Set up GUI
wa = 750;
ha = 750;

margins_hor = [125,125];
margins_ver = [125,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);

H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[min(x_ROI)+5e4 max(x_ROI)-5e4],'ylim',[min(y_ROI)+5e4 max(y_ROI)-5e4],'fontsize',24,'xgrid','on','ygrid','on');
hold on
contourf(x_ROI,y_ROI,(Hi_ROI(:,:,1).*maskHi0_ROI(:,:,1))',20,'LineColor','none');
cbar2=colorbar;
t=title('Ice thickness initial state (m)');
t.Units='normalized';
clim([0 (max(max(Hi_ROI(:,:,1)))+max(max(Hi_ROI(:,:,end))))/2]);

plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');

if allow_save_plots
    print([path_save,output_folder,'_Hi_t0_',ROI],'-dpng','-r300')
end

% now final state
H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[min(x_ROI)+5e4 max(x_ROI)-5e4],'ylim',[min(y_ROI)+5e4 max(y_ROI)-5e4],'fontsize',24,'xgrid','on','ygrid','on');
hold on
contourf(x_ROI,y_ROI,(Hi_ROI(:,:,end).*maskHi0_ROI(:,:,end))',20,'LineColor','none');
cbar2=colorbar;
t=title('Ice thickness final state (m)');
t.Units='normalized';
clim([0 (max(max(Hi_ROI(:,:,1)))+max(max(Hi_ROI(:,:,end))))/2]);

if allow_mesh_update
    plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
    plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
    plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black','LineStyle','-.');
    plot(RAISED_20ka_small.X,RAISED_20ka_small.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
    plot(RAISED_20ka_medium.X,RAISED_20ka_medium.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
    plot(RAISED_20ka_big.X,RAISED_20ka_big.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
else
    plot(CF2(:,1),CF2(:,2),'LineWidth',2,'Color','red');
    plot(CL2(:,1),CL2(:,2),'LineWidth',2,'Color','blue');
    plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','green');
    plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black','LineStyle','-.');
    plot(RAISED_20ka_small.X,RAISED_20ka_small.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
    plot(RAISED_20ka_medium.X,RAISED_20ka_medium.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
    plot(RAISED_20ka_big.X,RAISED_20ka_big.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
end
if allow_save_plots
    print([path_save,output_folder,'_Hi_tf_',ROI],'-dpng','-r300')
end

end % if ROI_set
%% Implement a plot with a zoom into a desired regions, without ROI.

% Plot range for Antarctic Peninsula
xmid = -2200e3;
ymid =  1200e3;
wx    = 600e3;
wy    = 800e3;

xmin = xmid - wx;
xmax = xmid + wx;
ymin = ymid - wy;
ymax = ymid + wy;

wa = 750;
ha = 750;

margins_hor = [125,125];
margins_ver = [125,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);
if plot_Voronoi
edgecolor = 'none';
mesh=mesh_updated;
d=mesh_updated.Hi(:,length(time_slice2)).*maskHi0_mesh_updated(:,length(time_slice2));
H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[xmin,xmax],'ylim',[ymin,ymax],'fontsize',24,'xgrid','on','ygrid','on');

if isfield( mesh, 'VVor')
  H.Patch = patch('vertices',mesh.Vor,'faces',changem(double(mesh.VVor),NaN),...
    'facecolor','flat','facevertexcdata',d,'edgecolor',edgecolor);
else
  H.Patch = patch('vertices',mesh.V( 1:mesh.nV,:),'faces',mesh.Tri( 1:mesh.nTri,:),...
    'facecolor','interp','facevertexcdata',d,'edgecolor',edgecolor);
end

pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);
set( H.Ax,'units','normalized');
hold on
plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');

else
H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[xmin,xmax],'ylim',[ymin,ymax],'fontsize',24,'xgrid','on','ygrid','on');
hold on
contourf(x,y,Hi_fix(:,:,end)',20,'LineColor','none');
end
plot(RAISED_20ka_small.X,RAISED_20ka_small.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
plot(RAISED_20ka_medium.X,RAISED_20ka_medium.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
plot(RAISED_20ka_big.X,RAISED_20ka_big.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black','LineStyle','-.');

cbar2=colorbar;
clim([0 2500]);
t=title('Ice thickness - AWIESM1 (m)');
t.Units='normalized';
hold off
if allow_save_plots
    print([path_save,output_folder,'_Hi_tf_AP2'],'-dpng','-r300')
end
uv_step=1:10:length(x);
%quiver(x(uv_step),y(uv_step),(u_surf(uv_step,uv_step,end).*maskHi0(uv_step,uv_step,end))',...
%    (v_surf(uv_step,uv_step,end).*maskHi0(uv_step,uv_step,end))',10,'LineStyle','-',...
%    'color','black');

% repeat same plot with uabs
H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[xmin,xmax],'ylim',[ymin,ymax],'fontsize',24,'xgrid','on','ygrid','on');
hold on
contourf(x,y,(uabs_surf(:,:,end).*maskHi0(:,:,end))',20,'LineColor','none');
plot(RAISED_20ka_small.X,RAISED_20ka_small.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
plot(RAISED_20ka_medium.X,RAISED_20ka_medium.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
plot(RAISED_20ka_big.X,RAISED_20ka_big.Y,'LineWidth',2,'Color','magenta','LineStyle','-.');
plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black','LineStyle','-.');
cbar2=colorbar;
clim([0 1000]);
t=title('Surface ice velocity - AWIESM1 (m/yr)');
t.Units='normalized';
if allow_save_plots
    print([path_save,output_folder,'_uabs_tf_AP'],'-dpng','-r300')
end
%% create an artifical mask for ice shelf ISMIP6 experiment, just for testing
% create a mask from grid BMB, basically were it exist BMB put a 1 bcs is
% an ice shelf
faketime=[2000:10:2300];
shelf_mask=ones([size(BMB(:,:,1)) length(faketime)]);
% I need to stack this with time variable... because the input is expected
% to have time!
mask_ismip6=shelf_mask;
mask_ismip6(:,:,1)=shelf_mask(:,:,1).*BMB(:,:,1);
pos_val1=find(mask_ismip6(:,:,1) > 0 | mask_ismip6(:,:,1)<0);
mask_ismip6(pos_val1)=1;
for k=1:length(faketime)
    mask_ismip6(:,:,k)=mask_ismip6(:,:,1);
end
% save netcdf file using x,y from file
ncid = netcdf.create(['/Users/frre9931/Documents/PhD/ANT_UFEMISM/mask_ismip6_test.nc'],'CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,1));
dim_y = netcdf.defDim(ncid,'y',size(y,1));
dim_time = netcdf.defDim(ncid,'time',size(faketime,2));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
id_time = netcdf.defVar(ncid,'time','double',dim_time);
id_mask_ismip6  = netcdf.defVar(ncid,'mask','double',[dim_x, dim_y, dim_time]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x'); 
netcdf.putAtt(ncid,id_y,'standard_name','y');

netcdf.putAtt(ncid,id_mask_ismip6,'standard_name','Ice shelf mask to melt'); 
netcdf.putAtt(ncid,id_time,'standard_name','time');
netcdf.putAtt(ncid,id_time,'units','years');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);
netcdf.putVar(ncid,id_time,faketime);

netcdf.putVar(ncid,id_mask_ismip6,mask_ismip6);

% Close file
% ==========

netcdf.close(ncid);














%%
% figure of Hb difference
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