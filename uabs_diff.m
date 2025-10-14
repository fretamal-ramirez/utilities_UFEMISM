% plot to create differences between control simulation and another run
% IMPORTANT: both should have the same mesh!

clear all; clc;

%========= PATH TO OUTPUT UFEMISM DIRECTORY ==========
control_folder = 'results_ant_PD_smooth_Hb_ctrl2500_PMP';
output_folder = 'results_ant_PD_smooth_Hb_retreat_mask_PMP';
folder_path= '/Volumes/One Touch/results_UFEMISM/tetralith_results/';
%ufe_folder_path=['/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/', output_folder];
%ufe_folder_path=['/Users/frre9931/Desktop/UFEMISM2.0_porting/', output_folder];
%ufe_folder_path=['/Users/frre9931/Desktop/tetralith_results/', output_folder];
ufe_ctrl_path=[folder_path, control_folder];
ufe_folder_path=[folder_path, output_folder];
allow_plot_mesh = true; % if we want to plot the mesh
allow_save_plots = false;
path_save = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/';

%========= END OF CONFIGURATION ======================

% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

% add files from MEaSUREs
coast_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Coastline_Antarctica_v02.shp');
% load the basins from MEaSUREs
basins_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Basins_Antarctica_v02.shp');

%% load the mesh
mesh_path_ctrl= [ufe_ctrl_path, '/main_output_ANT_00001.nc'];
mesh_path_first= [ufe_folder_path, '/main_output_ANT_00001.nc']; %initial state

% read_mesh_from_file + CL, GL and CF
mesh_ctrl=read_mesh_from_file(mesh_path_ctrl);
mesh_first=read_mesh_from_file(mesh_path_first);

time_slice_init=1 ; % first time
time_slice1=ncread(mesh_path_first,'time'); % load the time of the mesh

CL=ncread(mesh_path_first,'coastline',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
GL=ncread(mesh_path_first,'grounding_line',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
CF=ncread(mesh_path_first,'calving_front',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
IM=ncread(mesh_path_first,'ice_margin',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);

%add variables of the mesh during the "last" time,
CL2=ncread(mesh_path_first,'coastline',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
GL2=ncread(mesh_path_first,'grounding_line',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
CF2=ncread(mesh_path_first,'calving_front',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
IM2=ncread(mesh_path_first,'ice_margin',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);

%% add plots with Voronois cells
mesh_ctrl.SMB            = ncread( mesh_path_ctrl,'SMB');
mesh_ctrl.Hi            = ncread( mesh_path_ctrl,'Hi');
mesh_ctrl.T2m            = ncread( mesh_path_ctrl,'T2m');
mesh_ctrl.Precip            = ncread( mesh_path_ctrl,'Precip');
mesh_ctrl.uabs              = ncread( mesh_path_ctrl,'uabs_surf');
mesh_ctrl.BMB               = ncread( mesh_path_ctrl,'BMB');

mesh_first.SMB            = ncread( mesh_path_first,'SMB');
mesh_first.Hi            = ncread( mesh_path_first,'Hi');
mesh_first.T2m            = ncread( mesh_path_first,'T2m');
mesh_first.Precip            = ncread( mesh_path_first,'Precip');
mesh_first.uabs              = ncread( mesh_path_first,'uabs_surf');
mesh_first.BMB               = ncread( mesh_path_first,'BMB');

% ====================
%  plot for velocities
% ====================
if allow_plot_mesh
    plot_mesh_data_b_RLIS(mesh_first,mesh_ctrl.uabs(:,end)-mesh_first.uabs(:,end));
    hold on
    %plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',2,'Color','black'); %R-LIS
    %plot(basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y,'LineWidth',2,'Color','black'); % Brunt
    plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
    %plot(coast_MEaSUREs.X,coast_MEaSUREs.Y,'LineWidth',2,'Color','black'); % Brunt
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',false,'ncol',256);
    clim([-500 500])
    cb = colorbar;
    cb.Label.String = 'Surface ice velocity diffence (m/yr)';
    if allow_save_plots
        %print([path_save,output_folder,'_diff-control_uabs'],'-dpng','-r300')
    end
end