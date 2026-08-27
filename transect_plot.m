% code to plot transects using UFEMISM outputs
clear all; clc;

%========= PATH TO OUTPUT UFEMISM DIRECTORY ==========
output_folder = 'results_ant_PD_maxphi_Hb-2000to-250m_SHR_ctrl2500_ocndT_1e1_N';
ufe_folder_path=['/Users/frre9931/Desktop/arrhenius_results/', output_folder];

allow_plot_mesh = true; % if we want to plot the mesh
allow_save_plots = false;
path_save = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/';
allow_mesh_update = false; % if remeshing is allowed in simulation
% number pointing the file with updated mesh
number_mesh ='1'; % i.e. 2 for main_output_ANT_00002.nc
plot_Voronoi=true;
%========= END OF CONFIGURATION ======================

% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

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

% CL3=ncread(mesh_path_first,'coastline',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
% GL3=ncread(mesh_path_first,'grounding_line',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
% CF3=ncread(mesh_path_first,'calving_front',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
% IM3=ncread(mesh_path_first,'ice_margin',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);


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
        %print([path_save,output_folder,'_mesh_1'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh(mesh_updated);
        hold on
        plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
        plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
        hold off
        if allow_save_plots
            %print([path_save,output_folder,'_mesh_',number_mesh],'-dpng','-r300')
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

% add variables on Voronois cells
mesh_first.Hi            = ncread( mesh_path_first,'Hi');
mesh_first.uabs              = ncread( mesh_path_first,'uabs_surf');
mesh_first.Hib               = ncread( mesh_path_first,'Hib');
mesh_first.Hs               = ncread( mesh_path_first,'Hs');
mesh_first.Hb                 = ncread( mesh_path_first,'Hb');
[Hi_fix_mesh, maskHi0_mesh]= Hi0_to_NaN_mesh(mesh_first.Hi);
mesh_first.mask = ncread( mesh_path_first, 'mask');
mesh_first.till_friction_angle = ncread( mesh_path_first,'till_friction_angle');

% calculate Hi differences
mesh_first.Hi_diff=(mesh_first.Hi(:,end))-(mesh_first.Hi(:,1));
% polygon for catchments
ice_boundaries=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/IceBoundaries_Antarctica_v02.shp');
GL_time = ncread(mesh_path_first,'grounding_line',[1,1,1], [size(mesh_first.E,1),2,11]);
% need to interpolate coordinates of GL in the transect for the plot
%% load transect using .trs file
transect_xy = [-3.4185e+05 1.6987e+06
-3.5170e+05 1.6972e+06
-3.6647e+05 1.6967e+06
-3.8469e+05 1.6994e+06
-4.0242e+05 1.7006e+06
-4.2042e+05 1.7047e+06
-4.4065e+05 1.7039e+06
-4.5440e+05 1.6974e+06
-4.6290e+05 1.6865e+06
-4.6937e+05 1.6820e+06]; % copied and pasted from VeststraumenA for now

% Distance between consecutive points
dx = diff(transect_xy(:,1));
dy = diff(transect_xy(:,2));

ds = sqrt(dx.^2 + dy.^2);
spacing = 5000; % 5km

% Cumulative distance from the beginning of the transect
distance_original = [0; cumsum(ds)];

% Total length of transect
total_distance = distance_original(end);

% New regularly spaced distance coordinates
distance = (0:spacing:total_distance)';

% Make sure the final point is included
if distance(end) < total_distance
    distance = [distance; total_distance];
end

% Interpolate x and y coordinates onto the regular distance vector
transect_x = interp1(distance_original,transect_xy(:,1),distance,'linear');
transect_y = interp1(distance_original,transect_xy(:,2),distance,'linear');

% New transect coordinates
transect_xy_reg = [transect_x transect_y];

for i=1:length(time_slice1) % loop on time
interpolant_hi = scatteredInterpolant(mesh_first.V(:,1),mesh_first.V(:,2),mesh_first.Hi(:,i),"linear",'none');
interpolant_hib = scatteredInterpolant(mesh_first.V(:,1),mesh_first.V(:,2),mesh_first.Hib(:,i),"linear",'none');
interpolant_hs = scatteredInterpolant(mesh_first.V(:,1),mesh_first.V(:,2),mesh_first.Hs(:,i),"linear",'none');
interpolant_vel = scatteredInterpolant(mesh_first.TriGC(:,1),mesh_first.TriGC(:,2),mesh_first.uabs(:,i),"linear",'none');
transect_hi(:,i)= interpolant_hi(transect_xy_reg(:,1),transect_xy_reg(:,2));
transect_hib(:,i)= interpolant_hib(transect_xy_reg(:,1),transect_xy_reg(:,2));
transect_hs(:,i)= interpolant_hs(transect_xy_reg(:,1),transect_xy_reg(:,2));
transect_uabs(:,i) = interpolant_vel(transect_xy_reg(:,1),transect_xy_reg(:,2));
end
distance=distance/1000; % convert to km for plots
%% create a video for the transect
% ==== VIDEO SETUP ====
video_name = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/transect_Veststraumen_50yr.mp4';
v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 1;
v.Quality = 100;
open(v);
fig=figure('Position',[0 0 800 800]);

for i=1:length(time_slice1) % loop on time
clf(fig)
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(distance,transect_hib(:,i))
hold on
plot(distance,transect_hs(:,i))
hold off
sgtitle(['OC 1.0 - Year ',num2str(time_slice1(i))])
set(gca, 'XDir','reverse')
xlabel('Distance along transect (km)');
ylabel('Elevation (m)');
%ylim([0 1600])
ylim([-800 1200])
xlim([0 max(distance)])
grid on
legend('Hib','Hs','Location','northwest')

% velocities
nexttile
plot(distance,transect_uabs(:,i))
xlabel('Distance along transect (km)');
ylabel('Surface velocity (m/a)');
ylim([0 1000])
xlim([0 max(distance)])
grid on
set(gca, 'XDir','reverse')

drawnow
frame = getframe(fig);
writeVideo(v,frame);
end
close(v);
close(fig);