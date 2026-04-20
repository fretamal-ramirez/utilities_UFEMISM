%% MULTIPANEL COMPARISON PLOTS (shared colorbars)
% Plot to zoom in RLIS and show GLs of PD and simulation
clear all; close all; clc;

% ==== DEFINE OUTPUTS ====
outputs = { ...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_max30_SHR_new'...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_M11_Hb-2000to-250m_SHR',...
};
tile_size = 300; % pixels for each panel
nCols = numel(outputs);
nRows = 2; % uabs_diff, Hi_diff

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

xmin = -878893.0; xmax = 131107.0;
ymin = 1137850.0; ymax = 2217850.0;

fig=figure('Units','pixels','Position',[100 100 fig_width+100 fig_height],'Visible','off');
tiledlayout(nRows,nCols,"TileSpacing","compact","Padding","compact");

% ==== PATHS ====
%basepath = '/Volumes/One Touch/results_UFEMISM/tetralith_results/';
basepath = '/Users/frre9931/Desktop/tetralith_results/';
colormaps.devon = '/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt';
plot_titles = {'Northings (m)', 'Northings (m)'};
letters_for_plots = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'}; % 12 plot max for now

% ==== Add extra functions ====
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
%path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
%path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

% ===== Data to complement plots ====
rock_outcrops=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ADD_RockOutcrops_RLIS.shp');

% Preallocate axes for later shared colorbars
ax_all = gobjects(nRows,nCols);
for k = 1:nCols
output_folder=outputs{k};
filepath = [basepath,output_folder,'/main_output_ANT_00001.nc'];

% Load MEaSUREs data
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

% === Load mesh ===
mesh = read_mesh_from_file(filepath);
mesh.uabs = ncread(filepath,'uabs_surf');
mesh.BMB  = ncread(filepath,'BMB');
mesh.Hi   = ncread(filepath,'Hi');
mesh.mask = ncread(filepath,'mask');
mesh.Hb   = ncread(filepath, 'Hb');
mesh.tfa  = ncread(filepath, 'till_friction_angle');
[Hi_fix, maskHi0] = Hi0_to_NaN_mesh(mesh.Hi);
mesh.Hi_diff = mesh.Hi(:,end) - mesh.Hi(:,1);
mesh.Hb_diff = mesh.Hb(:,end) - mesh.Hb(:,1);

for i=1:size(mesh.BMB,1)
    for j=1:size(mesh.BMB,2)
        if mesh.BMB(i,j)==0
            mesh.BMB(i,j)=NaN;
        end
    end
end
    % === Grounding & ice margins ===
nE = size(mesh.E,1);
nt = length(ncread(filepath,'time'));
GL1 = ncread(filepath,'grounding_line',[1,1,1],[nE,2,1]);
GL2 = ncread(filepath,'grounding_line',[1,1,nt],[nE,2,1]);
IM2 = ncread(filepath,'ice_margin',[1,1,nt],[nE,2,1]);

% info from ROI file
uabs_ufe=ncread([basepath,output_folder, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'uabs_surf');
x_ufe=ncread([basepath,output_folder, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'x');
y_ufe=ncread([basepath,output_folder, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'y');
Hi_ufe=ncread([basepath,output_folder, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'Hi');
BMB_ufe=ncread([basepath,output_folder, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'BMB');
x_ufe_min=min(x_ufe); x_ufe_max= max(x_ufe);
y_ufe_min=min(y_ufe); y_ufe_max= max(y_ufe);

% create a mask where there is ice, to mask the velocity later
[~, maskHi0_ufe]= Hi0_to_NaN(Hi_ufe);
[xx_measures,yy_measures]=meshgrid(x_MEaSUREs,y_MEaSUREs);
[xx_ufe,yy_ufe]=meshgrid(x_ufe,y_ufe);
% interpolate MEaSUREs to same resolution as UFEMISM output
uabs_measures_interp=interp2(xx_measures,yy_measures,uabs_MEaSUREs',xx_ufe,yy_ufe);
%uabs_diff=uabs_measures_interp-(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))';
uabs_diff=(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))'-uabs_measures_interp;

% fill 0s from BMB with NaNs in UFE data
for i=1:size(BMB_ufe,1)
    for j=1:size(BMB_ufe,2)
        if BMB_ufe(i,j,end) == 0
            BMB_ufe(i,j,end)=NaN;
        end
    end
end

% BMB observation data
BMB_Davison=ncread('/Users/frre9931/Downloads/Davison_ice_shelf_mass_budget/data/basal_melt/basal_melt_map_racmo_firn_air_corrected.nc','Band1');
x_Davison=ncread('/Users/frre9931/Downloads/Davison_ice_shelf_mass_budget/data/basal_melt/basal_melt_map_racmo_firn_air_corrected.nc','x');
y_Davison=ncread('/Users/frre9931/Downloads/Davison_ice_shelf_mass_budget/data/basal_melt/basal_melt_map_racmo_firn_air_corrected.nc','y');
[xx_Davison,yy_Davison]=meshgrid(x_Davison,y_Davison);
% interpolate Davison'bmb to same resolution as UFEMISM output
bmb_davison_interp=interp2(xx_Davison,yy_Davison,BMB_Davison',xx_ufe,yy_ufe);
bmb_diff=-1*((BMB_ufe(:,:,end).*maskHi0_ufe(:,:,end))')-bmb_davison_interp;

% calculate statistics inside polygons

% load ROI polygon from QGIS
ROI_polygon=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ROI_for_R-LIS.shp');
% load shapefiles with basins
basins_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Basins_Antarctica_v02.shp');
% original file with ice shelves, that delineates RLIS and BIS, used just
% for plot
iceshelves_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/IceShelf_Antarctica_v02.shp');
% load simplified calving front shapefiles to not cut pixels in Ufe
shelf_RLIS= shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/simplified_shelf_RLIS.shp');
shelf_BIS= shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/simplified_shelf_BIS.shp');

% check points that are inside of each polygon for the gridded output
in_ROI=inpolygon(xx_ufe,yy_ufe,ROI_polygon.X,ROI_polygon.Y);
% basins_MEaSUREs 3 for BIS and 4 for RLIS
in_basin_RLIS=inpolygon(xx_ufe,yy_ufe,basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y);
in_basin_BIS=inpolygon(xx_ufe,yy_ufe,basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y);
in_shelf_RLIS=inpolygon(xx_ufe,yy_ufe,shelf_RLIS.X,shelf_RLIS.Y);
in_shelf_BIS=inpolygon(xx_ufe,yy_ufe,shelf_BIS.X,shelf_BIS.Y);

diff_uabs_in_ROI=uabs_diff.*in_ROI;
diff_uabs_in_basin_RLIS=uabs_diff.*in_basin_RLIS;
diff_uabs_in_basin_BIS=uabs_diff.*in_basin_BIS;
diff_uabs_in_shelf_RLIS=uabs_diff.*in_shelf_RLIS;
diff_uabs_in_shelf_BIS=uabs_diff.*in_shelf_BIS;

for i=1:size(in_ROI,1)
    for j=1:size(in_ROI,2)
        if diff_uabs_in_ROI(i,j) == 0
            diff_uabs_in_ROI(i,j)=NaN;
        elseif diff_uabs_in_basin_RLIS(i,j) == 0
            diff_uabs_in_basin_RLIS(i,j) = NaN;
        elseif diff_uabs_in_basin_BIS(i,j) == 0
            diff_uabs_in_basin_BIS(i,j) = NaN;
        elseif diff_uabs_in_shelf_RLIS(i,j) == 0
            diff_uabs_in_shelf_RLIS(i,j) = NaN;
        elseif diff_uabs_in_shelf_BIS(i,j) == 0
            diff_uabs_in_shelf_BIS(i,j) = NaN;
        end
    end
end
% same for vec/mesh for the Hi
in_ROI_mesh=inpolygon(mesh.V(:,1),mesh.V(:,2),ROI_polygon.X,ROI_polygon.Y);
in_basin_RLIS_mesh = inpolygon(mesh.V(:,1),mesh.V(:,2),basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y);
in_basin_BIS_mesh = inpolygon(mesh.V(:,1),mesh.V(:,2),basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y);
in_shelf_RLIS_mesh = inpolygon(mesh.V(:,1),mesh.V(:,2),shelf_RLIS.X,shelf_RLIS.Y);
in_shelf_BIS_mesh = inpolygon(mesh.V(:,1),mesh.V(:,2),shelf_BIS.X,shelf_BIS.Y);

mesh.Hi_diff_in_ROI=mesh.Hi_diff.*in_ROI_mesh;
Hi_diff_in_basin_RLIS = mesh.Hi_diff.*in_basin_RLIS_mesh;
Hi_diff_in_basin_BIS = mesh.Hi_diff.*in_basin_BIS_mesh;
Hi_diff_in_shelf_RLIS = mesh.Hi_diff.*in_shelf_RLIS_mesh;
Hi_diff_in_shelf_BIS = mesh.Hi_diff.*in_shelf_BIS_mesh;

for i=1:length(mesh.Hi_diff_in_ROI)
    if mesh.Hi_diff_in_ROI(i) == 0
        mesh.Hi_diff_in_ROI(i) = NaN;
    elseif Hi_diff_in_basin_RLIS(i) == 0
        Hi_diff_in_basin_RLIS(i) = NaN;
    elseif Hi_diff_in_basin_BIS(i) == 0
        Hi_diff_in_basin_BIS(i) = NaN;
    elseif Hi_diff_in_shelf_RLIS(i) == 0
        Hi_diff_in_shelf_RLIS(i) = NaN;
    elseif Hi_diff_in_shelf_BIS(i) == 0
        Hi_diff_in_shelf_BIS(i) = NaN;
    end
end

% calculate the RMSEs
% rmse for velocites
rmseU=rmse(uabs_measures_interp,(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))','all','omitnan');
rmseU_in_ROI=rmse(uabs_measures_interp,(uabs_ufe(:,:,end)'.*(diff_uabs_in_ROI./diff_uabs_in_ROI)),'all','omitnan');
rmseU_in_basin_RLIS = rmse(uabs_measures_interp,(uabs_ufe(:,:,end)'.*(diff_uabs_in_basin_RLIS./diff_uabs_in_basin_RLIS)),'all','omitnan');
rmseU_in_basin_BIS = rmse(uabs_measures_interp,(uabs_ufe(:,:,end)'.*(diff_uabs_in_basin_BIS./diff_uabs_in_basin_BIS)),'all','omitnan');
rmseU_in_shelf_RLIS = rmse(uabs_measures_interp,(uabs_ufe(:,:,end)'.*(diff_uabs_in_shelf_RLIS./diff_uabs_in_shelf_RLIS)),'all','omitnan');
rmseU_in_shelf_BIS = rmse(uabs_measures_interp,(uabs_ufe(:,:,end)'.*(diff_uabs_in_shelf_BIS./diff_uabs_in_shelf_BIS)),'all','omitnan');

% rmse for ice thickness
rmseH=sqrt(sum(mesh.Hi_diff(:,1).^2)/length(mesh.Hi_diff(:,1)));
rmseH_in_ROI=rmse(mesh.Hi(:,1),(mesh.Hi(:,end).*(mesh.Hi_diff_in_ROI./mesh.Hi_diff_in_ROI)),'all','omitnan');
rmseHi_in_basin_RLIS = rmse(mesh.Hi(:,1),(mesh.Hi(:,end).*(Hi_diff_in_basin_RLIS./Hi_diff_in_basin_RLIS)),'all','omitnan');
rmseHi_in_basin_BIS = rmse(mesh.Hi(:,1),(mesh.Hi(:,end).*(Hi_diff_in_basin_BIS./Hi_diff_in_basin_BIS)),'all','omitnan');
rmseHi_in_shelf_RLIS = rmse(mesh.Hi(:,1),(mesh.Hi(:,end).*(Hi_diff_in_shelf_RLIS./Hi_diff_in_shelf_RLIS)),'all','omitnan');
rmseHi_in_shelf_BIS = rmse(mesh.Hi(:,1),(mesh.Hi(:,end).*(Hi_diff_in_shelf_BIS./Hi_diff_in_shelf_BIS)),'all','omitnan');

% print the RMSEs as a table in the terminal
fprintf('RMSE for :    RL basin | RLIS | SW basin | BIS\n');
fprintf('Ice velocity  = %.0f | %.0f | %.0f | %.0f\n', rmseU_in_basin_RLIS, rmseU_in_shelf_RLIS, rmseU_in_basin_BIS, rmseU_in_shelf_BIS);
fprintf('Ice thickness = %.0f | %.0f | %.0f | %.0f\n', rmseHi_in_basin_RLIS, rmseHi_in_shelf_RLIS, rmseHi_in_basin_BIS, rmseHi_in_shelf_BIS);

% masks outside the ice sheet
mask_current = mesh.Hi(:,end) > 0;

% mask_combined = mask_current & ~mask_ROI;
% 
% % --- Apply mask ---
% Hi_diff(~mask_initial) = NaN;
% 
% uabs_initial = mesh.uabs(:,1);
% uabs_current = mesh.uabs(:,tid);
% 
% uabs_diff = uabs_current - uabs_initial;

% Project node mask to triangle centers
Fmask = scatteredInterpolant(mesh.V(:,1), mesh.V(:,2), ...
                     double(mask_current), ...
                     'nearest','none');

mask_vel = Fmask(mesh.TriGC(:,1),mesh.TriGC(:,2));
uabs_masked = mesh.uabs(:,end);
uabs_masked(mask_vel==0)=NaN;

GL_forced = ncread([basepath,'results_ant_PD_maxphi_Hb-2000to-250m_SHR_retreat_ocndT_2e1', '/main_output_ANT_00001.nc'],'grounding_line',[1,1,7],[nE,2,1]);
% %
% plot
% %

    % =================================
    % --- 1,1 VELOCITY ROW IN MESH ---
    % =================================
    ax_all(1,k) = nexttile(1 + (k-1)); % row 1

    % if velocities as they come
    % plot_mesh_data_b_onlyRLIS(mesh, log10(uabs_masked), ax_all(1,1));
    % hold on;
    % plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',0.8,'Color','green'); %R-LIS
    % plot(GL_forced(:,1),GL_forced(:,2),'r','LineWidth',0.8,'LineStyle','-');
    % plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % cptcmap(colormaps.devon,'flip',false,'ncol',100);
    % clim([log10(1) log10(2000)]);
    % cb=colorbar;
    % cb.Ticks = log10([1 10 50 100 500 1000 2000]);
    % cb.TickLabels = {'1','10','50','100','500','1000','2000'};
    % cb.Label.String = 'Velocity (m/yr)';
    % cb.Label.FontSize = 12;

    % if velocity difference
    contourf(x_ufe,y_ufe,uabs_diff,100,'LineColor','none');
    hold on;
    plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',0.8,'Color','green'); %R-LIS
    plot(GL_forced(:,1),GL_forced(:,2),'r','LineWidth',0.8,'LineStyle','-');
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',false,'ncol',256);
    clim([-800 800]);
    set(ax_all(1,k),'XLim',[xmin,xmax],'YLim',[ymin,ymax],...
        'XGrid','on','YGrid','on');
    box off
    %cb=colorbar;    
    %cb.Label.String = 'Velocity difference (m/yr)';
    %cb.Label.FontSize = 12;
    if k == 1
        ylabel(plot_titles{1},'FontWeight','bold');
    end  
    text(-7.8e5,2.13e6,letters_for_plots{k},'FontWeight','bold','FontSize',13);
    %text(-6.0e5,1.95e6,'(a)','FontWeight','bold','FontSize',13);
    % =================================
    % --- 2,1 Hi DIFF RLIS IN MESH ---
    % =================================
    ax_all(2,k) = nexttile(nCols + k); % row 2    
    plot_mesh_data_a_RLIS(mesh, mesh.Hi_diff(:,1).*maskHi0(:,end), ax_all(2,k));
    %plot_mesh_data_a_RLIS(mesh, Hi_diff_in_basin_RLIS, ax_all(2,2));
    hold on;
    plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',0.8,'Color','green'); %R-LIS
    plot(GL_forced(:,1),GL_forced(:,2),'LineWidth',0.8,'LineStyle','-','Color','red');
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    %plot(iceshelves_MEaSUREs(37).X,iceshelves_MEaSUREs(37).Y,'LineWidth',0.8,'Color','black'); %R-LIS
    %plot(iceshelves_MEaSUREs(36).X,iceshelves_MEaSUREs(36).Y,'LineWidth',0.8,'Color','black'); % Brunt
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',true,'ncol',256);
    clim([-500 500]);
    xlabel('Eastings (m)','FontWeight','bold');
    %cb=colorbar;
    %cb.Label.String = 'Ice thickness difference (m)';
    %cb.Label.FontSize = 12;
    text(-7.8e5,2.13e6,letters_for_plots{k+nCols},'FontWeight','bold','FontSize',13);
    if k == 1
        ylabel(plot_titles{2},'FontWeight','bold');
    end   
end

% ===============================
% === SHARED COLORBARS PER ROW ===
% ===============================
for r = 1:nRows
    row_axes = ax_all(r,:);  % all axes in row r
    % Set the same CLim for all axes in this row
    clim_values = cell2mat(get(row_axes,'CLim'));
    minCL = min(clim_values(:,1));
    maxCL = max(clim_values(:,2));
    set(row_axes,'CLim',[minCL maxCL]);

    % Add one colorbar to the last axes of the row
    cb = colorbar(row_axes(end),'eastoutside');
    switch r
        case 1
            %cb.Ticks = log10([1 10 50 100 500 1000 2000]);
            %cb.TickLabels = {'1','10','50','100','500','1000','2000'};
            cb.Label.String = 'Surface velocity difference (m/yr)';
            cb.Label.FontSize = 12;
        case 2
            %cb.Label.String = 'BMB (m/yr)';
            cb.Label.String = 'Ice thickness difference (m)';
            cb.Label.FontSize = 12;
            %cb.Label.String = 'Till friction angle (°)';
    end
end
print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/multipanel_RLIS.png','-dpng','-r300');
%% repeat some of the things above to create only the simulation towards 2300
clear all;
outputs = { ...
    'results_ant_PD_maxphi_Hb-2000to-250m_SHR_retreat_ocndT_2e1',...
};

tile_size = 300; % pixels for each panel
nCols = 1;
nRows = 1; % Hi_diff
timeslice = 2300; % time to be plotted

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

% ==== PATHS ====
%basepath = '/Volumes/One Touch/results_UFEMISM/tetralith_results/';
basepath = '/Users/frre9931/Desktop/tetralith_results/';
colormaps.devon = '/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt';
plot_titles = {'Velocity (m/yr)', 'Basal melt rate (m/yr)', 'ΔIce thickness (m)'};
%plot_titles = {'Velocity (m/yr)', 'Basal melt rate (m/yr)', 'ΔIce thickness (m)', 'Till friction angle (°)'};
%plot_titles = {'Velocity (m/yr)', 'Till friction angle (°)', 'ΔIce thickness (m)', 'Final ice thickness (m)'};

% ==== Add extra functions ====
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
%path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
%path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

% ===== Data to complement plots ====
rock_outcrops=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ADD_RockOutcrops_RLIS.shp');
% load shapefiles with basins
basins_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Basins_Antarctica_v02.shp');

% Preallocate axes for later shared colorbars
ax_all = gobjects(nRows,nCols);
output_folder=outputs{1};
filepath = [basepath,output_folder,'/main_output_ANT_00001.nc'];

% Load MEaSUREs data
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

modeltime = ncread(filepath, 'time');
time_to_plot = find(modeltime==timeslice);
% === Load mesh ===
mesh = read_mesh_from_file(filepath);
mesh.uabs = ncread(filepath,'uabs_surf');
mesh.BMB  = ncread(filepath,'BMB');
mesh.Hi   = ncread(filepath,'Hi');
mesh.mask = ncread(filepath,'mask');
mesh.Hb   = ncread(filepath, 'Hb');
mesh.tfa  = ncread(filepath, 'till_friction_angle');
[Hi_fix, maskHi0] = Hi0_to_NaN_mesh(mesh.Hi);
mesh.Hi_diff = mesh.Hi(:,time_to_plot) - mesh.Hi(:,1);
mesh.Hb_diff = mesh.Hb(:,time_to_plot) - mesh.Hb(:,1);


for i=1:size(mesh.BMB,1)
    for j=1:size(mesh.BMB,2)
        if mesh.BMB(i,j)==0
            mesh.BMB(i,j)=NaN;
        end
    end
end
    % === Grounding & ice margins ===
nE = size(mesh.E,1);
nt = length(ncread(filepath,'time'));
GL1 = ncread(filepath,'grounding_line',[1,1,1],[nE,2,1]);
GL2 = ncread(filepath,'grounding_line',[1,1,nt],[nE,2,1]);
GL2300 = ncread(filepath,'grounding_line',[1,1,time_to_plot],[nE,2,1]);
IM2 = ncread(filepath,'ice_margin',[1,1,nt],[nE,2,1]);
IM2300 = ncread(filepath,'ice_margin',[1,1,time_to_plot],[nE,2,1]);

% masks outside the ice sheet
mask_init    = mesh.Hi(:,1) > 0;
mask_current = mesh.Hi(:,time_to_plot) > 0;
mesh.Hi_diff(~mask_init)=NaN;

%% plot
xmin = -650000.0; xmax = -100000.0;
ymin = 1550000.0; ymax = 1990000.0;

fig=figure('Units','pixels','Position',[100 100 fig_width+100 fig_height],'Visible','off');
tiledlayout(nRows,nCols,"TileSpacing","compact","Padding","compact");

    % =================================
    % --- 1,1 Hi DIFF RLIS IN MESH ---
    % =================================
    ax_all(1,1) = nexttile(1); % row 3
    plot_mesh_data_a_onlyRLIS(mesh, mesh.Hi_diff(:,1), ax_all(1,1));
    %plot_mesh_data_a_RLIS(mesh, Hi_diff_in_basin_RLIS, ax_all(2,2));
    hold on;
    plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',0.8,'Color','green'); %R-LIS
    plot(GL1(:,1),GL1(:,2),'LineWidth',0.8,'LineStyle','-','Color','red');
    plot(GL2300(:,1),GL2300(:,2),'k','LineWidth',0.8);
    plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    plot(IM2300(:,1),IM2300(:,2),'k','LineWidth',0.8);
    %plot(iceshelves_MEaSUREs(37).X,iceshelves_MEaSUREs(37).Y,'LineWidth',0.8,'Color','black'); %R-LIS
    %plot(iceshelves_MEaSUREs(36).X,iceshelves_MEaSUREs(36).Y,'LineWidth',0.8,'Color','black'); % Brunt
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',true,'ncol',256);
    clim([-1000 1000]);
    ylabel('Northings (m)','FontWeight','bold');
    xlabel('Eastings (m)','FontWeight','bold');
    title('Year 2300 relative to 2000')
    %text(-6.0e5,1.95e6,'(b)','FontWeight','bold','FontSize',13);
    %text(-8e5,2.1e6,['RMSE',string(round(rmseH_in_ROI))]);
    cb=colorbar;
    cb.Label.String = 'Ice thickness difference (m)';
    cb.Label.FontSize = 12;
    text(-7.8e5,2.13e6,'(d)','FontWeight','bold','FontSize',13);
   
    print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/Hidiff_RLIS.png','-dpng','-r300');
