%% MULTIPANEL COMPARISON PLOTS of ice thickness (shared colorbars)
clear all; close all; clc;

% ==== DEFINE OUTPUTS ====
outputs = { ...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_M11_Hb-2000to-250m_SHR_trans',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_max30_SHR_trans',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_max40_SHR',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_M11_Hb-2000to-250m_Budd',...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma99_PMP_roughness_max30_Budd',...
};
titles_name = { ...
    'BCAP',...
    'STD',...
    'STD mphi40',...
    'BCAP SL2',...
    'STD SL2',...
};
tile_size = 300; % pixels for each panel
%nCols = numel(outputs);
nCols = 3;
nRows = 2; % Hi_diff

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

% ==== PATHS ====
%basepath = '/Volumes/One Touch/results_UFEMISM/tetralith_results/outputs';
%basepath = '/Users/frre9931/Desktop/tetralith_results/';
basepath = '/Users/frre9931/Desktop/arrhenius_results/';
colormaps.devon = '/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt';
%plot_titles = {'Velocity (m/yr)', 'Basal melt rate (m/yr)', 'ΔIce thickness (m)'};
plot_titles = {'Northings (m)', 'Northings (m)'};
plot_font_size = 8;
letters_for_plots = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'}; % 12 plot max for now
%plot_titles = {'Velocity (m/yr)', 'Basal melt rate (m/yr)', 'ΔIce thickness (m)', 'Till friction angle (°)'};
%plot_titles = {'Velocity (m/yr)', 'Till friction angle (°)', 'ΔIce thickness (m)', 'Final ice thickness (m)'};

% palette for simple colors
palette_jfly = [
    102 102 102;
    230 159   0;
    86  180 233;
    0   158 115;
    240 226 66;
    0   114 178;
    213  94   0;
    204 121 167;
] / 255;

% ==== Add extra functions ====
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
%path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
%path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

% ===== Data to complement plots ====
ROI_polygon=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ROI_for_R-LIS.shp');
rock_outcrops=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ADD_RockOutcrops_RLIS.shp');

fig = figure('Visible','off', ...
   'Units','centimeters', ...
   'Position',[2 2 17.9 25.4]);

%fig = figure('Visible','off');

set(fig,'PaperUnits','centimeters');
set(fig,'PaperSize',[17.9 25.4]);
set(fig,'PaperPosition',[0 0 17.9 25.4]);

tiledlayout(nRows,nCols,"TileSpacing","compact","Padding","compact");
basins_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Basins_Antarctica_v02.shp');

% Number of panels

nPanels = numel(outputs);

assert(nPanels <= nRows*nCols, ...
    'There are more simulations than available tiles.');
assert(numel(titles_name) == nPanels, ...
    'Each simulation needs a title.');

fig = figure('Visible','off');
tl = tiledlayout(nRows,nCols, ...
    'TileSpacing','compact','Padding','compact');

ax_all = gobjects(nPanels,1);

for i = 1:nPanels

    % Automatically determine row, column, and tile
    r = ceil(i/nCols);
    c = mod(i-1,nCols)+1;
    tile = (r-1)*nCols+c;

    output_folder = outputs{i};
    filepath = fullfile(basepath,output_folder, ...
        'main_output_ANT_00001.nc');

    % Load mesh and ice thickness
    mesh = read_mesh_from_file(filepath);
    mesh.Hi = ncread(filepath,'Hi');
    
    % load all components of strain
    du_dx = ncread(filepath,'du_dx_3D');
    du_dy = ncread(filepath,'du_dy_3D');
    du_dz = ncread(filepath,'du_dz_3D');
    dv_dx = ncread(filepath,'dv_dx_3D');
    dv_dy = ncread(filepath,'dv_dy_3D');
    dv_dz = ncread(filepath,'dv_dz_3D');
    dw_dx = ncread(filepath,'dw_dx_3D');
    dw_dy = ncread(filepath,'dw_dy_3D');
    dw_dz = ncread(filepath,'dw_dz_3D');

    % diagonal tensor for surface and last time
    exx = du_dx(:,:,end);
    eyy = dv_dy(:,:,end);
    ezz = - (exx + eyy);
    %ezz = dw_dz(:,1,end);
    % 
    exy = 0.5 * (du_dy(:,:,end) + dv_dx(:,:,end));
    exz = 0.5 * (du_dz(:,:,end) + dw_dx(:,:,end));
    eyz = 0.5 * (dv_dz(:,:,end) + dw_dy(:,:,end));

    eps_eff = sqrt(0.5 .*(exx.^2 + eyy.^2 + ezz.^2) + exy.^2 + exz.^2 + eyz.^2 );
    %eps_eff = mean(eps_eff,2);
    eps_eff = eps_eff(:,end);%

    mesh.mask = ncread(filepath,'mask');
    % vertically averaged velocities from DIVA
    u_vav = ncread(filepath,'u_vav');
    v_vav = ncread(filepath,'v_vav');

    % interpolate velocities to voronoi mesh
    u_vav_to_vi = scatteredInterpolant(mesh.TriGC(:,1),mesh.TriGC(:,2), u_vav,'nearest', 'none');
    v_vav_to_vi = scatteredInterpolant(mesh.TriGC(:,1),mesh.TriGC(:,2), v_vav,'nearest', 'none');
    u_vav_vi = u_vav_to_vi(mesh.V(:,1), mesh.V(:,2));
    v_vav_vi = v_vav_to_vi(mesh.V(:,1), mesh.V(:,2));

    %calculate grounded ice using the mesh mask
    %grounded_idx=find(mesh.mask(:,end)~=3 | mesh.mask(:,end)~=5 | mesh.mask(:,end)~=7);
    % C%type_icefree_land                        = 1
    % C%type_icefree_ocean                       = 2
    % C%type_grounded_ice                        = 3
    % C%type_floating_ice                        = 4
    % C%type_groundingline_gr                    = 5
    % C%type_groundingline_fl                    = 6
    % C%type_calvingfront_gr                     = 7
    % C%type_calvingfront_fl                     = 8
    % C%type_margin                              = 9
    % C%type_coastline                           = 10
    % Identify grounded ice (including grounding line and grounded calving front)
    grounded_idx = ismember(mesh.mask(:,end), [3, 5, 7]);
    
    % calculate the horizontal shear strain for the last time, calculated
    % as 1/2 * (du/dy + dv/dx), for the last vertical level, which does not
    % matter as it is a 2D variable (same in every vertical level)
    horizontal_shear_strain = 0.5 * (du_dy(:,end,end) + dv_dx(:,end,end));
    
    % mask out floating ice
    horizontal_shear_strain(~grounded_idx)=NaN;
    u_vav_vi(~grounded_idx,end)=NaN;
    v_vav_vi(~grounded_idx,end)=NaN;
    eps_eff(~grounded_idx) = NaN;

    % Difference between final and initial thickness
    [~,maskHi_ROI] = Hi0_to_NaN_mesh(mesh.Hi(:,end));

    % Read grounding line and ice margin
    nE = size(mesh.E,1);
    nt = numel(ncread(filepath,'time'));

    GL1 = ncread(filepath,'grounding_line', ...
        [1,1,1],[nE,2,1]);
    GL2 = ncread(filepath,'grounding_line', ...
        [1,1,nt],[nE,2,1]);
    IM2 = ncread(filepath,'ice_margin', ...
        [1,1,nt],[nE,2,1]);

    % Plot in automatically selected tile
    ax_all(i) = nexttile(tile);
    %plot_mesh_data_a_Veststraumen(mesh,horizontal_shear_strain,ax_all(i),plot_font_size);
    plot_mesh_data_a_Veststraumen(mesh,eps_eff,ax_all(i),plot_font_size);
    hold(ax_all(i),'on');

    % QUIVER PLOT
    % quiver plot in Tri mesh
    %quiver(mesh.TriGC(1:10:end,1),mesh.TriGC(1:10:end,2),u_vav(1:10:end,end),v_vav(1:10:end,end),...
    %    'AutoScaleFactor',1.5,'LineWidth',0.3);
    % quiver plot in voronoi mesh
    quiver(mesh.V(1:5:end,1),mesh.V(1:5:end,2),u_vav_vi(1:5:end,end),v_vav_vi(1:5:end,end),...
        'AutoScaleFactor',3,'LineWidth',0.3);

    %plot(ax_all(i),GL2(:,1),GL2(:,2),'LineWidth',0.7,'Color',palette_jfly(4,:));
    plot(ax_all(i),GL2(:,1),GL2(:,2),'LineWidth',0.7,'Color','black');
    plot(ax_all(i),IM2(:,1),IM2(:,2), ...
        'k','LineWidth',0.7);

    plot(ax_all(i),basins_MEaSUREs(3).X, ...
        basins_MEaSUREs(3).Y, ...
        'LineWidth',0.7,'Color',palette_jfly(8,:));
    plot(ax_all(i),basins_MEaSUREs(4).X, ...
        basins_MEaSUREs(4).Y, ...
        'LineWidth',0.7,'Color',palette_jfly(8,:));

    plot(ax_all(i),rock_outcrops.X,rock_outcrops.Y, ...
        'LineWidth',0.7,'Color',[0.25 0.25 0.25], ...
        'LineStyle',':');

    % palette from red to blue
    % cptcmap( ...
    %     '/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt', ...
    %     'flip',true,'ncol',256);
    % clim(ax_all(i),[-0.03 0.03]);
    
    % palette from yellow to black
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/lajolla/lajolla.cpt'...
            ,'flip',true,'ncol',256);
    clim([0 0.1]);

    title(ax_all(i),titles_name{i}, ...
        'Interpreter','none','FontSize',10);

    text(ax_all(i),-7.8e5,2.13e6,letters_for_plots{i}, ...
        'FontWeight','bold','FontSize',10);

    if r == nRows
        xlabel(ax_all(i),'Eastings (m)', ...
            'FontWeight','bold','FontSize',10);
    end
    if c == 1
        ylabel(ax_all(i),'Northings (m)', ...
            'FontWeight','bold','FontSize',10);
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';  % Place it to the right of all panels
cb.FontSize = 9;
cb.LineWidth = 0.8;

% One shared colorbar for the complete figure
%cb = colorbar(ax_all(end),'eastoutside');
%cb.Label.FontSize = 10;
cb.Label.String = 'Effective strain rate (a^{-1})';
%cb.Ticks = -500:100:500;

print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/multipanel_plot_shear_init.pdf','-dpdf','-vector');