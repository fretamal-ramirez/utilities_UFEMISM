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
    mesh.mask = ncread(filepath,'mask');

    % Difference between final and initial thickness
    Hi_diff = mesh.Hi(:,end)-mesh.Hi(:,1);
    [~,maskHi_ROI] = Hi0_to_NaN_mesh(Hi_diff);

    % Geographic ROI mask on the unstructured mesh vertices
    in_ROI_mesh = inpolygon( ...
        mesh.V(:,1), mesh.V(:,2), ...
        ROI_polygon.X, ROI_polygon.Y);

    % Include vertices with ice at either the initial or final time
    has_ice = mesh.Hi(:,1) > 0 | mesh.Hi(:,end) > 0;

    % Valid vertices within the ROI
    valid = in_ROI_mesh & has_ice & isfinite(Hi_diff);

    % RMS magnitude of ice thickness change (m)
    Hi_RMS(i) = sqrt(mean(Hi_diff(valid).^2));

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
    plot_mesh_data_a_RLIS(mesh, ...
        Hi_diff.*maskHi_ROI,ax_all(i),plot_font_size);
    hold(ax_all(i),'on');

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

    cptcmap( ...
        '/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt', ...
        'flip',true,'ncol',256);

    clim(ax_all(i),[-500 500]);
    title(ax_all(i),titles_name{i}, ...
        'Interpreter','none','FontSize',10);

    text(ax_all(i),-7.8e5,2.13e6,letters_for_plots{i}, ...
        'FontWeight','bold','FontSize',10);
    text(ax_all(i), -7.8e5, 2.00e6, ...
        sprintf('E = %.0f m', Hi_RMS(i)), ...
        'FontSize',8,'FontWeight','bold');

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
cb.Label.String = 'Ice thickness difference (m)';
cb.Ticks = -500:100:500;

print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/multipanel_plotH_init.pdf','-dpdf','-vector');