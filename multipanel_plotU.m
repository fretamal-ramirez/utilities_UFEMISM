%% MULTIPANEL DIFFERENCE SURFACE VELOCITY

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

pageWidth  = 179;   % mm
pageHeight = 254;   % mm

tile_size = 300; % pixels for each panel
nCols = 3;
nRows = 2; % uabs, BMB, Hi_diff

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

plot_titles = {'Northings (m)', 'Northings (m)'};
plot_font_size = 8;
letters_for_plots = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'}; % 12 plot max for now

% ==== PATHS ====
basepath = '/Users/frre9931/Desktop/arrhenius_results/';

colormaps.devon = '/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt';

% ==== Add extra functions ====
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

% ===== Data to complement plots ====
rock_outcrops=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ADD_RockOutcrops_RLIS.shp');

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

nPanels = numel(outputs);

assert(nPanels <= nRows*nCols, ...
    'There are more simulations than available tiles.');
assert(numel(titles_name) == nPanels, ...
    'Each simulation needs a title.');

% create figure
fig = figure('Visible','off', ...
   'Units','centimeters', ...
   'Position',[2 2 17.9 25.4]);

set(fig,'PaperUnits','centimeters');
set(fig,'PaperSize',[17.9 25.4]);
set(fig,'PaperPosition',[0 0 17.9 25.4]);

tl = tiledlayout(nRows,nCols, ...
    'TileSpacing','compact','Padding','compact');

ax_all = gobjects(nPanels,1);

for i = 1:nPanels

    % Automatically determine row, column, and tile
    r = ceil(i/nCols);
    c = mod(i-1,nCols)+1;
    tile = (r-1)*nCols+c;

    output_folder = outputs{i};
    filepath = [basepath,output_folder,'/main_output_ANT_00001.nc'];

% === Load mesh ===
mesh = read_mesh_from_file(filepath);
mesh.uabs = ncread(filepath,'uabs_surf');
mesh.Hi   = ncread(filepath,'Hi');
mesh.mask = ncread(filepath,'mask');
[Hi_fix, maskHi0] = Hi0_to_NaN_mesh(mesh.Hi);
mesh.Hi_diff = mesh.Hi(:,end) - mesh.Hi(:,1);

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

for k=1:size(in_ROI,1)
    for j=1:size(in_ROI,2)
        if diff_uabs_in_ROI(k,j) == 0
            diff_uabs_in_ROI(k,j)=NaN;
        elseif diff_uabs_in_basin_RLIS(k,j) == 0
            diff_uabs_in_basin_RLIS(k,j) = NaN;
        elseif diff_uabs_in_basin_BIS(k,j) == 0
            diff_uabs_in_basin_BIS(k,j) = NaN;
        elseif diff_uabs_in_shelf_RLIS(k,j) == 0
            diff_uabs_in_shelf_RLIS(k,j) = NaN;
        elseif diff_uabs_in_shelf_BIS(k,j) == 0
            diff_uabs_in_shelf_BIS(k,j) = NaN;
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

for k=1:length(mesh.Hi_diff_in_ROI)
    if mesh.Hi_diff_in_ROI(k) == 0
        mesh.Hi_diff_in_ROI(k) = NaN;
    elseif Hi_diff_in_basin_RLIS(k) == 0
        Hi_diff_in_basin_RLIS(k) = NaN;
    elseif Hi_diff_in_basin_BIS(k) == 0
        Hi_diff_in_basin_BIS(k) = NaN;
    elseif Hi_diff_in_shelf_RLIS(k) == 0
        Hi_diff_in_shelf_RLIS(k) = NaN;
    elseif Hi_diff_in_shelf_BIS(k) == 0
        Hi_diff_in_shelf_BIS(k) = NaN;
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

%% plot
    % =================================
    % --- 1,2 VELOCITY DIFF ---
    % =================================
    ax_all(i) = nexttile(tile); % row 1
    %contourf(x_ufe,y_ufe,uabs_diff,50,'LineColor','none','EdgeColor','none');
    %contourf(x_ufe,y_ufe,diff_uabs_in_ROI,100,'LineColor','none');
    h=imagesc(ax_all(i),x_ufe,y_ufe,uabs_diff);
    hold(ax_all(i),'on');
    plot(ax_all(i),basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',0.8,'Color','black'); %R-LIS it was in 2.0 before
    plot(ax_all(i),basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y,'LineWidth',0.8,'Color','black'); % Brunt
    plot(ax_all(i),iceshelves_MEaSUREs(37).X,iceshelves_MEaSUREs(37).Y,'LineWidth',0.8,'Color','black'); %R-LIS
    plot(ax_all(i),iceshelves_MEaSUREs(36).X,iceshelves_MEaSUREs(36).Y,'LineWidth',0.8,'Color','black'); % Brunt
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',false,'ncol',256);
    clim([-800 800]);
    set(ax_all(i),'XLim',[x_ufe_min,x_ufe_max],'YLim',[y_ufe_min,y_ufe_max],...
        'XGrid','on','YGrid','on');
    box off
    set(h,'AlphaData', ~isnan(uabs_diff)); axis xy; % make NaN transparent and flip image
    
    title(ax_all(i),titles_name{i}, ...
        'Interpreter','none','FontSize',10);

    text(ax_all(i),-7.8e5,2.13e6,letters_for_plots{i}, ...
        'FontWeight','bold','FontSize',10);
    text(ax_all(i), -7.8e5, 2.00e6, ...
        sprintf('E = %.0f m', round(rmseU_in_ROI)), ...
        'FontSize',8,'FontWeight','bold');

    if r == nRows
        xlabel(ax_all(i),'Eastings (m)', ...
            'FontWeight','bold','FontSize',10);
    end
    if c == 1
        ylabel(ax_all(i),'Northings (m)', ...
            'FontWeight','bold','FontSize',10);
    end
    %text(-8e5,2.1e6,['RMSE',string(round(rmseU_in_ROI))]);
end

cb = colorbar;
cb.Layout.Tile = 'east';  % Place it to the right of all panels
cb.FontSize = 9;
cb.LineWidth = 0.8;

% One shared colorbar for the complete figure
%cb = colorbar(ax_all(end),'eastoutside');
%cb.Label.FontSize = 10;
cb.Label.String = 'Surface velocity difference (m a^{-1})';
%cb.Ticks = -500:100:500;

print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/multipanel_plot_uabs_init.pdf','-dpdf','-vector');