%% MULTIPANEL COMPARISON PLOTS (shared colorbars)
% This plot has only 3 columns, 1. simulation 2. observation 3. difference
% first two columns share the colorbar, but think about what is worth to
% show, do I need to show the observation? maybe just modelled and
% difference? for example Hi is worth to show only the difference, same
% with velocities. Maybe not so with BMB, so maybe BMB model with
% parameterisation vs observed? then it will be a plot of 2x2... think
clear all; close all; clc;

% ==== DEFINE OUTPUTS ====
outputs = { ...
    'results_ant_PD_inversion_dHdt_init_R-LIS_gamma40_PMP_roughness_max20_HR',...
};
titles_name = { ...
    'HR control', ...
    'Observed', ...
    'Difference'
};
tile_size = 300; % pixels for each panel
nCols = 2;
nRows = 3; % uabs, BMB, Hi_diff

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

% ==== PATHS ====
basepath = '/Volumes/One Touch/results_UFEMISM/tetralith_results/';
%basepath = '/Users/frre9931/Desktop/tetralith_results/';
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

% Preallocate axes for later shared colorbars
ax_all = gobjects(nRows,nCols);
output_folder=outputs{1};
filepath = [basepath,output_folder,'/main_output_ANT_00001.nc'];

% Load MEaSUREs data
basins_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Basins_Antarctica_v02.shp');
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

% create a mask where there is ice, to mask the velocity later
[~, maskHi0_ufe]= Hi0_to_NaN(Hi_ufe);
[xx_measures,yy_measures]=meshgrid(x_MEaSUREs,y_MEaSUREs);
[xx_ufe,yy_ufe]=meshgrid(x_ufe,y_ufe);
% interpolate MEaSUREs to same resolution as UFEMISM output
uabs_measures_interp=interp2(xx_measures,yy_measures,uabs_MEaSUREs',xx_ufe,yy_ufe);
uabs_diff=uabs_measures_interp-(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))';

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
bmb_diff=bmb_davison_interp-abs((BMB_ufe(:,:,end).*maskHi0_ufe(:,:,end))');

%% plot
fig=figure('Units','pixels','Position',[100 100 fig_width+100 fig_height],'Visible','off');
tiledlayout(nRows,nCols,"TileSpacing","compact","Padding","compact");

    % =================================
    % --- 1,1 VELOCITY ROW IN MESH ---
    % =================================
    ax_all(1,1) = nexttile(1); % row 1
    plot_mesh_data_b_AIS(mesh, log10(mesh.uabs(:,end)), ax_all(1,1));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    cptcmap(colormaps.devon,'flip',false,'ncol',100);
    clim([log10(1) log10(2000)]);
    cb=colorbar;
    cb.Ticks = log10([1 10 50 100 500 1000 2000]);
    cb.TickLabels = {'1','10','50','100','500','1000','2000'};
    cb.Label.String = 'Velocity (m/yr)';
    cb.Label.FontSize = 12;
    ylabel('Northings (m)','FontWeight','bold');

    % =================================
    % --- 1,2 VELOCITY DIFF ---
    % =================================
    ax_all(1,2) = nexttile(2); % row 1
    rmseU=rmse(uabs_measures_interp,(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))','all','omitnan');
    contourf(x_ufe,y_ufe,uabs_diff,100,'LineColor','none');
    hold on;
    plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',2,'Color','black'); %R-LIS
    plot(basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y,'LineWidth',2,'Color','black'); % Brunt
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',false,'ncol',256);
    clim([-800 800]);
    cb=colorbar;    
    cb.Label.String = 'Velocity difference (m/yr)';
    cb.Label.FontSize = 12;
    text(-8e5,2.1e6,['RMSE',string(round(rmseU))]);
    % =================================
    % --- 2,1 Hi AIS IN MESH ---
    % =================================
    ax_all(2,1) = nexttile(3); % row 3
    rmseH=sqrt(sum(mesh.Hi_diff(:,1).^2)/length(mesh.Hi_diff(:,1)));
    plot_mesh_data_a_AIS(mesh, mesh.Hi_diff(:,1), ax_all(2,1));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-500 500]);
    ylabel('Northings (m)','FontWeight','bold');
    % =================================
    % --- 2,2 Hi RLIS IN MESH ---
    % =================================
    ax_all(2,2) = nexttile(4); % row 3
    plot_mesh_data_a_RLIS(mesh, mesh.Hi_diff(:,1), ax_all(2,2));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-500 500]);
    text(-8e5,2.1e6,['RMSE',string(round(rmseH))]);
    cb=colorbar;
    cb.Label.String = 'Ice thickness difference (m)';
    cb.Label.FontSize = 12;
    % =================================
    % --- 3,1 BMB RLIS IN MESH ---
    % =================================
    ax_all(3,1) = nexttile(5); % row 2
    plot_mesh_data_a_RLIS(mesh, mesh.BMB(:,end), ax_all(3,1));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    plot(GL1(:,1),GL1(:,2),'LineWidth',0.8,'Color','green','linestyle','-.')
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-2 2]);
    ylabel('Northings (m)','FontWeight','bold');
    xlabel('Eastings (m)','FontWeight','bold');
    cb=colorbar;
    cb.Label.String = 'Ice shelf BMB (m/yr)';
    cb.Label.FontSize = 12;
    % =================================
    % --- 3,2 BMB DIFFERENCE ---
    % =================================
    ax_all(3,2) = nexttile(6); % row 1
    %rmseBMB=rmse(uabs_measures_interp,(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))','all','omitnan');
    contourf(x_ufe,y_ufe,bmb_diff,100,'LineColor','none');
    hold on;
    plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',2,'Color','black'); %R-LIS
    plot(basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y,'LineWidth',2,'Color','black'); % Brunt
    %cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
    %        ,'flip',false,'ncol',256);
    clim([-5 5]);
    xlabel('Eastings (m)','FontWeight','bold');
    cb=colorbar;
    cb.Label.String = 'Basal melt rate difference (m/yr)';
    cb.Label.FontSize = 12;
    print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/multipanel_test.png','-dpng','-r300');
%%
for i = 1:nCols
    output_folder = outputs{i};
    filepath = fullfile(basepath, output_folder, 'main_output_ANT_00001.nc');

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

    % === Grounding & ice margins ===
    nE = size(mesh.E,1);
    nt = length(ncread(filepath,'time'));
    GL1 = ncread(filepath,'grounding_line',[1,1,1],[nE,2,1]);
    GL2 = ncread(filepath,'grounding_line',[1,1,nt],[nE,2,1]);
    IM2 = ncread(filepath,'ice_margin',[1,1,nt],[nE,2,1]);

    % ======================
    % --- 1. VELOCITY ROW ---
    % ======================
    ax_all(1,i) = nexttile(1 + (i-1)); % row 1
    plot_mesh_data_b_RLIS(mesh, log10(mesh.uabs(:,end)), ax_all(1,i));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    cptcmap(colormaps.devon,'flip',false,'ncol',100);
    clim([log10(1) log10(2000)]);
    title(titles_name{i},'Interpreter','none');
    %title(output_folder,'Interpreter','none');
    if i == 1
        ylabel(plot_titles{1},'FontWeight','bold');
    end

    % ===================
    % --- 2. BMB ROW ---
    % ===================
    ax_all(2,i) = nexttile(nCols + i); % row 2
    plot_mesh_data_a_RLIS(mesh, mesh.BMB(:,end), ax_all(2,i));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    plot(GL1(:,1),GL1(:,2),'LineWidth',0.8,'Color','green','linestyle','-.')
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-2 2]);
    if i == 1
        ylabel(plot_titles{2},'FontWeight','bold');
    end

    % ===========================
    % --- 3. HI DIFFERENCE ROW ---
    % ===========================
    ax_all(3,i) = nexttile(2*nCols + i); % row 3
    plot_mesh_data_a_RLIS(mesh, mesh.Hi_diff(:,1), ax_all(3,i));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-1000 1000]);
    if i == 1
        ylabel(plot_titles{3},'FontWeight','bold');
    end

    % % ===========================
    % % --- 4. FINAL Hi row ---
    % % ===========================
    % ax_all(4,i) = nexttile(3*nCols + i); % row 3
    % plot_mesh_data_a_RLIS(mesh, mesh.Hi(:,end).*maskHi0(:,end), ax_all(4,i));
    % hold on;
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    % %cptcmap('GMT_polar','flip',true,'ncol',100);
    % clim([0 3000]);
    % if i == 1
    %     ylabel(plot_titles{4},'FontWeight','bold');
    % end

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
            cb.Ticks = log10([1 10 50 100 500 1000 2000]);
            cb.TickLabels = {'1','10','50','100','500','1000','2000'};
            cb.Label.String = 'Velocity (m/yr)';
        case 2
            cb.Label.String = 'BMB (m/yr)';
            %cb.Label.String = 'Till friction angle (°)';
        case 3
            cb.Label.String = 'ΔHi (m)';
         case 4
            %cb.Label.String = 'Hi (m)';
            cb.Label.String = 'Till friction angle (°)';
    end
end

