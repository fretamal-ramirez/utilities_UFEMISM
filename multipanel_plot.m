%% MULTIPANEL COMPARISON PLOTS (shared colorbars)
clear all; clc;

% ==== DEFINE OUTPUTS ====
outputs = { ...
    'results_ant_PD_smooth_Hb_retreat_mask_NMP', ...
    'results_ant_PD_smooth_Hb_retreat_mask_FCMP', ...
    'results_ant_PD_smooth_Hb_retreat_mask_PMP' ...
};
tile_size = 300; % pixels for each panel
nCols = numel(outputs);
nRows = 3; % uabs, BMB, Hi_diff

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

% ==== PATHS ====
%basepath = '/Volumes/One Touch/results_UFEMISM/tetralith_results/';
basepath = '/Users/frre9931/Desktop/tetralith_results/';
colormaps.devon = '/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt';
plot_titles = {'Velocity (m/yr)', 'Basal melt rate (m/yr)', 'ΔIce thickness (m)'};

% ==== Add extra functions ====
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
%path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
%path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

fig=figure('Units','pixels','Position',[100 100 fig_width fig_height],'Visible','off');
tiledlayout(nRows,nCols,"TileSpacing","compact","Padding","compact");

% Preallocate axes for later shared colorbars
ax_all = gobjects(nRows,nCols);

for i = 1:nCols
    output_folder = outputs{i};
    filepath = fullfile(basepath, output_folder, 'main_output_ANT_00001.nc');

    % === Load mesh ===
    mesh = read_mesh_from_file(filepath);
    mesh.uabs = ncread(filepath,'uabs_surf');
    mesh.BMB  = ncread(filepath,'BMB');
    mesh.Hi   = ncread(filepath,'Hi');
    mesh.mask = ncread(filepath,'mask');
    [Hi_fix, maskHi0] = Hi0_to_NaN_mesh(mesh.Hi);
    mesh.Hi_diff = mesh.Hi(:,end) - mesh.Hi(:,1);

    % === Grounding & ice margins ===
    nE = size(mesh.E,1);
    nt = length(ncread(filepath,'time'));
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
    title(output_folder,'Interpreter','none');
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
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-250 250]);
    if i == 1
        ylabel(plot_titles{3},'FontWeight','bold');
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
            cb.Ticks = log10([1 10 50 100 500 1000 2000]);
            cb.TickLabels = {'1','10','50','100','500','1000','2000'};
            cb.Label.String = 'Velocity (m/yr)';
        case 2
            cb.Label.String = 'BMB (m/yr)';
        case 3
            cb.Label.String = 'ΔHi (m)';
    end
end

print(fig,'multipanel_plot.png','-dpng','-r300');