%% MULTIPANEL VIDEO: Velocity + ΔHi every 100 yrs relative to 2000
clear all; close all; clc;

% ==== DEFINE OUTPUTS ====
outputs = { ...
    %'results_ant_PD_maxphi_30_SHR_ctrl2500',...
    'results_ant_PD_maxphi_30_SHR_ctrl2500_new',...
    'results_ant_PD_maxphi_Hb-2000to-250m_SHR_ctrl2500',...
    %'results_ant_PD_maxphi_30_SHR_retreat',...
    'results_ant_PD_maxphi_30_SHR_retreat_new',...
    'results_ant_PD_maxphi_Hb-2000to-250m_SHR_retreat',...
};

titles_name = { ...
    %'control', ...
    'control', ...
    'control Hb-2000to-250m', ...
    %'no-shelf', ...
    'no-shelf', ...
    'no-shelf Hb-2000to-250m',...
};

% ==== TIME SETTINGS ====
time_indices = 2:1:11;         % it was 3:2:11
years        = 2050:50:2500;   % it was 2100:100:2500

% ==== FIGURE SETTINGS ====
tile_size = 300;
nCols = numel(outputs);
nRows = 2;                     % velocity + Hi
fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

% ==== PATHS ====
basepath = '/Users/frre9931/Desktop/tetralith_results/';

path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

rock_outcrops = shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ADD_RockOutcrops_RLIS.shp');

% ==== VIDEO SETUP ====
video_name = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/Velocity_Hi_diff_100yr.mp4';
v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 1;
open(v);

% ==== CREATE FIGURE ====
fig = figure('Visible','off','Units','pixels','Position',[100 100 fig_width fig_height]);
tiledlayout(nRows,nCols,"TileSpacing","compact","Padding","compact");

% ==========================================================
% ===================== TIME LOOP ==========================
% ==========================================================

for t = 1:length(time_indices)

    tid = time_indices(t);

    ax_all = gobjects(nRows,nCols);

    for i = 1:nCols

        output_folder = outputs{i};
        filepath = fullfile(basepath, output_folder, 'main_output_ANT_00001.nc');

        % --- Load mesh and variables ---
        mesh = read_mesh_from_file(filepath);
        mesh.Hi   = ncread(filepath,'Hi');
        mesh.uabs = ncread(filepath,'uabs_surf');

        % --- Differences relative to 2000 ---
        %Hi_diff   = mesh.Hi(:,tid)   - mesh.Hi(:,1);
        %uabs_diff = mesh.uabs(:,tid) - mesh.uabs(:,1);

        %masks
        % --- Compute differences ---
        Hi_initial = mesh.Hi(:,1);
        Hi_current = mesh.Hi(:,tid);

        Hi_diff = Hi_current - Hi_initial;

        % --- Build masks ---
        mask_initial = Hi_initial > 0;
        mask_current = Hi_current > 0;
        mask_ROI     = abs(Hi_diff) < 1e-4; % No change imposed outside ROI

        mask_combined = mask_current & ~mask_ROI;

        % --- Apply mask ---
        Hi_diff(~mask_initial) = NaN;

        uabs_initial = mesh.uabs(:,1);
        uabs_current = mesh.uabs(:,tid);

        uabs_diff = uabs_current - uabs_initial;

        % Project node mask to triangle centers
        Fmask = scatteredInterpolant(mesh.V(:,1), mesh.V(:,2), ...
                             double(mask_combined), ...
                             'nearest','none');

        mask_vel = Fmask(mesh.TriGC(:,1),mesh.TriGC(:,2));

        uabs_diff(mask_vel==0) = NaN;
        Hi_diff(mask_ROI) = NaN;

        nE = size(mesh.E,1);
        GLt = ncread(filepath,'grounding_line',[1,1,tid],[nE,2,1]);
        IMt = ncread(filepath,'ice_margin',[1,1,tid],[nE,2,1]);

        % ===================================================
        % ================= VELOCITY ROW ====================
        % ===================================================

        ax_all(1,i) = nexttile(i);
        cla(ax_all(1,i))

        plot_mesh_data_b_RLIS(mesh, uabs_diff .* mask_vel, ax_all(1,i));
        hold on;

        plot(GLt(:,1),GLt(:,2),'k','LineWidth',0.8);
        plot(IMt(:,1),IMt(:,2),'k','LineWidth',0.8);
        plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25 0.25 0.25]);

        cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/lajolla/lajolla.cpt',...
                'flip',true,'ncol',256);

        clim([0 500]);

        title(titles_name{i},'Interpreter','none');

        if i == 1
            ylabel('Surface velocity difference (m/yr)');
        end

        % ===================================================
        % =================== HI ROW ========================
        % ===================================================

        ax_all(2,i) = nexttile(nCols + i);
        cla(ax_all(2,i))

        plot_mesh_data_a_RLIS(mesh, Hi_diff, ax_all(2,i));
        hold on;

        plot(GLt(:,1),GLt(:,2),'k','LineWidth',0.8);
        plot(IMt(:,1),IMt(:,2),'k','LineWidth',0.8);
        plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25 0.25 0.25]);

        cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt',...
                'flip',true,'ncol',256);

        clim([-500 500]);

        if i == 1
            ylabel('Ice thickness difference (m)');
        end

    end

    % ===================================================
    % ============ SHARED COLORBARS ====================
    % ===================================================

    for r = 1:nRows
        row_axes = ax_all(r,:);
        clim_values = cell2mat(get(row_axes,'CLim'));
        minCL = min(clim_values(:,1));
        maxCL = max(clim_values(:,2));
        set(row_axes,'CLim',[minCL maxCL]);

        cb = colorbar(row_axes(end),'eastoutside');

        switch r
            case 1
                cb.Label.String = 'Surface velocity difference (m/yr)';
            case 2
                cb.Label.String = 'Ice thickness difference (m)';
        end
        cb.Label.FontSize = 12;
    end

    sgtitle(['Relative to year 2000 — Year ' num2str(years(t))],...
            'FontSize',14,'FontWeight','bold');

    drawnow

    frame = getframe(fig);
    writeVideo(v,frame);

end

close(v);

disp('Video successfully created.')