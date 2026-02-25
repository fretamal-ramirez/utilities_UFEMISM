%% MULTIPANEL COMPARISON PLOTS (shared colorbars)
clear all; close all; clc;

% ==== DEFINE OUTPUTS ====
outputs = { ...
    'results_ant_PD_maxphi_30_SHR_retreat_new',...
    'results_ant_PD_maxphi_30_SHR_retreat_ocndT_5e-1',...
    'results_ant_PD_maxphi_30_SHR_retreat_ocndT_2e1',...
};
titles_name = { ...
    %'max ϕ = 20', ...
    'no-shelf', ...
    'no-shelf & ocnT+0.5°C', ...
    'no-shelf & ocnT+2.0°C',...
};
tile_size = 300; % pixels for each panel
nCols = numel(outputs);
nRows = 2; % uabs, BMB, Hi_diff

fig_width  = tile_size * nCols;
fig_height = tile_size * nRows;

% ==== PATHS ====
%basepath = '/Volumes/One Touch/results_UFEMISM/tetralith_results/';
basepath = '/Users/frre9931/Desktop/tetralith_results/';
colormaps.devon = '/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt';
%plot_titles = {'Velocity (m/yr)', 'Basal melt rate (m/yr)', 'ΔIce thickness (m)'};
plot_titles = {'Northings (m)', 'Northings (m)'};
letters_for_plots = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'}; % 12 plot max for now
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
    mesh.Hb   = ncread(filepath, 'Hb');
    mesh.tfa  = ncread(filepath, 'till_friction_angle');
    [Hi_fix, maskHi0] = Hi0_to_NaN_mesh(mesh.Hi);
    mesh.Hi_diff = mesh.Hi(:,end) - mesh.Hi(:,1);
    mesh.Hb_diff = mesh.Hb(:,end) - mesh.Hb(:,1);
    mesh.uabs_diff = mesh.uabs(:,end) - mesh.uabs(:,1);
    [~, maskHi_ROI] = Hi0_to_NaN_mesh(mesh.Hi_diff);
    
    % create a mask from Hi mesh to uabs mesh
    Fmask = scatteredInterpolant(mesh.V(:,1), mesh.V(:,2), maskHi0(:,end).*maskHi_ROI,'nearest', 'none');
    mask_vel = Fmask(mesh.TriGC(:,1),mesh.TriGC(:,2));

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
    %plot_mesh_data_b_RLIS(mesh, log10(mesh.uabs(:,end)), ax_all(1,i));
    plot_mesh_data_b_RLIS(mesh, mesh.uabs_diff.*mask_vel, ax_all(1,i));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/lajolla/lajolla.cpt'...
            ,'flip',true,'ncol',256);
    clim([0 1000]);
    %cptcmap(colormaps.devon,'flip',false,'ncol',100);
    %clim([log10(1) log10(2000)]);
    title(titles_name{i},'Interpreter','none');
    text(-7.8e5,2.13e6,letters_for_plots{i},'FontWeight','bold','FontSize',13);
    %title(output_folder,'Interpreter','none');
    if i == 1
        ylabel(plot_titles{1});
    end

    % % ===================
    % % --- 2. PHI ROW ---
    % % ===================
    % ax_all(2,i) = nexttile(nCols + i); % row 2
    % plot_mesh_data_a_RLIS(mesh, mesh.tfa(:,end), ax_all(2,i));
    % hold on;
    % plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % clim([0 30]);
    % if i == 1
    %     ylabel(plot_titles{2},'FontWeight','bold');
    % end

    % ===========================
    % --- 2. HI DIFFERENCE ROW ---
    % ===========================
    ax_all(2,i) = nexttile(nCols + i); % row 2
    plot_mesh_data_a_RLIS(mesh, mesh.Hi_diff(:,1).*maskHi_ROI, ax_all(2,i));
    hold on;
    plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    plot(GL1(:,1),GL1(:,2),'LineWidth',1.0,'Color','green','linestyle','-.')
    %cptcmap('GMT_polar','flip',true,'ncol',100);
    cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',true,'ncol',256);
    clim([-1000 1000]);
    xlabel('Eastings (m)')
    text(-7.8e5,2.13e6,letters_for_plots{i+nCols},'FontWeight','bold','FontSize',13);
    if i == 1
        ylabel(plot_titles{2});
    end

    % % ===================
    % % --- 2. BMB ROW ---
    % % ===================
    % ax_all(2,i) = nexttile(nCols + i); % row 2
    % plot_mesh_data_a_RLIS(mesh, mesh.BMB(:,end), ax_all(2,i));
    % hold on;
    % plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % plot(GL1(:,1),GL1(:,2),'LineWidth',0.8,'Color','green','linestyle','-.')
    % cptcmap('GMT_polar','flip',true,'ncol',100);
    % clim([-2 2]);
    % if i == 1
    %     ylabel(plot_titles{2},'FontWeight','bold');
    % end

    % % ===========================
    % % --- 3. HI DIFFERENCE ROW ---
    % % ===========================
    % ax_all(3,i) = nexttile(2*nCols + i); % row 3
    % plot_mesh_data_a_RLIS(mesh, mesh.Hi_diff(:,1), ax_all(3,i));
    % hold on;
    % plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % plot(rock_outcrops.X,rock_outcrops.Y,'LineWidth',0.8,'color',[0.25, 0.25, 0.25]);
    % cptcmap('GMT_polar','flip',true,'ncol',100);
    % clim([-1000 1000]);
    % if i == 1
    %     ylabel(plot_titles{3},'FontWeight','bold');
    % end

    % % ===========================
    % % --- 4. Hb DIFFERENCE ROW ---
    % % ===========================
    % ax_all(4,i) = nexttile(3*nCols + i); % row 3
    % plot_mesh_data_a_RLIS(mesh, mesh.Hb_diff(:,1), ax_all(4,i));
    % hold on;
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % %cptcmap('GMT_polar','flip',true,'ncol',100);
    % clim([0 100]);
    % if i == 1
    %     ylabel(plot_titles{4},'FontWeight','bold');
    % end

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

    % % ===================
    % % --- 4. PHI ROW ---
    % % ===================
    % ax_all(4,i) = nexttile(3*nCols + i); % row 2
    % plot_mesh_data_a_RLIS(mesh, mesh.tfa(:,end), ax_all(4,i));
    % hold on;
    % plot(GL2(:,1),GL2(:,2),'k','LineWidth',0.8);
    % plot(IM2(:,1),IM2(:,2),'k','LineWidth',0.8);
    % clim([0 30]);
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
            %cb.Ticks = log10([1 10 50 100 500 1000 2000]);
            %cb.TickLabels = {'1','10','50','100','500','1000','2000'};
            cb.Label.String = 'Surface velocity difference (m/yr)';
            cb.Label.FontSize = 12;
        case 2
            %cb.Label.String = 'BMB (m/yr)';
            cb.Label.String = 'Ice thickness difference (m)';
            cb.Label.FontSize = 12;
            %cb.Label.String = 'Till friction angle (°)';
        case 3
            cb.Label.String = 'ΔHi (m)';
         case 4
            %cb.Label.String = 'Hi (m)';
            cb.Label.String = 'Till friction angle (°)';
    end
end

print(fig,'/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/multipanel/multipanel_plot.png','-dpng','-r300');