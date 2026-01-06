function H = plot_mesh_data_b_RLIS(mesh, d, ax)
% Plot UFEMISM mesh data (Riiser-Larsen region)
% Optional input:
%   ax — axes handle to plot into. If not provided, creates new figure.

% Plot range for Riiser-Larsen
xmid = (-878893.0 + 131107.0)/2;
ymid = (1137850.0 + 2217850.0)/2;
wx = 600e3;
wy = 800e3;

xmin = xmid - wx; xmax = xmid + wx;
ymin = ymid - wy; ymax = ymid + wy;
xmin = -878893.0; xmax = 131107.0;
ymin = 1137850.0; ymax = 2217850.0;
edgecolor = 'none';

% === Check if an axes handle was provided ===
if nargin < 3 || isempty(ax)
    % Standalone figure mode
    wa = 750; ha = 750;
    margins_hor = [125,125];
    margins_ver = [125,50];
    wf = sum(margins_hor) + wa;
    hf = sum(margins_ver) + ha;

    H.Fig = figure('Position',[200,200,wf,hf],'Color','w');
    H.Ax = axes('Parent',H.Fig,'Units','pixels',...
        'Position',[margins_hor(1),margins_ver(1),wa,ha],...
        'XLim',[xmin,xmax],'YLim',[ymin,ymax],...
        'FontSize',24,'XGrid','on','YGrid','on');
else
    % Multi-panel mode
    H.Fig = ancestor(ax,'figure');
    H.Ax = ax;
    set(H.Ax,'XLim',[xmin,xmax],'YLim',[ymin,ymax],...
        'XGrid','on','YGrid','on','Layer','top');
end

% === Plot patch ===
H.Patch = patch('Parent',H.Ax,...
    'Vertices',mesh.V(1:mesh.nV,:),...
    'Faces',mesh.Tri(1:mesh.nTri,:),...
    'FaceColor','flat','FaceVertexCData',d,...
    'EdgeColor',edgecolor);

% === Colorbar only for standalone mode ===
if nargin < 3
    pos = get(H.Ax,'Position');
    H.Cbar = colorbar(H.Ax,'Location','eastoutside');
    set(H.Ax,'Position',pos);
end

set(H.Ax,'Units','normalized');
end
