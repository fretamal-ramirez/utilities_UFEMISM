function H = plot_mesh_data_a_RLIS(mesh, d, ax)
% Plot UFEMISM mesh data (Riiser-Larsen region)
% Optional input:
%   ax — axes handle to plot into. If not provided, creates new figure.

% Plot range for Riiser-Larsen
xmid = (-878893.0+131107.0)/2;
ymid = (1137850.0+2217850.0)/2;
wx = 600e3;
wy = 800e3;

xmin = xmid - wx; xmax = xmid + wx;
ymin = ymid - wy; ymax = ymid + wy;

edgecolor = 'none';

% === Check if an axes handle was provided ===
if nargin < 3 || isempty(ax)
    % Standalone figure mode (default)
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
        'XGrid','on','YGrid','on');
end

% === Plot patch ===
if isfield(mesh, 'VVor')
    VVor_NaN = double(mesh.VVor);
    VVor_NaN(VVor_NaN == 0) = NaN;
    H.Patch = patch('Parent',H.Ax,...
        'Vertices',mesh.Vor,'Faces',VVor_NaN,...
        'FaceColor','flat','FaceVertexCData',d,...
        'EdgeColor',edgecolor);
else
    H.Patch = patch('Parent',H.Ax,...
        'Vertices',mesh.V(1:mesh.nV,:),...
        'Faces',mesh.Tri(1:mesh.nTri,:),...
        'FaceColor','interp','FaceVertexCData',d,...
        'EdgeColor',edgecolor);
end

if nargin < 3
    % only add colorbar if it's a standalone plot
    pos = get(H.Ax,'Position');
    H.Cbar = colorbar(H.Ax,'Location','eastoutside');
    set(H.Ax,'Position',pos);
end

set(H.Ax,'Units','normalized');
end
