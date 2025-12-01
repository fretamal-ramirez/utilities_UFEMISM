function H = plot_mesh_data_b_AIS( mesh, d, ax)
% Optional input:
%   ax — axes handle to plot into. If not provided, creates new figure.

edgecolor = 'none';
% edgecolor = 'k';

% === Check if an axes handle was provided ===
if nargin < 3 || isempty(ax)

    % Axes and figure size
    xw = mesh.xmax - mesh.xmin;
    yw = mesh.ymax - mesh.ymin;
    if xw >= yw
    wa = 800;
     ha = wa * yw / xw;
    else
    ha = 800;
    wa = ha * xw / yw;
    end
    wf = 25 + wa + 100;
    hf = 25 + ha + 50;

    H.Fig = figure('position',[200,200,wf,hf],'color','w');
    H.Ax  = axes('parent',H.Fig,'units','pixels','position',[25,25,wa,ha],'fontsize',24,...
    'xtick',[],'ytick',[],'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
else
    % Multi-panel mode
    H.Fig = ancestor(ax,'figure');
    H.Ax = ax;
    set(H.Ax,'XLim',[mesh.xmin,mesh.xmax],'YLim',[mesh.ymin,mesh.ymax],...
        'XGrid','on','YGrid','on');
end
H.Patch = patch('vertices',mesh.V( 1:mesh.nV,:),'faces',mesh.Tri( 1:mesh.nTri,:),...
     'facecolor','flat','facevertexcdata',d,'edgecolor',edgecolor);

% === Colorbar only for standalone mode ===
if nargin < 3
pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);
end

set( H.Ax,'units','normalized');

end