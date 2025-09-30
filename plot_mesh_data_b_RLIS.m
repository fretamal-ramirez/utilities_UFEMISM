function H = plot_mesh_data_b_RLIS( mesh, d)

% Plot range for Riiser-Larsen
xmid = (-878893.0+131107.0)/2;
ymid = (1137850.0+2217850.0)/2;
wx    = 600e3;
wy    = 800e3;

xmin = xmid - wx;
xmax = xmid + wx;
ymin = ymid - wy;
ymax = ymid + wy;

wa = 750;
ha = 750;

margins_hor = [125,125];
margins_ver = [125,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);

edgecolor = 'none';
% edgecolor = 'k';

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[xmin,xmax],'ylim',[ymin,ymax],'fontsize',24,'xgrid','on','ygrid','on');
H.Patch = patch('vertices',mesh.V( 1:mesh.nV,:),'faces',mesh.Tri( 1:mesh.nTri,:),...
  'facecolor','flat','facevertexcdata',d,'edgecolor',edgecolor);

pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);

set( H.Ax,'units','normalized');

end