% code to load and read variables from UFEMISM simulations
clear all; clc;

%========= PATH TO OUTPUT UFEMISM DIRECTORY ==========
output_folder = 'results_ant_PD_inversion_dHdt_init_R-LIS_gamma80_PMP_roughness_max30_SHR';
%ufe_folder_path=['/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/', output_folder];
%ufe_folder_path=['/Users/frre9931/Desktop/UFEMISM2.0_porting/', output_folder];
ufe_folder_path=['/Users/frre9931/Desktop/tetralith_results/', output_folder];
%ufe_folder_path=['/Volumes/One Touch/results_UFEMISM/tetralith_results/', output_folder];
allow_plot_mesh = true; % if we want to plot the mesh
allow_save_plots = false;
path_save = '/Users/frre9931/Documents/PhD/ANT_UFEMISM/plots_ant/Riiser-Larsen/';
allow_mesh_update = false; % if remeshing is allowed in simulation
% number pointing the file with updated mesh
number_mesh ='2'; % i.e. 2 for main_output_ANT_00002.nc
plot_Voronoi=true;
%========= END OF CONFIGURATION ======================

% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));
path(path,genpath('/Users/frre9931/Documents/PhD/m_map'));
path(path,genpath('/Users/frre9931/Documents/PhD/Antarctic-Mapping-Tools-main'));
path(path,genpath('/Users/frre9931/Documents/PhD/cptcmap-pkg/cptcmap'));

% load the coastline from MEaSUREs
coast_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Coastline_Antarctica_v02.shp');
%function from Antarctic Maping Tools to project from x,y to lon,lat
[coast_lat,coast_lon]=ps2ll(coast_MEaSUREs.X,coast_MEaSUREs.Y);
% load the basins from MEaSUREs
basins_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Basins_Antarctica_v02.shp');
%u_MEaSUREs=ncread('/Users/frre9931/Documents/PhD/MEaSUREs/antarctica_ice_velocity_450m_v2.nc','VX');
%v_MEaSUREs=ncread('/Users/frre9931/Documents/PhD/MEaSUREs/antarctica_ice_velocity_450m_v2.nc','VY');
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

%% load the mesh

mesh_path_first= [ufe_folder_path, '/main_output_ANT_00001.nc']; %initial state

% read_mesh_from_file + CL, GL and CF
mesh_first=read_mesh_from_file(mesh_path_first);

time_slice_init=1 ; % first time
time_slice1=ncread(mesh_path_first,'time'); % load the time of the mesh

CL=ncread(mesh_path_first,'coastline',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
GL=ncread(mesh_path_first,'grounding_line',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
CF=ncread(mesh_path_first,'calving_front',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
IM=ncread(mesh_path_first,'ice_margin',[1,1,1], [size(mesh_first.E,1),2,time_slice_init]);
%GIC=ncread(mesh_path_first,'grounded_ice_contour',[1,1,1],[11788,2,time_slice_init]); % not really useful

%add variables of the mesh during the "last" time,
CL2=ncread(mesh_path_first,'coastline',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
GL2=ncread(mesh_path_first,'grounding_line',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
CF2=ncread(mesh_path_first,'calving_front',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);
IM2=ncread(mesh_path_first,'ice_margin',[1,1,length(time_slice1)], [size(mesh_first.E,1),2,1]);

% CL3=ncread(mesh_path_first,'coastline',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
% GL3=ncread(mesh_path_first,'grounding_line',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
% CF3=ncread(mesh_path_first,'calving_front',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);
% IM3=ncread(mesh_path_first,'ice_margin',[1,1,length(time_slice1)-10], [size(mesh_first.E,1),2,1]);


if allow_mesh_update
    % path to the new mesh, for now just changing the number here
    mesh_path_update = [ufe_folder_path, '/main_output_ANT_0000',number_mesh,'.nc']; 
    mesh_updated=read_mesh_from_file(mesh_path_update);
    time_slice2=ncread(mesh_path_update,'time');
    CL_mesh2=ncread(mesh_path_update,'coastline',[1,1,1], [11788,2,length(time_slice2)]);
    GL_mesh2=ncread(mesh_path_update,'grounding_line',[1,1,1], [11788,2,length(time_slice2)]);
    CF_mesh2=ncread(mesh_path_update,'calving_front',[1,1,1], [11788,2,length(time_slice2)]);
end
if allow_plot_mesh
    plot_mesh(mesh_first);
    hold on
    plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
    hold off
    if allow_save_plots
        %print([path_save,output_folder,'_mesh_1'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh(mesh_updated);
        hold on
        plot(CF_mesh2(:,1),CF_mesh2(:,2),'LineWidth',2,'Color','red');
        plot(CL_mesh2(:,1),CL_mesh2(:,2),'LineWidth',2,'Color','blue');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','green');
        hold off
        if allow_save_plots
            %print([path_save,output_folder,'_mesh_',number_mesh],'-dpng','-r300')
        end
    else
        plot_mesh(mesh_first);
        hold on
        plot(CF2(:,1),CF2(:,2),'LineWidth',2,'Color','red');
        plot(CL2(:,1),CL2(:,2),'LineWidth',2,'Color','blue');
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','green');
        hold off
    end    
end
%% add plots with Voronois cells
mesh_first.SMB            = ncread( mesh_path_first,'SMB');
mesh_first.Hi            = ncread( mesh_path_first,'Hi');
mesh_first.T2m            = ncread( mesh_path_first,'T2m');
mesh_first.Precip            = ncread( mesh_path_first,'Precip');
mesh_first.uabs              = ncread( mesh_path_first,'uabs_surf');
mesh_first.BMB               = ncread( mesh_path_first,'BMB');
mesh_first.Hb                 = ncread( mesh_path_first,'Hb');
[Hi_fix_mesh, maskHi0_mesh]= Hi0_to_NaN_mesh(mesh_first.Hi);
mesh_first.mask = ncread( mesh_path_first, 'mask');

if allow_mesh_update
    % path to the new mesh, for now just changing the number here
    mesh_path_update = [ufe_folder_path, '/main_output_ANT_0000',number_mesh,'.nc']; 
    mesh_updated=read_mesh_from_file(mesh_path_update);
    time_slice2=ncread(mesh_path_update,'time');
    CL_mesh2=ncread(mesh_path_update,'coastline',[1,1,1], [11788,2,length(time_slice2)]);
    GL_mesh2=ncread(mesh_path_update,'grounding_line',[1,1,1], [11788,2,length(time_slice2)]);
    CF_mesh2=ncread(mesh_path_update,'calving_front',[1,1,1], [11788,2,length(time_slice2)]);
    IM_mesh2=ncread(mesh_path_update,'ice_margin',[1,1,1], [11788,2,length(time_slice2)]);
    mesh_updated.SMB    = ncread(mesh_path_update, 'SMB');
    mesh_updated.Hi     = ncread(mesh_path_update, 'Hi');
    mesh_updated.T2m    = ncread(mesh_path_update, 'T2m');
    mesh_updated.Precip = ncread(mesh_path_update, 'Precip');
    mesh_updated.uabs   = ncread(mesh_path_update,'uabs_surf');
    mesh_updated.BMB    = ncread(mesh_path_update, 'BMB');
    [Hi_fix_mesh_updated, maskHi0_mesh_updated]= Hi0_to_NaN_mesh(mesh_updated.Hi);
end
% calculate Hi differences
mesh_first.Hi_diff=(mesh_first.Hi(:,end))-(mesh_first.Hi(:,1));
% polygon for catchments
ice_boundaries=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/IceBoundaries_Antarctica_v02.shp');

% % calculate mask out of the ice using TriGC
% clean_idx1=~isnan(IM2(:,1));
% IM2_nonan=IM2(clean_idx1,:);
% clean_idx2=~isnan(IM2_nonan(:,2));
% IM2_nonan=IM2_nonan(clean_idx2,:);
% Tri_inside_IM2=inpolygon(mesh_first.TriGC(:,1),mesh_first.TriGC(:,2),IM2(:,1),IM2(:,2));

% plots
% ice thickness
if allow_plot_mesh
    plot_mesh_data_a_RLIS(mesh_first,mesh_first.Hi(:,1).*maskHi0_mesh(:,1));
    clim([0 4000]);
    hold on
    %patch('vertices',V,'faces',Tri,'facecolor','none','edgecolor','r','linewidth',3,'marker','o');
    %plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    %plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','red');
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
    t=title([output_folder,' t0'],'Interpreter','none');
    t.Units='normalized';
    %plot(ice_boundaries(5).X,ice_boundaries(5).Y,'LineWidth',2)
    if allow_save_plots
        print([path_save,output_folder,'_mesh_1_Hi_t0'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh_data_a_RLIS(mesh_updated,mesh_updated.Hi(:,length(time_slice2)).*maskHi0_mesh_updated(:,length(time_slice2)));
        clim([0 4000]);
        hold on
        plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');
        t=title(output_folder);
        t.Units='normalized';
    else
        plot_mesh_data_a_RLIS(mesh_first,mesh_first.Hi(:,end).*maskHi0_mesh(:,end));
        clim([0 4000]);
        hold on
        %plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
        %plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','red');
        plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
        t=title([output_folder,' tf'],'Interpreter','none');
        t.Units='normalized';
        if allow_save_plots
            print([path_save,output_folder,'_mesh_1_Hi_tf'],'-dpng','-r300')
        end
    end
end
% ====================
%  plot for velocities
% ====================
if allow_plot_mesh
    % plot_mesh_data_b_RLIS(mesh_first,log(mesh_first.uabs(:,1)));
    % hold on
    % plot(CF(:,1),CF(:,2),'LineWidth',2,'Color','red');
    % plot(CL(:,1),CL(:,2),'LineWidth',2,'Color','blue');
    % plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green');
    % set(gca,'ColorScale','log')
    % clim([10^-1 10^1]);
    if allow_save_plots
        %print([path_save,output_folder,'_mesh_1_uabs_t0'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh_data(mesh_updated,log(mesh_updated.uabs(:,length(time_slice2))));
        hold on
        plot(IM_mesh2(:,1),IM_mesh2(:,2),'LineWidth',2,'Color','black');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');
        set(gca,'ColorScale','log')
        clim([10^-1 10^1]);
    else
        %plot_mesh_data_b_RLIS(mesh_first,log(mesh_first.uabs(:,end)));
        plot_mesh_data_b_RLIS(mesh_first,log10(mesh_first.uabs(:,end)));
        hold on
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','black');
        plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
        cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/devon/devon.cpt'...
            ,'flip',false,'ncol',100);
        clim([log10(1) log10(2000)])
        cb = colorbar;
        cb.Ticks = log10([1 10 50 100 500 1000 2000]); % positions in log10 space
        cb.TickLabels = {'1','10','50','100','500','1000','2000'}; % readable labels
        cb.Label.String = 'Velocity (m/yr)';
        t=title([output_folder,' tf'],'Interpreter','none');
        t.Units='normalized';
        %set(gca,'ColorScale','log')
        %clim([10^0 10^3]);
        if allow_save_plots
            print([path_save,output_folder,'_mesh_1_uabs_tf'],'-dpng','-r300')
        end
    end
end
%========================
% plot for BMB
% =======================
if allow_plot_mesh
    plot_mesh_data_a_RLIS(mesh_first,mesh_first.BMB(:,1));
    hold on
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','black');
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-2 2]);
    t=title([output_folder,' t0'],'Interpreter','none');
    t.Units='normalized';
    if allow_save_plots
        print([path_save,output_folder,'_mesh_1_BMB_t0'],'-dpng','-r300')
    end
    if allow_mesh_update
        plot_mesh_data_a_RLIS(mesh_updated,mesh_updated.BMB(:,length(time_slice2)));
        hold on
        plot(IM_mesh2(:,1),IM_mesh2(:,2),'LineWidth',2,'Color','black');
        plot(GL_mesh2(:,1),GL_mesh2(:,2),'LineWidth',2,'Color','red');
        cptcmap('GMT_polar','flip',true,'ncol',100);
        clim([-2 2]);
        if allow_save_plots
          print([path_save,output_folder,'_mesh_2_BMB_tf'],'-dpng','-r300')
        end
    else
        plot_mesh_data_a_RLIS(mesh_first,mesh_first.BMB(:,end));
        hold on
        plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','black');
        plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
        plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','green','linestyle', ...
            '-.');
        %plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','red','LineStyle','--');
        cptcmap('GMT_polar','flip',true,'ncol',100);
        clim([-2 2]);
        t=title([output_folder,' tf'],'Interpreter','none');
        t.Units='normalized';
        if allow_save_plots
            print([path_save,output_folder,'_mesh_1_BMB_tf'],'-dpng','-r300')
        end
    end
end
% ========== 
% == plot for Hi_tf - Hi_t0
% ==========
if allow_plot_mesh
    plot_mesh_data_a_RLIS(mesh_first,mesh_first.Hi_diff(:,1));
    hold on
    %cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/bukavu/bukavu.cpt'...
    %        ,'flip',false,'ncol',100);
    colorbar;
    cptcmap('GMT_polar','flip',true,'ncol',100);
    clim([-1000 1000]);
    plot(GL2(:,1),GL2(:,2),'LineWidth',2,'Color','black');
    plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
    t=title([output_folder,' Hi(tf) - Hi(t0)'],'Interpreter','none');
    t.Units='normalized';
    if allow_save_plots
        print([path_save,output_folder,'_mesh_1_Hi_diff'],'-dpng','-r300')
    end
end
% ========== 
% == plot for mask
% ==========
if allow_plot_mesh
    plot_mesh_data_a_RLIS(mesh_first,mesh_first.mask(:,1));
    hold on
    %cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/bukavu/bukavu.cpt'...
    %        ,'flip',false,'ncol',100);
    colorbar;
    plot(GL(:,1),GL(:,2),'LineWidth',2,'Color','black');
    plot(IM(:,1),IM(:,2),'LineWidth',2,'Color','black');
    t=title([output_folder,' Mask'],'Interpreter','none');
    t.Units='normalized';
end

%% different section of the code to plot differences between modelled and observed
% I will use the data from MEaSUREs and grid output from the ROI

uabs_ufe=ncread([ufe_folder_path, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'uabs_surf');
x_ufe=ncread([ufe_folder_path, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'x');
y_ufe=ncread([ufe_folder_path, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'y');
Hi_ufe=ncread([ufe_folder_path, '/main_output_ANT_grid_ROI_RiiserLarsen.nc'],'Hi');
% create a mask where there is ice, to mask the velocity later
[~, maskHi0_ufe]= Hi0_to_NaN(Hi_ufe);
[xx_measures,yy_measures]=meshgrid(x_MEaSUREs,y_MEaSUREs);
[xx_ufe,yy_ufe]=meshgrid(x_ufe,y_ufe);
% interpolate MEaSUREs to same resolution as UFEMISM output
uabs_measures_interp=interp2(xx_measures,yy_measures,uabs_MEaSUREs',xx_ufe,yy_ufe);
uabs_diff=uabs_measures_interp-(uabs_ufe(:,:,end).*maskHi0_ufe(:,:,end))';

% Set up GUI
wa = 750;
ha = 750;

margins_hor = [125,125];
margins_ver = [125,50];

wf = margins_hor( 1) + wa + margins_hor( 2);
hf = margins_ver( 1) + ha + margins_ver( 2);

H.Fig = figure( 'position',[100,100,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[margins_hor(1),margins_ver(1),wa,ha],...
  'xlim',[min(x_ufe)+5e4 max(x_ufe)-5e4],'ylim',[min(y_ufe)+5e4 max(y_ufe)-5e4],'fontsize',24,'xgrid','on','ygrid','on');
hold on
contourf(x_ufe,y_ufe,uabs_diff,100,'LineColor','none');
plot(basins_MEaSUREs(4).X,basins_MEaSUREs(4).Y,'LineWidth',2,'Color','black'); %R-LIS
plot(basins_MEaSUREs(3).X,basins_MEaSUREs(3).Y,'LineWidth',2,'Color','black'); % Brunt
%plot(coast_MEaSUREs.X,coast_MEaSUREs.Y,'LineWidth',2,'Color','black'); % Brunt
%plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');
cptcmap('/Users/frre9931/Documents/PhD/ScientificColourMaps8/vik/vik.cpt'...
            ,'flip',false,'ncol',256);
clim([-800 800])
cb = colorbar;
cb.Label.String = 'Surface ice velocity diffence (m/yr)';
t=title(output_folder,'Interpreter','none');
t.Units='normalized';

%% code to plot the difference between two simulations
output_folder2= 'results_ant_PD_maxphi_20_retreat_mask_code_SMB_and_phi_50percent'; 
ufe_folder_path2=['/Users/frre9931/Desktop/tetralith_results/', output_folder2];

mesh_path_sim2= [ufe_folder_path2, '/main_output_ANT_00001.nc']; %initial state

% read_mesh_from_file
mesh_sim2=read_mesh_from_file(mesh_path_sim2);
mesh_sim2.Hb            = ncread( mesh_path_sim2,'Hb');

% now both simulations have a value for Hb
% first is stored in mesh_first.Hb and latter mesh_sim2.Hb
% as we are working with same meshs it is fine to just substract

% first - sim2 for the (end) time-step
Hb_difference = mesh_first.Hb(:,end) - mesh_sim2.Hb(:,end);
plot_mesh_data_a_RLIS(mesh_first, Hb_difference);
hold on
plot(IM2(:,1),IM2(:,2),'LineWidth',2,'Color','black');



