% Code to eliminate the ice rises in R-LIS from the bathymetry of Bedmachine
% add functions from UFEMISM library
path(path,genpath('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/tools/matlab'));

% load Bedmachine
Hb=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc','Hb');
Hi=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc','Hi');
Hs=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc','Hs');

x=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc','x');
y=ncread('/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc','y');
Hb_smoothed=smoothdata2(Hb',"movmean",50); %I do not remember if was with 50 or 70 the one that worked...
Hs_smoothed=smoothdata2(Hs',"movmean",10);
Hi_smoothed=smoothdata2(Hi',"movmean",100);
Hb_smoothed200=smoothdata2(Hb',"movmean",200); %estaba con 50
[x_grid,y_grid]=meshgrid(x,y);

% load polygon from R-LIS as reference in plots
RLIS=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ROI_for_R-LIS.shp');
% load GL from measures
GL=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/GroundingLine_Antarctica_v02.shp');
% load polygon created from the intersection of grounded ice and ice rises
% from the ROI
ice_rises=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ice_rises_intersected.shp');
% load the same file buth with a buffer of 1 km to select more values for
% smoothing
ice_rises_buffer=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/ice_rises_buffer5km.shp');
% another file with a more extreme case smoothing almost everything
ice_rises_half=shaperead('/Users/frre9931/Documents/PhD/RiiserLarsen/section_and_IR.shp');

figure()
contourf(x,y,Hb')
hold on
plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red'); % plot ROI
plot(GL.X,GL.Y,'LineWidth',2,'Color','red'); % plot grounding line
plot(ice_rises.X,ice_rises.Y,'LineWidth',2,'Color','cyan');

% check if is inside polygon ice_rises
%[in,on]=inpolygon(x_grid,y_grid,ice_rises_buffer.X,ice_rises_buffer.Y);
[in,on]=inpolygon(x_grid,y_grid,ice_rises_half.X,ice_rises_half.Y);

%[in_half,on_half]=inpolygon(x_grid,y_grid,ice_rises_half.X,ice_rises_buffer.half);
% assign value of Hb' to Hb_modified
Hb_modified=Hb';
Hb_modified200=Hb';
Hb_modified_section=Hb';
% do the same for Hs and Hi
Hs_modified=Hs';
Hi_modified=Hi';
% assign a smoothed value to the points inside the polygon
Hb_modified(in)=Hb_smoothed(in);
Hb_modified200(in)=Hb_smoothed200(in);
Hs_modified(in)=Hs_smoothed(in);
Hi_modified(in)=Hi_smoothed(in);

% plot modified Hb
figure()
contourf(x,y,Hb_modified)
hold on
plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
plot(GL.X,GL.Y,'LineWidth',2,'Color','red');

% plot of grid points inside polygon to double check that is not empty
figure()
contourf(x_grid,y_grid,in)
hold on
plot(GL.X,GL.Y,'LineWidth',2,'Color','red');

% load the GL from UFEMISM to check if it modify all the needed values
% as the resolution of Bedmachine is 5 km not all the ice rises are
% selected during the inpolygon

mesh_first=read_mesh_from_file("/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20/main_output_ANT_00001.nc");
GL_ufe=ncread("/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20/main_output_ANT_00001.nc",'grounding_line',[1,1,1], [size(mesh_first.E,1),2,1]);
plot(GL_ufe(:,1),GL_ufe(:,2),'LineWidth',2,'Color','cyan')

figure()
contourf(x,y,Hs-Hs_modified)
%contourf(x,y,Hb'-Hb_new)
hold on
plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
plot(GL.X,GL.Y,'LineWidth',2,'Color','red');

Hb_new=Hb'-abs(Hb'-Hb_modified)-abs(Hb'-Hb_modified200);
figure()
contourf(x,y,Hb'-Hb_new)
hold on
plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
plot(GL.X,GL.Y,'LineWidth',2,'Color','red');
% % repeat plot for Hi and Hs and see
% figure()
% contourf(x,y,Hs')
% hold on
% plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
% plot(GL.X,GL.Y,'LineWidth',2,'Color','red');
% 
% figure()
% contourf(x,y,Hi')
% hold on
% plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
% plot(GL.X,GL.Y,'LineWidth',2,'Color','red');

% is not working! first think if there is a way to find if the ice is
% grounded or not and modify it from there. Or other option is to add a
% buffer in the shapefile? so I can smooth more points? even tho both
% should work together
%Hb_new=Hb_new';
mask_grounded=zeros(size(Hb_new));
Hs=Hs';
Hi=Hi';
for i=1:size(Hb,1)
    for j=1:size(Hb,2)
        if Hs(i,j)-Hs_modified(i,j)<0
            Hs_modified(i,j)=Hs(i,j);
        end
        if Hi(i,j)-Hi_modified(i,j)<0
        %if Hi(i,j)-(Hs(i,j)-Hs_modified(i,j))<0
            Hi_modified(i,j)=Hi(i,j);
        end
        if Hb_new(i,j)>=Hs(i,j)-Hi(i,j)
            mask_grounded(i,j)=1;
        end
    end
end
figure()
contourf(x,y,mask_grounded)

figure()
contourf(x,y,Hs_modified)%abs(Hs'-Hs_modified))
hold on
plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
plot(GL.X,GL.Y,'LineWidth',2,'Color','red');

figure()
contourf(x,y,Hi-Hi_modified)
hold on
plot(RLIS.X,RLIS.Y,'LineWidth',2,'Color','red');
plot(GL.X,GL.Y,'LineWidth',2,'Color','red');

Hs_new=Hs-abs(Hs-Hs_modified);
Hi_new=Hi-abs(Hi-Hi_modified); %apply same rate of change to thickness to be sure is consistent
%% last idea to be sure this is working or not.. 
% interpolate the data to a 10km dataset as UFE
x_10km=x(1:2:end);
y_10km=y(1:2:end);
[xx_10km,yy_10km]=meshgrid(x_10km,y_10km);
Hb10km=interp2(x_grid,y_grid,Hb_new,xx_10km,yy_10km);
Hi10km=interp2(x_grid,y_grid,Hi_modified,xx_10km,yy_10km);
Hs10km=interp2(x_grid,y_grid,Hs_modified,xx_10km,yy_10km);
mask_grounded2=zeros(size(Hb10km));
for i=1:size(Hb10km,1)
    for j=1:size(Hb10km,2)
        if Hb10km(i,j)>=Hs10km(i,j)-Hi10km(i,j)
            mask_grounded2(i,j)=1;
        end
    end
end
figure()
contourf(x_10km,y_10km,mask_grounded2)
%% save the new file with all the same variables as Bedmachine
% need to load Hi, Hs and idk what else..
% save it as Bedmachine_smoothed_Hb_ice_rises.nc

% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/RiiserLarsen/Bedmachine_v3_smoothed_Hb_ice_rises.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(x,1));
dim_y = netcdf.defDim(ncid,'y',size(y,1));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);

id_Hs  = netcdf.defVar(ncid,'Hs','double',[dim_x, dim_y]);
id_Hb  = netcdf.defVar(ncid,'Hb','double',[dim_x, dim_y]);
id_Hi  = netcdf.defVar(ncid,'Hi','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 

netcdf.putAtt(ncid,id_Hs,'standard_name','Surface topography'); 
netcdf.putAtt(ncid,id_Hs,'units','meters');

netcdf.putAtt(ncid,id_Hb,'standard_name','Bedrock elevation'); 
netcdf.putAtt(ncid,id_Hb,'units','meters');

netcdf.putAtt(ncid,id_Hi,'standard_name','Ice thickness'); 
netcdf.putAtt(ncid,id_Hi,'units','meters');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,x);
netcdf.putVar(ncid,id_y,y);

netcdf.putVar(ncid,id_Hs,Hs');
netcdf.putVar(ncid,id_Hb,Hb_new');
netcdf.putVar(ncid,id_Hi,Hi');

% Close file
% ==========

netcdf.close(ncid);