% code to plot transects outputs from UFEMISM
% example for VeststraumenA

path_to_folder = '/Users/frre9931/Desktop/arrhenius_results/';

% name of the run
output_folder = 'results_R-LIS_init_BCAP_H1_O1_SMB1_GHF2_2km';

% load netcdf file
transectA.Hi = ncread([path_to_folder,output_folder,'/transect_VeststraumenA.nc'],'Hi');
transectA.Hb = ncread([path_to_folder,output_folder,'/transect_VeststraumenA.nc'],'Hb');
transectA.Hs = ncread([path_to_folder,output_folder,'/transect_VeststraumenA.nc'],'Hs');
transectA.V = ncread([path_to_folder,output_folder,'/transect_VeststraumenA.nc'],'V'); % transect vertex coordinate

% calculate the distance from Vertex 1 to Vertex LAST
for i = 1:length(transectA.V(:,1))
    distA(i) = sqrt( (transectA.V(i,1)-transectA.V(1,1))^2 ...
                      + (transectA.V(i,2)-transectA.V(1,2))^2  );
end

% plot for the last time
figure()
plot(distA/1000,transectA.Hb(:,end))
hold on
plot(distA/1000,transectA.Hs(:,end))
xlabel('Distance of profile (km)')

% continue working on velocities and strain