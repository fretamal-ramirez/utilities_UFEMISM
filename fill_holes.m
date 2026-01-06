% code to fill holes in the dataset of Bedmachine

% the idea is to fill the holes where Hi==0 and Hb<0, using the second min
% value of the neighbours.
path_bedmachine='/Users/frre9931/Desktop/UFEMISM2.0_main/UFEMISM2.0/data/Bedmachine_Antarctica/Bedmachine_v3_Antarctica_5km.nc';
hi=ncread(path_bedmachine,'Hi');
hb=ncread(path_bedmachine,'Hb');
hs=ncread(path_bedmachine,'Hs');
xi=ncread(path_bedmachine,'x');
yi=ncread(path_bedmachine,'y');
[x_grid,y_grid]=meshgrid(xi,yi);

hi_filled = hi; % copy to modify
hs_filled = hs;

[nx, ny] = size(hi);

mask_filledvalues=zeros(size(hi));
for i = 1:nx
    for j = 1:ny
        if hi(i,j) ==0 && hb(i,j) < 0

            % Define local window boundaries (safe for edges)
            i_min = max(i-1, 1);
            i_max = min(i+1, nx);
            j_min = max(j-1, 1);
            j_max = min(j+1, ny);

            % Extract neighborhood values and coordinates
            neigh_hi = hi(i_min:i_max, j_min:j_max);
            neigh_hs = hs(i_min:i_max, j_min:j_max); % surface elevation array
            [rows_local, cols_local] = find(neigh_hi > 0);

            if numel(rows_local) >= 6 % at least 6 valid ice neighbors

                % Get neighbor thickness values
                vals = neigh_hi(sub2ind(size(neigh_hi), rows_local, cols_local));

                % Sort by value and get the second minimum (or first)
                [sorted_vals, sort_idx] = sort(vals);
                fill_val = sorted_vals(1); % use first (min) or 2nd (sorted_vals(2))
                chosen_idx = sort_idx(1);

                % Map chosen local index back to global coordinates
                local_row = rows_local(chosen_idx);
                local_col = cols_local(chosen_idx);
                global_row = i_min + local_row - 1;
                global_col = j_min + local_col - 1;

                % Assign values
                hi_filled(i,j) = fill_val;
                hs_filled(i,j) = hs(global_row, global_col); % elevation of chosen neighbor
                mask_filledvalues(i,j) = 1;
            end
        end
    end
end

% repeat the same if Hs>0 and Hs<0.11
for i = 1:nx
    for j = 1:ny
        if hs(i,j)<1.0 && hb(i,j) < 0

            % Define local window boundaries (safe for edges)
            i_min = max(i-1, 1);
            i_max = min(i+1, nx);
            j_min = max(j-1, 1);
            j_max = min(j+1, ny);

            % Extract neighborhood values and coordinates
            neigh_hi = hi(i_min:i_max, j_min:j_max);
            neigh_hs = hs(i_min:i_max, j_min:j_max); % surface elevation array
            [rows_local, cols_local] = find(neigh_hs > 1.0);

            if numel(rows_local) >= 6 % at least 6 valid ice neighbors

                % Get neighbor thickness values
                vals = neigh_hs(sub2ind(size(neigh_hs), rows_local, cols_local));

                % Sort by value and get the second minimum (or first)
                [sorted_vals, sort_idx] = sort(vals);
                if isempty(sorted_vals)
                    %do nothing
                    %hs_filled(i,j) = 0.1; % keep it as it is
                else
                    fill_val = sorted_vals(1); % use first (min) or 2nd (sorted_vals(2))
                    chosen_idx = sort_idx(1);

                    % Map chosen local index back to global coordinates
                    local_row = rows_local(chosen_idx);
                    local_col = cols_local(chosen_idx);
                    global_row = i_min + local_row - 1;
                    global_col = j_min + local_col - 1;

                    % Assign values
                    hs_filled(i,j) = fill_val;
                    hi_filled(i,j) = hi_filled(global_row, global_col); % elevation of chosen neighbor
                    mask_filledvalues(i,j) = 2;
                 end
            end
        end
    end
end


% 
% 
% 
% for i = 1:nx
%     for j = 1:ny
%         if hs(i,j) > 0 && hs(i,j) < 0.11
% 
%             % Define local window boundaries (safe for edges)
%             i_min = max(i-1, 1);
%             i_max = min(i+1, nx);
%             j_min = max(j-1, 1);
%             j_max = min(j+1, ny);
% 
%             % Extract neighborhood
%             neigh = hs(i_min:i_max, j_min:j_max);
%             neigh = neigh(:);
%             neigh_nozero = neigh(neigh > 0);
% 
%             % Count nonzero neighbors
%             num_nonzero = numel(neigh_nozero);
% 
%             % Only fill if surrounded by ice (e.g. 7 or 8 ice neighbors)
%             if num_nonzero >= 6
%                 neigh_sorted = sort(neigh_nozero(neigh_nozero > 0.11));
%                 if isempty(neigh_sorted)
%                     hs_filled(i,j) = 0.1; % keep it as it is
%                 else
%                     hs_filled(i,j) = neigh_sorted(1); % second minimum
%                     mask_filledvalues(i,j)=2; % store in the mask
%                 end
%             end
%         end
%     end
% end
% 
% mask_Hi=zeros(size(hi));
% mask_Hs=zeros(size(hs));
% hs_mod=hs;
% for i=1:nx
%     for j=1:ny
%         if mask_filledvalues(i,j)==1
%             hs_mod(i,j)=hi_filled(i,j);
%         elseif mask_filledvalues(i,j)==2
%             hs_mod(i,j)=hs_filled(i,j);
%         end
%         if hi_filled(i,j)>0
%             mask_Hi(i,j)=1;
%         end
%         if hs_mod(i,j)>=0.11 % 0.1 is like non value...
%             mask_Hs(i,j)=1;
%         end
%     end
% end

%% plots to check if everything is allright
coast_MEaSUREs=shaperead('/Users/frre9931/Documents/PhD/MEaSUREs/Coastline_Antarctica_v02.shp');
GL_UFE_HR=ncread('/Users/frre9931/Desktop/tetralith_results/results_ant_PD_inversion_dHdt_init_R-LIS_gamma20_PMP_roughness_max20_HR/main_output_ANT_00001.nc','grounding_line',[1,1,1],[109694,2,1]);
figure()
contourf(x_grid,y_grid,(mask_filledvalues)');
hold on
plot(coast_MEaSUREs.X,coast_MEaSUREs.Y,'red','LineWidth',2)
plot(GL_UFE_HR(:,1),GL_UFE_HR(:,2),'black','LineWidth',2)

figure()
contourf(x_grid,y_grid,(hi_filled)',50,'linecolor','none')
hold on
plot(GL_UFE_HR(:,1),GL_UFE_HR(:,2),'black','LineWidth',2)

figure()
contourf(x_grid,y_grid,(hs_filled)',50,'linecolor','none')
hold on
plot(GL_UFE_HR(:,1),GL_UFE_HR(:,2),'black','LineWidth',2)

figure()
contourf(x_grid,y_grid,(hb)',50,'linecolor','none')
hold on
plot(GL_UFE_HR(:,1),GL_UFE_HR(:,2),'black','LineWidth',2)

%% save the dataset hi_filled

% create netcdf file
ncid = netcdf.create('/Users/frre9931/Documents/PhD/BedMachine/Bedmachine_5km_filled.nc','CLOBBER');

% Define dimensions
% =================

dim_x = netcdf.defDim(ncid,'x',size(xi,1));
dim_y = netcdf.defDim(ncid,'y',size(yi,1));

% Define variables
% ================

id_x = netcdf.defVar(ncid,'x','double',dim_x);
id_y = netcdf.defVar(ncid,'y','double',dim_y);
id_Hi  = netcdf.defVar(ncid,'Hi','double',[dim_x, dim_y]);
id_Hb  = netcdf.defVar(ncid,'Hb','double',[dim_x, dim_y]);
id_Hs  = netcdf.defVar(ncid,'Hs','double',[dim_x, dim_y]);

% Add information to variables
% ==================

netcdf.putAtt(ncid,id_x,'standard_name','x-axis distance from center of projection'); 
netcdf.putAtt(ncid,id_x,'units','m'); 
netcdf.putAtt(ncid,id_y,'standard_name','y-axis distance from center of projection');
netcdf.putAtt(ncid,id_y,'units','m'); 

netcdf.putAtt(ncid,id_Hi,'standard_name','Ice thickness'); 
netcdf.putAtt(ncid,id_Hi,'units','m');
netcdf.putAtt(ncid,id_Hb,'standard_name','Bedrock elevation'); 
netcdf.putAtt(ncid,id_Hb,'units','m');
netcdf.putAtt(ncid,id_Hs,'standard_name','Surface elevation'); 
netcdf.putAtt(ncid,id_Hs,'units','m');

% End definition mode
% ===================

netcdf.endDef(ncid)

% Save data
% =========

netcdf.putVar(ncid,id_x,xi);
netcdf.putVar(ncid,id_y,yi);

netcdf.putVar(ncid,id_Hi,hi_filled);
netcdf.putVar(ncid,id_Hb,hb);
netcdf.putVar(ncid,id_Hs,hs_filled);


% Close file
% ==========

netcdf.close(ncid);
