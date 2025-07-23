function data_filled = fillNaNwithNearest(data)
    % Replace NaNs in a 3D matrix with the nearest non-NaN value
    
    % Logical mask of NaNs
    nan_mask = isnan(data);
    
    % Early exit if there are no NaNs
    if ~any(nan_mask(:))
        data_filled = data;
        return
    end
    
    % Replace NaNs with 0 temporarily to avoid errors
    temp_data = data;
    temp_data(nan_mask) = 0;
    
    % Compute distance transform and nearest neighbor indices
    [~, idx] = bwdist(~nan_mask);  % indices of nearest non-NaN values
    
    % Fill NaNs with nearest neighbor values
    data_filled = data;
    data_filled(nan_mask) = data(idx(nan_mask));
end