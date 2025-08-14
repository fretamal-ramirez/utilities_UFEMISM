function [field2d, mask]= Hi0_to_NaN_mesh(Hi)

% Hi is the ice thickness in a 2D matrix vi,time
% outputs
% field2d is Hi with 0 as NaNs
% mask, ice=1, no ice=NaN, to later on multiply it to any field for that
% time

% fill in with NaNs the Hi=0
mask=ones(size(Hi));
field2d=Hi;

for i=1:size(field2d,1)
    for j=1:size(field2d,2)
        if field2d(i,j)==0
            field2d(i,j)=NaN;
            mask(i,j)=NaN;
        end
    end
end

end