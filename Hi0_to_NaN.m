function [field3d, mask]= Hi0_to_NaN(Hi)

% Hi is the ice thickness in a 3D matrix x,y,val
% outputs
% field3d is Hi with 0 as NaNs
% mask, ice=1, no ice=NaN, to later on multiply it to any field for that
% time

% fill in with NaNs the Hi=0
mask=ones(size(Hi));
field3d=Hi;

for i=1:size(field3d,1)
    for j=1:size(field3d,2)
        for m=1:size(field3d,3)
            if field3d(i,j,m)==0
                field3d(i,j,m)=NaN;
                mask(i,j,m)=NaN;
            end
        end
    end
end

end

