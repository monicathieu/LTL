function [out] = asRow(vec)
%UNTITLED4 makes sure a vector is a column vector

sz = size(vec);
nonSingletonDims = sz>1;

if sum(nonSingletonDims>1)
    error('data must only have one non-singleton dimension')
end

vecLength = sz(nonSingletonDims);

out = reshape(vec,1,vecLength);

end
