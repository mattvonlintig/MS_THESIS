function [outMat] = shiftMat(inMat,ds,dim)
%SHIFTMAT Shifts indvidual rows/columns of a data matrix
%
%   INPUTS
%       inMat   -   data matrix [m by n]
%       ds      -   vector of offset values... SIZE(ds) == m,n, or 1
%       dim     -   dimension to shift on
%
%   OUTPUTS
%       outMat  -   shifted data matrix
%

if dim == 2
    inMat = inMat';
end

for ii = 1:numel(ds)
    inMat(:,ii) = circshift(inMat(:,ii),ds(ii));
end

if dim == 2
    outMat = inMat';
else
    outMat = inMat;
end

end