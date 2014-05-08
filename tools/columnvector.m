function [out] = columnvector(in)
% columnvector transposes the input vector to a column vector if it is not
% already a column vector.
%
% joerg buescher
%

dim = size(in) ;
if dim(1) < dim(2)
    out = in' ;
else
    out = in ;
end
