function [out] = rowvector(in)
% rowvector transposes the input vector to a row vector if it is not
% already a row vector.
%
% joerg buescher
%

dim = size(in) ;
if dim(1) > dim(2)
    out = in' ;
else
    out = in ;
end
