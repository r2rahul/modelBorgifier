function out = ismatlab
out = false ;
if ~isempty(strfind(version,'R20'))
    out = true ;
end