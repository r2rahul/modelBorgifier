% this file is based on Matlab Support Solution ID 1-D782JW
% 
function outHtml = colText(inText, inColor)
% return a HTML string with colored font

if isnumeric(inColor)
    if numel(inColor) ~= 3
        error('colText: color must be 3-value RGB or html-compatible word')
    end
    inColor = ['#' dec2hex(round(inColor(1)),2) dec2hex(round(inColor(2)),2) dec2hex(round(inColor(3)),2)] ;
end
outHtml = ['<html><font color="', ...
    inColor, ...
    '">', ...
    inText, ...
    '</font></html>'];