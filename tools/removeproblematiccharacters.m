% this file is published under Creative Commons BY-NC-SA
% 
% Assimilating genome-scale metabolic reconstructions with modelBorgifier
% in preparation
%
% John T. Sauls and Joerg M. Buescher
% BRAIN Aktiengesellschaft
% Microbial Production Technologies Unit
% Quantitative Biology and Sequencing Platform
% Darmstaeter Str. 34-36
% 64673 Zwingenberg, Germany
% www.brain-biotech.de
% jrb@brain-biotech.de
% 
%
function out = removeproblematiccharacters(in,varargin)
% replaces all problematic characters ins input string in with substitute
% sub
if isempty(in)
    out = in ;
    return
end
if isempty(varargin)
    sub = '_' ;
else
    sub = varargin{1} ;
end

probpos =  [strfind(in,' ') ...
            strfind(in,'@') ...
            strfind(in,'"') ...
            strfind(in,'�') ...
            strfind(in,'%') ...
            strfind(in,'&') ...
            strfind(in,'/') ...
            strfind(in,'(') ...
            strfind(in,')') ...
            strfind(in,'=') ...
            strfind(in,'?') ...
            strfind(in,'�') ...
            strfind(in,'!') ...
            strfind(in,'#') ...
            strfind(in,'*') ...
            strfind(in,'+') ...
            strfind(in,'-') ...
            strfind(in,'~') ...
            strfind(in,'"') ...
            strfind(in,'<') ...
            strfind(in,'>') ...
            strfind(in,',') ...
            strfind(in,';') ...
            strfind(in,'.') ...
            strfind(in,':') ...
            strfind(in,'\') ] ;
out = in ;
if isempty(sub)
    out(probpos) = [] ;
else
    if length(sub) == 1
        out(probpos) = sub ;
    else
        error('removeproblematiccharacters: substitute string too long')
    end
end

probpos = [strfind(out,'�') strfind(out,'�')] ;
if ~isempty(probpos) ; out = [out(1:probpos-1) 'oe' out(probpos+1:end)] ; end
probpos = [strfind(out,'�') strfind(out,'�')] ;
if ~isempty(probpos) ; out = [out(1:probpos-1) 'ue' out(probpos+1:end)] ; end
probpos = [strfind(out,'�') strfind(out,'�')] ;
if ~isempty(probpos) ; out = [out(1:probpos-1) 'ae' out(probpos+1:end)] ; end
probpos = strfind(out,'�') ;
if ~isempty(probpos) ; out = [out(1:probpos-1) 'ss' out(probpos+1:end)] ; end
probpos = strfind(out,'�') ;
if ~isempty(probpos) ; out = [out(1:probpos-1) 'mu' out(probpos+1:end)] ; end


