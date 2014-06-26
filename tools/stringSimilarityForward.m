% stringSimilarityForward compares two strings and outputs a simple
% similarity score
%
%  score = stringSimilarityForward(input1, input2, wordsize)
%
%INPUT
%  input1:   first string to compare
%  input2:	 second string to compare
%  wordsize: integer that defines the number of letters that are compared at a time (3 is a good value)   
%
%CALLED BY
%  findMetMatch
%  compareCbModels
%
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
function score = stringSimilarityForward(input1, input2, wordsize)
% compares the similarity of two strings and returns a score between 0 (not
% similar and 1 (identical)

% shortcut if input 1 and input 2 are identical
if strcmp(input1, input2)
    score = 1 ;
    return
end

score = zeros(length(input1)-wordsize,1) ;
for i1 = 1:length(input1)-wordsize
    for i2 = 1:length(input2)-wordsize
        score(i1) = score(i1) + strcmpi(input1(i1:i1+wordsize), input2(i2:i2+wordsize)) ;
    end
end
score = mean(score ./ max(score)) ;

score(isnan(score)) = 0 ;