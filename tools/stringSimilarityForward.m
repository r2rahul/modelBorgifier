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
length1 = length(input1) ;
length2 = length(input2) ;

% reduce word size for short inputs
wordsize = max([ min([wordsize, length1-1, length2-1]) , 2]) ;

% % quick and dirty implementations, finds 'b' == ' ' and 'y' == 'd' etc...
% v1 = sqrt(double(input1)' * double(input2)) ;
% v2 = (v1 - floor(v1)) < 1e-9 ;

% % nicer implementation, only slightly slower but more correct
v2 = (double(input1)' * ones(1,length2)) == (ones(length1,1) * double(input2)) ;

v3 = v2(wordsize:end,wordsize:end) .* v2(1:(end-wordsize+1),1:(end-wordsize+1)) ;
%score = sum(sum(v3)) / ((length1-wordsize) * (length2-wordsize)) ;
% better weighting
score = min([1, sum(sum(v3)) / (max([length1 length2]) + abs(length1 -length2) ) ]) ;

score(isnan(score)) = 0 ;

% % for testing only
% for i1 = 1:length(input1)
%     for i2 = 1:length(input2)
%         t{i1,i2} = [input1(i1) input2(i2)] ;
%     end
% end