function [ varargout ] = svmpredict( varargin )
%SVMPREDICT Summary of this function goes here
%   [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model [, 'libsvm_options']);
%
%        -testing_label_vector:
%            An m by 1 vector of prediction labels. If labels of test
%            data are unknown, simply use any random values. (type must be double)
%        -testing_instance_matrix:
%            An m by n matrix of m testing instances with n features.
%            It can be dense or sparse. (type must be double)
%        -model:
%            The output of svmtrain.
%        -libsvm_options:
%            A string of testing options in the same format as that of LIBSVM.

varargout = svmpredict(varargin) ;

end

