function [ model ] = svmtrain( varargin )
%SVMTRAIN 
% model = svmtrain(training_label_vector, training_instance_matrix [, 'libsvm_options']);
%
%        -training_label_vector:
%            An m by 1 vector of training labels (type must be double).
%        -training_instance_matrix:
%            An m by n matrix of m training instances with n features.
%            It can be dense or sparse (type must be double).
%        -libsvm_options:
%            A string of training options in the same format as that of LIBSVM.
%           Commonly used options are:
%                -v n    n-fold cross-validation
%                -t 0    linear kernel
%                -t 2    radial basis (default)
%                -s 0    SVC type: 0 = C-SVC 3 = epsilon-SVM
%                -C x    C-parameter value = x (default = 1)
%                -g x    gamma parameter value = x
%                -p x    epsilon parameter in loss function = x (try 0.1)
%

model = svmtrain(varargin) ;

end

