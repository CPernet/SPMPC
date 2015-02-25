function [out1, out2]=spmpc_sample(in1,in2)

% do some resampling of the data
% this is bootstrap not permutation
% this means once we take some data
% one puts it back in the sample
% the same data point can be taken
% several times
%
% output index (out1) and new sample (out2)
% from the nb of data points (in1)
% and the data themselves (in2)

out1 = ceil( rand(in1,1)*length(in2) );
out2 = in2(out1);
