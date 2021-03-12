function [D] = discrete_rel_entropy(p,q)
%DISCRETE_REL_ENTROPY D(p||q) for discrete pmfs
%   INPUTS: Two n-D arrays of same shape, with all elements summing to 1.
%   OUTPUT: Relative entropy in bits

mask = (p==0); %causes NaN issues
summands = p.*reallog(p./q); %reallog warns on invalid inputs
summands(mask) = 0; %0log(0) = 0 by convention
D = sum(summands, 'all');
D = D/log(2); %from nats to bits
end

