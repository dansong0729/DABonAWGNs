function [pYX] = getawgnqtransition(xsupport, q, N)
%GETAWGNQTRANSITION Get transition matrix of AWGN-QO channel.
%   INPUTS:
%   xsupport: column vector of input distribution support points
%   q: column vector of boundaries, ascending
%   N: channel noise power

%   OUTPUTS: dim(q)-1 by size(xsupport) matrix, each column being a
%   probability vector
%   If q is incomplete (range is not (-Inf, Inf)), only part of channel
%   transition matrix is output

%row - col makes a matrix of the differences

pYX = normcdf(q-xsupport.', 0, N);
%subtract cdf at left boundary from cdf at right boundary for bin
%probability
pYX = pYX(2:end, :) - pYX(1:end-1, :);
end

