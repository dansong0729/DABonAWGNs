function [pYX] = getawgnqtransition(xsupport, q, N)
%GETAWGNQTRANSITION Get transition matrix of AWGN-QO channel.
%   INPUTS:
%   xsupport: input pmf support points, as column vector
%   q: bin thresholds in ascending order, as column vector
%   N: channel noise power

%   OUTPUTS: dim(q)-1 by size(xsupport) matrix, each column being a
%   probability vector
%   
%   If q is incomplete (range is not (-Inf, Inf)), only part of channel
%   transition matrix is output

%col-row makes matrix of pairwise differences
offsets = q-xsupport.';
left = offsets(1:end-1, :); %@ left boundary of each bin
right = offsets(2:end, :); %@ right boundary of each bin
pYX = left; %get output matrix of right size
%subtract cdf of left from right for bin probability
%better precision when using normcdf at negative values
negate = (left>=0)&(right>=0); %should flip when both positive offsets
pYX(~negate) = normcdf(right(~negate), 0, N) - normcdf(left(~negate), 0, N);
pYX(negate) = normcdf(-left(negate), 0, N) - normcdf(-right(negate), 0, N);
end