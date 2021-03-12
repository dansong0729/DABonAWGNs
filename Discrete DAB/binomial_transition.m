function [pYX] = binomial_transition(xsupp,n)
%BINOMIAL_TRANSITION Transition matrix of binomial channel
%   Just a thin wrapper around binopdf, for convenience
%
%   INPUTS:
%   xsupp: Input (X) pmf support points, aka values of p, as col. Must be
%   in [0,1].
%   n: Scalar parameter of binomial channel; number of trials
%
%   OUTPUTS: n+1 by length(xsupp) matrix, each column being a
%   probability vector

pYX = zeros(n+1, length(xsupp));
for i=0:n
    pYX(i+1, :) = binopdf(i, n, xsupp);
end

