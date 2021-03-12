function [derivative] = binomial_derivative(n, k, p)
%BINOMIAL_DERIVATIVE d/dp of binomial(n,k,p)
%   Thin wrapper around MATLAB functions for convenience. Used in computing
%   the derivative of p(Y|X) of binomial channel w.r.t X support point.
%   
%   For precision issues, p<1/2 is preferred. Use symmetry to enforce this.
%   
%   WARNING: For now, outputs Inf/NaN when p = 0,1
%   TODO: Testing

derivative = binopdf(k, n, p).*(k-n.*p)./p./(1-p);

end

