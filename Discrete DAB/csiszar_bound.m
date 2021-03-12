function [bound, xstar] = csiszar_bound(pX, xsupp, pYX_func, dpYX_func, y_range, A)
%CSISZAR_BOUND Finds Csiszar min-max upper bound and x that achieves it.
%   Detailed explanation goes here
%   INPUTS
%   pX: Input pmf, as column probability vector
%   xsupp: Input pmf support points in ASCENDING order, as column vector
%   pYX_func: Function handle taking (y,x) as input and returns p_Y|X
%   dpYX_func: Function handle taking (y,x) as input and returns d/dx p_Y|X
%   y_range: Vector of possible y values
%   A: Amplitude constraint


Q = get_transition(xsupp, y_range, pYX_func);
pY = Q*pX;

%Get local maxima
%Find d/dx sign changes
num_samples = length(pX)*25;
samples = linspace(-A, A, num_samples);
%compute derivative
derivative = discrete_rel_derivative(samples, pY, pYX_func, dpYX_func, y_range);
%index of sample preceding + to - sign change (assumes not exactly 0)
indices = find(sign(derivative(1:end-1))==1&sign(derivative(2:end))==-1);
%find exact zeros of derivative
fun = @(x) discrete_rel_entropy(PY, PYgivenX, x, xsupp, N);
maxes = zeros(size(indices)); %maximum values
xes = zeros(size(indices)); %x values
for i = 1:length(maxes)
    xval = fzero(fun, [samples(indices(i)), samples(indices(i)+1)]);
    maxes(i) = discrete_rel_entropy(pY, Q(:, i));
    xes(i) = xval;
end

%get global max
[bound, idx] = max(maxes);
xstar = xes(idx);
end


