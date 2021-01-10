function [q_new] = trisect(q)
%TRISECT Add thresholds as in "Progressive Reads" see reference
%   INPUT: Original thresholds column q (with +- Inf) (length>=4)
%   OUTPUT: Threshold column vector with one new threshold on each side of
%   each original threshold, equally space (not optimized yet)

n = length(q);

%left/right boundaries for finite thresholds
left = q(2:end-2)';
right = q(3:end-1)';

%new thresholds
lnew = 2/3*left + right/3;
rnew = left/3 + 2/3*right;

%interleave
stacked = [left; lnew; rnew];

q_new = zeros(3*n-4, 1);
q_new(3:end-3) = stacked(:);
q_new(end-2) = right(end);
q_new(1) = -Inf;
q_new(end) = Inf;

%new "outer" thresholds
q_new(2) = 2*left(1) - stacked(2);
q_new(end-1) = 2*right(end)-stacked(end);
end