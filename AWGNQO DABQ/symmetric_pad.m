function [padded] = symmetric_pad(x,m)
%SYMMETRIC_PAD Pad x to length m with zeros, such that result is symmetric 
%about the middle
%   If length(x) and m don't have the same parity, more zeros will be added
%   to the end. Works for column vectors only.

diff = m-length(x);
padded = [zeros(floor(diff/2), 1); x; zeros(ceil(diff/2), 1)];
end

