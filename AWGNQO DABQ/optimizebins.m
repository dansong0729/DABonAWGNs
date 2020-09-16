function [q] = optimizebins(pX, xsupport, q_init, N, isSymmetric)
%OPTIMIZEBINS Find optimal bins for AWGNQO channel given input pmf
%   Performs alternating optimization on each bin threshold using fminbnd
%   and terminates on sufficiently small mutual information gain between
%   iterations.
%   
%   INPUTS:
%   pX: input pmf, as column probability vector
%   xsupport: input pmf support points in ASCENDING order, as column vector
%   q_init: initial guess for q
%   N: channel noise power
%   isSymmetric: set true to exploit xsupport symmetry about 0
%   
%   OUTPUT:
%   q: vector of optimal bin thresholds in ascending order, as column
%   vector and including -Inf and Inf
%
%TODO: Characterize and fix precision issues.


endgap = 1e-4; %1e-4 %threshold of MI change for convergence
maxIter = 1000;

%initialize q to a guess, TODO: improve this (MLE?)
%use midpoints
q = q_init;

qsize = size(q,1); %used a lot when isSymmetric

options = optimset('FunValCheck','on'); %for safety

for iter = 1:maxIter
    oldMI = discreteMI(pX, getawgnqtransition(xsupport, q, N));
    indices = 2:qsize-1;
    %only do half if symmetric
    if isSymmetric
        indices = 2:floor(qsize/2);
        if mod(qsize,2)
            %middle element is 0 when symmetric and odd
            q(ceil(qsize/2)) = 0; 
        end
    end
    %indices = indices(randperm(qsize)); %random order
    for i = indices
        %function to minimize:
        %negative of q(i) conntribution to mutual information
        objective = @(qi) -discreteMI(pX, getawgnqtransition(xsupport, [q(i-1); qi; q(i+1)], N));
        %deal with infinites in bounds
        lower = max(q(i-1), q(i)-abs(q(i))-endgap*10);
        upper = min(q(i+1), q(i)+abs(q(i))+endgap*10);
        if isSymmetric
            upper = min(upper, 0); %only do negatives when symmetric
        end
        %maximize MI over q(i)
        q(i) = fminbnd(objective, lower, upper, options);
        %assign symmetric point
        if isSymmetric
            q(qsize+1-i) = -q(i);
        end
    end
    newMI = discreteMI(pX, getawgnqtransition(xsupport, q, N));
    %check MI convergence
    if all(~((newMI-oldMI)>endgap)) %just stops if this produces NaN
        %if no boundary moved more that endgap
        break
    end
end
end

