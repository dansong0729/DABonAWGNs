function [q] = optimizebins(pX, xsupport, N)
%OPTIMIZEBINS Find optimal bins for AWGNQO channel given input pmf
%   Performs alternating optimization on each bin threshold using fminbnd
%   and terminates on sufficiently small mutual information gain between
%   iterations.
%   
%   INPUTS:
%   pX: input pmf, as column probability vector
%   xsupport: input pmf support points in ASCENDING order, as column vector
%   N: channel noise power
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
q = (xsupport(2:end) + xsupport(1:end-1))./2;
q = [-Inf; q; Inf];


options = optimset('FunValCheck','on'); %for debugging
% %=====DEBUG=====
% qhist = [];
% %=====END=====
for iter = 1:maxIter
    oldMI = discreteMI(pX, getawgnqtransition(xsupport, q, N));
    for i = 2:size(q,1)-1
        %function to minimize:
        %negative of q(i) conntribution to mutual information
        objective = @(qi) -discreteMI(pX, getawgnqtransition(xsupport, [q(i-1); qi; q(i+1)], N));
        %deal with infinites in bounds
        lower = max(q(i-1), q(i)-abs(q(i))-endgap*10);
        upper = min(q(i+1), q(i)+abs(q(i))+endgap*10);
        %maximize MI over q(i)
        q(i) = fminbnd(objective, lower, upper, options);
    end
    newMI = discreteMI(pX, getawgnqtransition(xsupport, q, N));
    %check MI convergence
    if all(~((newMI-oldMI)>endgap)) %just stops if this produces NaN
        %if no boundary moved more that endgap
        break
    end
%     %=====DEBUG=====
%     qhist = [qhist q];
%     %=====END=====
end
% %=====DEBUG=====
% if iter > 100
%     iter
%     pX = pX
%     xsupport = xsupport
%     N = N
%     figure
%     hold on
%     for qi = qhist'
%         plot([1:iter]', qi)
%     end
%     plot(zeros(size(xsupport), 1), xsupport)
% end
% %=====END=====

end

