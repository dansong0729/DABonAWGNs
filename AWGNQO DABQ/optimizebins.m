function [q] = optimizebins(pX, xsupport, N)
%OPTIMIZEBINS Summary of this function goes here
%   INPUTS:
%   pX: column probability vector
%   xsupport: column vector of input distribution support points ASCENDING
%   ORDER
%   N: channel noise power

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
        %maximize MI over q(i) using old q values
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

