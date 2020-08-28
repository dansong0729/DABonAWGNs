function [q] = optimizebins(pX, xsupport, N)
%OPTIMIZEBINS Summary of this function goes here
%   INPUTS:
%   pX: column probability vector
%   xsupport: column vector of input distribution support points ASCENDING
%   ORDER
%   N: channel noise power

endgap = 1e-4; %1e-4 %threshold of q change for convergence
maxIter = 1000;

%initialize q to a guess, TODO: improve this (MLE?)
%use midpoints
q = (xsupport(2:end) + xsupport(1:end-1))./2;
q = [-Inf; q; Inf];


options = optimset('FunValCheck','on'); %for debugging

for iter = 1:maxIter
    nochange = true; %flag to decide termination
    qindices = 2:size(q)-1; %skip infinities; maybe randomize order later
    %TODO parfor? must wait until end to change q for loop independence
    %qnew = q
    for i = qindices
        old = q(i);
        %function to minimize
        objective = @(qi) -discreteMI(pX, getawgnqtransition(xsupport, [q(i-1); qi; q(i+1)], N));
        %deal with infinites in bounds
        lower = max(q(i-1), q(i)-abs(q(i))-endgap*10);
        upper = min(q(i+1), q(i)+abs(q(i))+endgap*10);
        %maximize MI over q(i) using old q values
        q(i) = fminbnd(objective, lower, upper, options); %TODO change to qnew
        if abs(q(i)-old)>endgap
            nochange = false;
        end
    end
    if nochange
        %if no boundary moved more that endgap
        break
    end
    %q = qnew
end
iter

end

