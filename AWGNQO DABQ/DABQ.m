function [pX, xsupport, q, MI, s] = DABQ(N, E, m)
%DABQ Summary of this function goes here
%   Detailed explanation goes here

%settings
maxIter = 1e2;
ETolerance = 1e-5;
BAETolerance = 1e-5;
DABTolerance = 3e-5;
DABiter = 1e2;
tol = 1e-4;

%init
[xsupport, U] = equilattice(m, E);
pX = U;
%initialize q to a guess, TODO: improve this (MLE?)
%use midpoints
q = (xsupport(2:end) + xsupport(1:end-1))./2;
q = [-Inf; q; Inf];

s = 0.1;
oldMI = 0;

for i = 1:maxIter
    q = optimizebins(pX, xsupport, q, N, true);
    DAB_MI = oldMI;
    %discrete DAB
    for j = 1:DABiter
        %cardinality may have changed
        [~, U] = equilattice(length(pX), E); %uniform spacing
        %Blahut-Arimoto
        [pX, MI, ~, s] = BAE_search(U, xsupport, q, N, E, s, ETolerance, BAETolerance);
        %check for DAB convergence
        if MI - DAB_MI < DABTolerance
            break
        end
        [pX,xsupport] = remove_redundant(pX, xsupport, 1e-10, 1e-3);
        DAB_MI = MI;
        [pX, xsupport] = iteratesupport(pX, xsupport, q, N, E);
        [pX,xsupport] = remove_redundant(pX, xsupport, 1e-10, 1e-3);
    end
    %end discrete DAB
    %check for DABQ convergence
    if MI - oldMI < tol
        break
    end
    oldMI = MI;
end
end

