function [pX, xsupport, q, MI] = DABQ(N, E, m)
%DABQ Summary of this function goes here
%   Detailed explanation goes here

%settings
maxIter = 1e3;
ETolerance = 1e-5;
BAETolerance = 1e-5;
DABTolerance = 3e-5;
DABiter = 1e2;
tol = 1e-4;

%init
[xsupport, U] = equilattice(m, E);
pX = U;
s = 0.1;
oldMI = 0;

for i = 1:maxIter
    q = optimizebins(pX, xsupport, N, true);
    DAB_MI = oldMI;
    for j = 1:DABiter
        [pX, MI, ~, s] = BAE_search((pX+U)/2, xsupport, q, N, E, s, ETolerance, BAETolerance);
        if MI - DAB_MI < DABTolerance
            break
        end
        DAB_MI = MI;
        [xsupport, pX] = iteratesupport(pX, xsupport, q, N, E);
    end
    if MI - oldMI < tol
        break
    end
    oldMI = MI;
end
end

