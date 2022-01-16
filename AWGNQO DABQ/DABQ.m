function [pX, xsupport, q, MI, s] = DABQ(N, E, supp_init, p_init, q_init, varargin)
%DABQ Implementation of DAB for quantized AWGN channel
%   INPUTS:
%   N: channel noise power
%   E: input power constraint
%   supp_init: col to initialize support points
%   p_init: col to initialize probabilities
%   q_init: col to initialize bins
%
%   OUTPUTS:
%   pX: col of optimal probabilities
%   xsupport: col of optimal support points
%   q: col of optimal bin thresholds, including +- INF endpoints
%   MI: mutual information achieved
%   s: langrange multiplier for power constraint


%settings
maxIter = 1e2;
ETolerance = 1e-5;
BAETolerance = 1e-5;
DABTolerance = 3e-5;
prob_thresh = 1e-6; %0 to disable remove redundant
sup_thresh = 1e-4; %^ same
DABiter = 1e2;
tol = 1e-4;

%init
xsupport = supp_init;
U = p_init;
pX = U;
%initialize q to a guess
q = q_init;

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
        [pX,xsupport] = remove_redundant(pX, xsupport, prob_thresh, sup_thresh);
        DAB_MI = MI;
        [pX, xsupport] = iteratesupport(pX, xsupport, q, N, E);
        [pX,xsupport] = remove_redundant(pX, xsupport, prob_thresh, sup_thresh);
    end
    %end discrete DAB
    %check for DABQ convergence
    if MI - oldMI < tol
        break
    end
    oldMI = MI;
end
end

