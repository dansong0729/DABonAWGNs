function [p_star, MI, E] = BAE_discrete(Q, init, ej, s, tol)
%DISCRETE_BA Computes capacity and optimal input pmf of constrained DMC
%   Straight outta the Blahut '72 paper.
%   INPUTs:
%   Q: matrix defining p(y|x), with probability vector columns
%       Q_ij = p(y_i|x_j)
%   init: initial guess for input pmf pj, as column probability vector
%   ej: "cost" (e.g. sq. amplitude) of support points, same order as init 
%   s: parameter controlling power constraint, set s=0 for unconstrained
%   tol: max error for termination
%
%   OUTPUTS:
%   pX: Capacity achieving input pmf (pj in paper), as column probability vector
%   MI: Capacity, in bits (CE, in paper)
%   E: Cost of pX

maxIter = 1e6;

Qmask = (Q==0); %log(0) causes issues
pj = init;
for i=1:maxIter
    %compute matrix of summands for cj formula
    cj = log2(Q./(Q*pj));
    cj(Qmask) = 0; %0*log(0) conventionally taken as 0
    cj(pj==0, :) = 0; %same for where pj causes Inf
    cj = Q.*cj;
    
    cj = sum(cj, 1)'-s*ej; %turn cj into column before subtraction
    cj = 2.^(cj);
    
    denom = pj.'*cj; %reused scalar value
    
    %check convergence
    IL = log2(denom);
    IU = log2(max(cj));
    if IU - IL < tol
        E = pj.'*ej;
        MI = s*E+IL;
        p_star = pj;
        return
    end
    
    %update pj
    pj = pj.*cj/denom;
end
fprintf('maxIter reached in Blahut Arimoto\n')
E = pj.'*ej;
MI = s*E+IL;
p_star = pj;
end

