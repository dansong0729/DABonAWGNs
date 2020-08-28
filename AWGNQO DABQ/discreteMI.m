function MI = discreteMI(pX,Q)
%discreteMI: Computes Mutual information given input pmf and transistion
%matrix.
%   INPUTS:
%   pX: column probability vector
%   Q: a matrix with probability vector columns
%   OUTPUTS: Mutual information, in bits

%By passing only some rows of transition matrix as Q, the contribution of
%the corresponding output points is computed.

%compute H(Y)
pY = Q*pX;
mask = (pY ==0);
HY = reallog(pY)./log(2);
HY(mask) = 0; %kill -Inf since we want 0log0=0
HY = -pY.'*HY;

%compute H(Y|X);
mask = (Q ==0);
HYX = reallog(Q)./log(2);
HYX(mask) = 0; %kill -Inf since we want 0log0=0
HYX = sum(Q.*HYX, 1);
HYX = -HYX*pX;

MI = HY - HYX;
end

