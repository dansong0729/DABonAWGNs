function MI = discrete_MI(pX,Q)
%discreteMI: Computes Mutual information given input pmf and transistion
%matrix.
%   INPUTS:
%   pX: Input pmf, as column probability vector
%   Q: A matrix with probability vector columns; Q_y,x = p_Y|X(y, x)
%   OUTPUT:
%   MI: Mutual information, in bits

%By passing only some rows of transition matrix as Q, the contribution of
%the corresponding output points is computed.

pY = Q*pX;
%compute log(pY|X / pY)
quotient = reallog(Q./pY); %reallog warns on invalid inputs
%locate NaN issues
Qmask = (Q==0);
ymask = (pY==0);
%set 0log(0) and 0log(1/0) to 0
quotient(Qmask) = 0;
quotient(ymask, :) = 0;
vec = sum(Q.*quotient, 1); %sum over y
MI = vec*pX; %sum over x
MI = MI/log(2); %from nats to bits
end

