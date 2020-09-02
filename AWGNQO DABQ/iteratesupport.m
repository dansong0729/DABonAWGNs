function [x_star] = iteratesupport(pX, xsupport, q, N, E)
%ITERATESUPPORT Perform round of updating support points
%   Described by step 4 in paper. Overall structure is copied from
%   optimizebins; innermost loop logic same as FindMaximumI in PC-AWGN.
%   INPUTS:
%   pX: input pmf, as column probability vector
%   xsupport: input pmf support points in ascending order, as column vector
%   q: bin thresholds in ascending order, as  column vector
%   N: channel noise power
%   E: input power constraint
%   isSymmetric: set true to exploit xsupport symmetry about 0
%   
%   OUTPUTs:
%   x_star: vector of new support points in ascending order, as column vector
%
%TODO: Characterize and fix precision issues.

rounds = 1; %makes sense to be higher since BAE_search is so slow

%initialize
x_star = xsupport;
m = size(xsupport,1); %constellation cardinality

TolX = 1e-4; %default 1e-4, put into options otherwise
options = optimset('FunValCheck','on'); %for safety

for iter = 1:rounds
    %only do half
    indices = 1:floor(m/2);
    if mod(m,2)
        %middle element is 0 when symmetric and odd
        x_star(ceil(m/2)) = 0; 
    end
    for i = indices
        if i ~= 1
            %general case
            
            %compute bounds
            e = pX.*xsupport.^2;
            if m == 4 
                P_in = 0;
                E_in = 0;
            else
                P_in = sum(pX([i+1:m-i]));
                E_in = sum(e([i+1:m-i]));
            end
            xSearchLim = max(xsupport(i-1),-sqrt((E-E_in)/(1-P_in)));
            %jank fix
            if xSearchLim > xsupport(i)
                continue
            end
            %end jank
            %objective to minimize: negative mutual information with
            %rescaled xsupport 
            fun = @(Qx_star) -Find_symmetric_I_x_star_out(xsupport, pX, i, Qx_star, q, N, E);
            [x_star(i), ~, ~] = fminbnd(fun, xSearchLim, xsupport(i)+TolX, options);
        else
            %outmost point special case
            xSearchLim = q(2) - 5*sqrt(N); %5 std out from leftmost boundary
            %jank fix
            xSearchLim = min(xSearchLim, x_star(i) - 5*sqrt(N));
            %end jank
            if m == 3
                pX([1,3]) = pX([1,3])-1e-3;
                pX(2) = pX(2)+2e-3;
            end
            %objective to minimize: negative mutual information with
            %rescaled xsupport 
            fun = @(Qx_star) -Find_symmetric_I_x_star_in(xsupport, pX, i, Qx_star, q, N, E);
            [x_star(i), ~, ~] = fminbnd(fun, xSearchLim, xsupport(i)+TolX, options);
        end
        %move mirror point
        x_star(m+1-i) = -x_star(i);
    end
end
end
