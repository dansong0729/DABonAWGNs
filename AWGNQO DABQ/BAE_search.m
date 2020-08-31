function[p_star, MI, EE, s] = BAE_search(pX, xSupport, q, N, E, s_init, ETolerance, BAETolerance)
%BAE_SEARCH Compute compacity of AWGNQO channel with input power constraint
%   Searches over Lagrange multiplier s passed to constrained
%   Blahut-Arimoto algorithm (BAE_search.m), and outputs resulting
%   capacity, input distribution, and energy achieved.
%   
%   INPUTS:
%   init: initial guess for input pmf, as column probability vector
%   xsupport: column vector of input pmf support points, same order as init
%   q: bin thresholds in ascending order, as column vector
%   N: channel noise power
%   E: input power constraint
%   s_init: initial guess for s
%   ETolerance: tolerance on achieving meeting power constraint
%   BAETolerance: tolerance for Blahut-Arimoto termination, on MI
%   
%   OUTPUTS:
%   pX: Capacity achieving input pmf, as column probability vector
%   MI: Capacity, in bits
%   EE: Power of input distribution defined by pX and xsupport
%   s: Lagrange multiplier parameter value that met power constraint E

%compute constants
Q = getawgnqtransition(xsupport, q, N); %channel transition matrix
ej = xSupport.^2; %power cost of each point

p_star = pX;

%wrap BAE_discrete into function s |-> EE-E to pass into fzero
    function err = wrapper(ess)
        guess = p_star; %TODO better guess (average multiple pX?)
        %intentional side effect of writing to pX and MI
        %output pX will serve as initial guess for next evaulation
        [p_star, MI, tempE] = BAE_discrete(Q, guess, ej, ess, BAETolerance);
        err = tempE-E+ETolerance; %err on the side of not exceeding E
    end

%check if constraint too weak
s = 0;
EE = wrapper(s) + E - ETolerance;
if EE < E
    fprintf('BAE finds p(x) with E = %f, constraint irrelevant', EE);
    return
end

%get carried by matlab
options = optimset('TolX', ETolerance); %can turn on plotting here
[s, EE, ~, ~] = fzero(@wrapper, s_init, options); %TODO better intialization
%EE = wrapper(s); %just in case fzero's last call to wrapper wasn't at s
EE = EE + E - ETolerance;

end

