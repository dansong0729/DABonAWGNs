function[p_star, MI, EE, s] = BAE_search(pX, xsupport, q, N, E, s_init, ETolerance, BAETolerance)
%BAE_SEARCH Compute compacity of AWGNQO channel with input power constraint
%   Searches over Lagrange multiplier s passed to constrained
%   Blahut-Arimoto algorithm (BAE_search.m), and outputs resulting
%   capacity, input distribution, and energy achieved.
%   
%   INPUTS:
%   pX: initial guess for input pmf, as column probability vector
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
ej = xsupport.^2; %power cost of each point

p_star = pX;

%wrap BAE_discrete into function s |-> EE-E to pass into fzero
    function err = wrapper(ess)
        %average guess to avoid 0s in initializing (can't ever recover from
        %being zero in BAE)
        guess = (pX+7*p_star)/8; %TODO better guess (average multiple pX?)
        %intentional side effect of writing to pX and MI
        %output pX will serve as initial guess for next evaulation
        [p_star, MI, tempE] = BAE_discrete(Q, guess, ej, ess, BAETolerance);
        err = tempE-E;
    end

%OutputFcn to determine halting
    function stop = outfun(~, optimValues, ~, varargin)
        %can turn on plotting here
        %optimplotfval(x, optimValues, state);
        %halt only when under power constraint
        stop = false;
        %fval is empty on the first iteration (wtf?!) so check is needed
        if ~isempty(optimValues.fval)
            stop = (optimValues.fval<0)&(optimValues.fval>-ETolerance);
        end
    end

%check if constraint too weak
s = 0;
EE = wrapper(s) + E;
if EE < E
    fprintf('BAE finds p(x) with E = %f, constraint irrelevant\n', EE);
    return
end

%get carried by matlab
%outfun handles termination, so TolX is made as small as possible (eps)
options = optimset('TolX', eps, 'OutputFcn', {@outfun},...
                    'FunValCheck', 'on', 'Display', 'off');
[s, ~, ~, ~] = fzero(@wrapper, s_init, options);
%paranoia; in case fzero's last call to wrapper wasn't at s
EE = wrapper(s);
EE = EE + E;

end

