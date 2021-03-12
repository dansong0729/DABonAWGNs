function [Q] = get_transition(xsupp, y_range, pYX_func)
%GET_TRANSITION Convenience function to get transition matrix from function
%that only takes scalars.
%   Outputs length(y_range) by length(xsupp) matrix, populated with
%   pYX_func evaluated at the corresponding x and y values.
Q = zeros(length(y_range), length(xsupp));
for i=1:length(y_range)
    for j=1:length(xsupp)
        Q(i,j) = pYX_func(y_range(i), xsupp(j));
    end
end
end

