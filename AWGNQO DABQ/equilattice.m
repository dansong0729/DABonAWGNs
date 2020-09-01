function [xi, wi] = equilattice(m, E)
    wi = ones([1,m])'./m; %equiprobable
    %compute equally spaced xi
    spacing = sqrt(3/(m^2 - 1)); %same as delta m in paper
    lower = spacing*(1-m);
    upper = spacing*(m-1);
    xi = linspace(lower, upper, m)'*sqrt(E);
end