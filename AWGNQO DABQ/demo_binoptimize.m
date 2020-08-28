% should output very close to [-Inf; 0; Inf], as this is the m=2 case in
% Madhow
xsupport = [-1; 1];
pX = [.5; .5];
N = rand(1)*5;
optimizebins(pX, xsupport, N)
%% fuzzing, trying to fix negative probability
for i = 1:20
    n = randi(5, 1)+1;
    xsupport = sort(rand(n,1)*4-2);
    pX = rand(n, 1);
    pX = pX./sum(pX);
    N = rand(1)*2;

    optimizebins(pX, xsupport, N)
end