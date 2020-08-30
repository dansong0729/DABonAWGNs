% should output very close to [-Inf; 0; Inf], as this is the m=2 case in
% Madhow
xsupport = [-1; 1];
pX = [.5; .5];
N = rand(1)*5;
optimizebins(pX, xsupport, N, true)
%% fuzzing, trying to fix negative probability
for i = 1:200
    n = randi(5, 1)+1;
    xsupport = sort(rand(n,1)*4-2);
    pX = rand(n, 1);
    pX = pX./sum(pX);
    N = (rand(1)+.1)*2;

    optimizebins(pX, xsupport, N, false);
end
"end"
%% test symmetric
for i = 1:200
    n = randi(10, 1)+1;
    xsupport = sort(rand(floor(n/2),1)*4);
    if mod(n,2)
        xsupport = [sort(-xsupport); 0; xsupport];
    else
        xsupport = [sort(-xsupport); xsupport];
    end
    
    pX = rand(n, 1);
    pX = pX./sum(pX);
    N = (rand(1)+.1)*2;
    
    optimizebins(pX, xsupport, N, true);
end
"end"