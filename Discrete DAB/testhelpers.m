%% check MI, rel entropy with BSC
pX = [.5; .5];
Q = [1 0; 0 1];
change = [-1  1; 1 -1]; 

for i = 0:10
    p = 0.1*i
    Q_iter = Q+p*change;
    MI = discrete_MI(pX, Q_iter)
    diff = MI - discrete_rel_entropy(Q_iter.*pX.', pX.'.*(Q_iter*pX))
end
%% check MI, rel entropy with binary erasure channel
pX = [.5; .5];
Q = [1 0; 0 0; 0 1];
change = [-1  0; 1 1; 0 -1]; 

for i = 0:10
    p = 0.1*i
    Q_iter = Q+p*change;
    MI = discrete_MI(pX, Q_iter)
    diff = MI - discrete_rel_entropy(Q_iter.*pX.', pX.'.*(Q_iter*pX))
end
%% random pX into MI, just to make sure it doesn't error
for i = 0:1e4
    n = randi(40, 1);
    Q = eye(n);
    pX = rand(n, 1);
    pX = pX./sum(pX);
    assert(abs(pX.'*log2(pX) + discrete_MI(pX, Q))<1e-10);
end
"success"
%% random pX into rel entropy, make sure close to 0
for i = 0:1e4
    n = randi(40, 1);
    pX = rand(n, 1);
    pX = pX./sum(pX);
    assert(discrete_rel_entropy(pX, pX)<1e-10);
end
"success"