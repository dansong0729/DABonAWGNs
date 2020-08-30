%check MI using BSC
pX = [.5; .5];
Q = [1 0; 0 1];
change = [-1  1; 1 -1]; 

for i = 0:10
    p = 0.1*i
    discreteMI(pX, Q+p*change)
end
%% use erasure to test MI
pX = [.5; .5];
Q = [1 0; 0 0; 0 1];
change = [-1  0; 1 1; 0 -1]; 

for i = 0:10
    p = 0.1*i
    discreteMI(pX, Q+p*change)
end
%% random pX into MI, just to make sure it doesn't error
for i = 0:100
    n = randi(20, 1);
    Q = eye(n);
    pX = rand(n, 1);
    pX = pX./sum(pX);
    assert(abs(pX.'*log2(pX) + discreteMI(pX, Q))<1e-10);
end
"success"
%% awgntransition fuzzing, just make sure everything looks right
q = [-6:6]';
q = [-Inf; q; Inf];
n = randi(5, 1)+1
xsupport = rand(n,1)*4-2;
pYX = getawgnqtransition(xsupport, q, 3);
assert(all(abs(sum(pYX, 1)-1)==0)); %probabilities add to 1
assert(all(pYX>=0, 'all'));

figure
hold on
for col = pYX
    plot(col)
end
%% test unconstrained BA with symmetric channels

%BSC + erasure
n = 2;
pX = rand(n, 1);
pX = pX./sum(pX)
%channel parameters
Q = [1 0; 0 0; 0 1];
alpha = rand(1);
p = rand(1)*alpha;
alpha = alpha-p;

Q = Q + alpha*[-1 0; 1 1; 0 -1] + p*[-1 1; 0 0; 1 -1]

%make sure it's close to uniform
[pX, MI, E] = BAE_discrete(Q, pX, [1;1], 0, 1e-4);
pX
MI
discreteMI(pX, Q) %correct MI

%% random af channels, random constraints
n = randi(20, 1)+1;
pX = rand(n, 1);
pX = pX./sum(pX);
k = n+randi(ceil(n/2),1) - floor(n/4);
Q = rand(k, n);
Q = Q./sum(Q);
ej = randi(10, [n 1]) - 1;

ss = logspace(-2,1, 6);
for s=ss
    s = s
    [pX, MI, E] = BAE_discrete(Q, pX, ej, s, 1e-5);
    MI = MI
    discreteMI(pX, Q); %correct MI
end