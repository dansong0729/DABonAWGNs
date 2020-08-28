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