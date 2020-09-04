N = 1;
E = 10^(-.3);
m = 8;

[pX, xsupport, q, MI] = DABQ(N,E,m)
%%
N = 1;
m = 8;
dbs = [-20, -10, -5, 0, 3, 5, 7, 10, 12, 15, 17, 20];
%dbs = -20:20;
pXs = [];
supports = [];
qs = [];
MIs = [];
Es = [];
for db = dbs
    E = 10.^(db/10);
    [pX, xsupport, q, MI] = DABQ(N,E,m);
    Es = [Es 10*log10(pX'*xsupport.^2)];
    MIs = [MIs MI];
    pXs = [pXs pX];
    supports = [supports xsupport];
    qs = [qs q];
end
MIs = MIs