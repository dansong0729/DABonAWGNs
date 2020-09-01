N = 1;
E = 1;
m = 4;

[pX, xsupport, q, MI] = DABQ(N,E,m)
%%
N = 1;
m = 4;
dbs = [-20, -10, -5, 0, 3, 5, 7, 10, 12, 15, 17, 20];
MIs = []
for db = dbs
    db = db
    E = 10.^(db/10);
    [pX, xsupport, q, MI] = DABQ(N,E,m);
    MIs = [MIs MI];
end
MIs = MIs