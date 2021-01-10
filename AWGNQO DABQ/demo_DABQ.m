N = 1;
E = 10^(3.1);
m = 4;
[supp_init, p_init] = equilattice(m, E);
%use midpoints
q_init = (supp_init(2:end) + supp_init(1:end-1))./2;
q_init = [-Inf; q_init; Inf];
[pX, xsupport, q, MI] = DABQ(N,E, supp_init, p_init);
%% compute
N = 1;
m = 8;
% dBs = [-20, -10, -5, 0, 3, 5, 7, 10, 12, 15, 17, 20];
dBs = 1:20;
inputPMFs = [];
xSupports = [];
qs = [];
MIs = [];
Es = [];
ss = [];
for db = dBs
    db
    E = 10.^(db/10);
    [supp_init, p_init] = equilattice(m, E);
    q_init = (supp_init(2:end) + supp_init(1:end-1))./2;
    q_init = [-Inf; q_init; Inf];
    [pX, xsupport, q, MI, s] = DABQ(N,E, supp_init, p_init, q_init);
    Es = [Es 10*log10(pX'*xsupport.^2)];
    MIs = [MIs MI];
    inputPMFs = [inputPMFs symmetric_pad(pX, m)];
    xSupports = [xSupports symmetric_pad(xsupport, m)];
    qs = [qs q];
    ss = [ss s];
end

%% plot

%plot rates
close all
figure(1)
grid on
plot(dBs, MIs)
xlabel('SNR (dB)')
ylabel('Rate (bits)')
ylim([min(MIs)-.1,log2(m)+.1])

%plot distributions
figure(2)
xlim([min(dBs)-1,max(dBs)+1])
hold on
for i = 1:length(dBs)
    [pX_merged, xsupport_merged] = remove_redundant(inputPMFs(:, i), xSupports(:, i), 1e-5, 1e-2);
    for j = 1:length(pX_merged)
        %pmf points
        plot(dBs(i), xsupport_merged(j),'ko', 'MarkerSize', 30*sqrt(pX_merged(j))/2+1e-10, 'MarkerFaceColor', 'r')
    end
end
%plot bin thresholds
for j=1:m
    %make "staircase"
    temp = [qs(j, :);qs(j, :)];
    plot(sort([dBs+.5, dBs-.5]), temp(:)')
end
hold off
title('DABQ Optimized input PMFs with cardinality ' + string(m))
xlabel('SNR (dB)')
ylabel('Mass Point Locations')
drawnow

%% stemplot
idx = 30;
figure
[pX_merged, xsupport_merged] = remove_redundant(inputPMFs(:, idx), xSupports(:, idx), 1e-5, 1e-2);
stem(xsupport_merged, pX_merged)
hold on
stem(qs(2:end-1, idx), ones(m-1,1)*max(pX_merged)*1.05, 'Marker', 'none')
ylim ([0,max(pX_merged)*1.05])