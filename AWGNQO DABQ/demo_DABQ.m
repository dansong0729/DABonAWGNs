N = 1;
<<<<<<< HEAD
% E = 10^(2.6);
% m = 32;
E = 10^(-.3);
m = 11;

[pX, xsupport, q, MI] = DABQ(N,E,m);
%% compute
N = 1;
m = 11;
%dBs = [-20, -10, -5, 0, 3, 5, 7, 10, 12, 15, 17, 20];
dBs = -10:5;
=======
E = 10^(3.1);
m = 4;

[pX, xsupport, q, MI] = DABQ(N,E,m)
%% compute
N = 1;
m = 8;
% dBs = [-20, -10, -5, 0, 3, 5, 7, 10, 12, 15, 17, 20];
dBs = 1:35;
>>>>>>> b9499652121f2ecd4599c9ffa9620438ce0e48a8
inputPMFs = [];
xSupports = [];
qs = [];
MIs = [];
Es = [];
ss = [];
for db = dBs
    db
    E = 10.^(db/10);
    [pX, xsupport, q, MI, s] = DABQ(N,E,m);
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
<<<<<<< HEAD
    [pX_merged, xsupport_merged] = remove_redundant(inputPMFs(:, i), xSupports(:, i), 1e-4, 5e-2);
=======
    [pX_merged, xsupport_merged] = remove_redundant(inputPMFs(:, i), xSupports(:, i), 1e-5, 1e-2);
>>>>>>> b9499652121f2ecd4599c9ffa9620438ce0e48a8
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