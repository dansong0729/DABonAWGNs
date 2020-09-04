N = 1;
E = 10;
m = 8;

[pX, xsupport, q, MI] = DABQ(N,E,m)
%% compute
N = 1;
m = 8;
dBs = [-20, -10, -5, 0, 3, 5, 7, 10, 12, 15, 17, 20];
%dBs = [-20, -10, -5, 0, 1:30];
inputPMFs = [];
xSupports = [];
qs = [];
MIs = [];
Es = [];
for db = dBs
    db
    E = 10.^(db/10);
    [pX, xsupport, q, MI] = DABQ(N,E,m);
    Es = [Es 10*log10(pX'*xsupport.^2)];
    MIs = [MIs MI];
    inputPMFs = [inputPMFs pX];
    xSupports = [xSupports xsupport];
    qs = [qs q];
end
MIs = MIs

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
    for j = 1:m
        %pmf points
        plot(dBs(i), xSupports(j, i),'ko', 'MarkerSize', 30*sqrt(inputPMFs(j,i))/2+1e-10, 'MarkerFaceColor', 'r')
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