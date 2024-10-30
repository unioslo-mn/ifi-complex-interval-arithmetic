clear
% close all

%% Generate random rectangular intervals
N = 1;
M = 1;
rI = ciat.RectangularInterval(randn(N,M),randn(N,M),randn(N,M),randn(N,M));

%% Cast to polar

rpI = ciat.PolarInterval(rI);

%% Plot
% figure;clf
cla;hold on;axis equal
plot(0,0,'+')
rI.plot('k')
rpI.plot('b')