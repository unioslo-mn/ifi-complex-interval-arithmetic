clear
% close all

%%  Generate random arcs

M = 4;
N = 3;

% Create arc array
centers = 3*complex(rand(M,N),rand(M,N));
radii = 2*rand(M,N)-1;
angInf = 2*pi*rand(M,N);
angSup = 2*pi*rand(M,N);
arcs = ciat.Arc(centers, radii, ciat.RealInterval(angInf,angSup));

%% Cast arcs
rI = ciat.RectangularInterval(arcs);
% cI = ciat.CircularInterval(arcs);
pI = ciat.PolarInterval(arcs);


% Plot
% figure;clf
cla;hold on;axis equal
arcs.plot('k','linewidth',2);
% rI.plot('g');
% cI.plot('r');
pI.plot('b');