clear
close all

%%  Generate random arcs

M = 3;
N = 4;

% Create arc array
centers = 3*complex(rand(M,N),rand(M,N));
radii = 2*rand(M,N)-1;
angInf = 2*pi*rand(M,N);
angSup = 2*pi*rand(M,N);
arcs = ciat.Arc(centers, radii, ciat.RealInterval(angInf,angSup));



%% Calculate area and verify with polygonal approximation

% Calculate arc areas
arcArea = abs(arcs.Area);

arcSamples = arcs.sample(1e4);
arcPolygons = ciat.PolygonalInterval(arcSamples);
arcPolyArea = arcPolygons.Area;

arcArea - arcPolyArea

% Plot
% figure;
clf;hold on;axis equal
arcs.plot('b');
arcPolygons.plot('g');