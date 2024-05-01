% clear
% close all

%%  Generate random arcs

M = 2;
N = 1;

% Create arc array
centers = 3*complex(rand(M,N),rand(M,N));
radii = 2*rand(M,N)-1;
angInf = 2*pi*rand(M,N);
angSup = 2*pi*rand(M,N);
arcs = ciat.Arc(centers, radii, ciat.RealInterval(angInf,angSup));

% Plot
% figure;
clf;hold on;axis equal
arcs.plot('k');
arcs.plotMap(0,.1,'k');

%% Arc addition

% Add arcs
arcPlus = arcs(1) + arcs(2);

% Add a vertex to the arc
arcTranslate = arcs(1) + 1;

arcPlus.plot('b');
arcPlus.plotMap(0,.1,'b');
arcTranslate.plot('g');
arcTranslate.plotMap(0,.1,'g');

%%  Generate random edges

M = 1;
N = 1;

% Create edge array
startpoints = complex(rand(M,N),rand(M,N));
endpoints = complex(rand(M,N),rand(M,N));
edges = ciat.Edge(startpoints,endpoints);

% Plot
edges.plot('k');
edges.plotMap(0,.1,'k');


%% Edge addition

% Add an arc to the edge
edgePlus = edges(1) + arcs(1);

% Add a vertex to the edge
edgeTranslate = edges(1) + 1;

% Plot
edgePlus.plot('b')
edgePlus.plotMap(0,.1,'b');
edgeTranslate.plot('g')
edgeTranslate.plotMap(0,.1,'g');


