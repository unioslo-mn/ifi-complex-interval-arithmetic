clear
close all

%% Define two arcs

arc1 = ciat.Arc(1+1i,1,ciat.RealInterval(1,2));
arc2 = ciat.Arc(2+2i,2,ciat.RealInterval(3,4)); 

% Plot
figure;clf;hold on;axis equal
arc1.plot;
arc2.plot;

%%  Arc array behavioour

M = 2;
N = 3;

% Create arc array
arcs = ciat.Arc(10*complex(rand(M,N),rand(M,N)), ...
                rand(M,N), ...
                2*pi*ciat.RealInterval(rand(M,N),rand(M,N)));

% Plot
arcs.plot;

%% Arc addition

% Add arcs
arcSum = arc1 + arc2;

% Add a vertex to the arc
arcTranslate = arc1 + 1;

arcSum.plot
arcTranslate.plot

%% Define two edges

edge1 = ciat.Edge(1+2i, 2+1i);
edge2 = ciat.Edge(0+3i, 2+0i);

edge1.plot
edge2.plot

%% Define an edge array

edgeArr = ciat.Edge(complex(rand(3,4),rand(3,4)), ...
                    complex(rand(3,4),rand(3,4)));

edgeArr.plot