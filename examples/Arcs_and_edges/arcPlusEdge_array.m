clear
% close all

%% Generate a rand(M,N)om arc and edge pair

M = 3;
N = 5;

% Random arc
center = 3*complex(rand(M,N),rand(M,N));
radius = 2*rand(M,N)-1;
angInf = 2*pi*rand(M,N);
angSup = 2*pi*rand(M,N);
arcAngle = ciat.RealInterval( angInf , angSup );
arc = ciat.Arc( center , radius , arcAngle);

% Random edge
startpoint = complex(rand(M,N),rand(M,N));
endpoint = complex(rand(M,N),rand(M,N));
edge = ciat.Edge( startpoint, endpoint );

% Add arc and edge
arcPlusEdge = arc + edge;

%% Plot

% Normalize arcs
arcNorm = arc + (-arc.Center);
edgeNorm = edge + (-edge.Midpoint);
arcPlusEdgeNorm = arcPlusEdge + (-edge.Midpoint);

% figure;clf
% Plot arcs
subplot(1,2,1);cla;hold on;axis equal;title('Original')
arc.plot('k','lineWidth',2);
edge.plot('k','lineWidth',2);
arcPlusEdge.plot('b','lineWidth',2);
% Plot normalized arcs
subplot(1,2,2);cla;hold on;axis equal;title('Normalized')
arcNorm.plot('k','lineWidth',2)
arcNorm.plotMap(0,.1,'k');
edgeNorm.plot('k','lineWidth',2)
edgeNorm.plotMap(0,.1,'k');
arcPlusEdgeNorm.plot('b','lineWidth',2)
arcPlusEdgeNorm.plotMap(0,.1,'b');


