clear
% close all

%% Define a specific and a random polyarc 

% % Specify polyarc
% center1 = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
% radius1 = [ 1 ; .1 ; -0.7 ; 0 ];
% arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
% arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
% arcs1 = ciat.Arc(center1,radius1,arcAngles);
% aI1 = ciat.PolyarcularInterval( arcs1 );

% Specify polyarc
center1 = [ 3+0.5i ;  4+2.5i ; 2+1.5i];
radius1 = [ 1 ; -0.7 ; 0 ];
arcAngles = [ [-.7,-.5] ; [.3,.5] ; [0 0] ] *pi;
arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
arcs1 = ciat.Arc(center1,radius1,arcAngles);
aI1 = ciat.PolyarcularInterval( arcs1 );

%%
% Polyarc from a random polar interval
pI = ciat.PolarInterval(4,5,0.3*pi,0.4*pi);
aI2 = ciat.PolyarcularInterval(pI);

%% Demonstrate detailed method of the plus function

aI1Arcs = [aI1.Arcs ; aI1.Vertices];
aI2Arcs = [aI2.Arcs ; aI2.Vertices];

aI1Edges = aI1.Edges;
aI2Edges = aI2.Edges;

aI3ArcPlusArc = aI1Arcs + aI2Arcs.';
aI3EdgePlusEdge = aI1Edges + aI2Edges.';
aI3ArcPlusEdge = aI1Arcs + aI2Edges.';
aI3EdgePlusArc = aI1Edges + aI2Arcs.';

% Select non-empty and non-vertex arcs and edges
aI3Arc = aI3ArcPlusArc(~isnan(aI3ArcPlusArc));
aI3Edge = [ aI3EdgePlusEdge(~isnan(aI3EdgePlusEdge));...
            aI3ArcPlusEdge(~isnan(aI3ArcPlusEdge));...
            aI3EdgePlusArc(~isnan(aI3EdgePlusArc))];
aI3Arc = aI3Arc(aI3Arc.Length>0);
aI3Edge = aI3Edge(aI3Edge.Length>0);

% Convert arcs and edges to rectangular intervals and check intersections
aI3ArcBox = ciat.RectangularInterval(aI3Arc);
aI3EdgeBox = ciat.RectangularInterval(aI3Edge);
capBox = cap( [aI3ArcBox;aI3EdgeBox] , [aI3ArcBox;aI3EdgeBox].');
capMask = tril(~isnan(capBox) & capBox.Area>0,-1);
[M,N] = find(capMask);

%% Run the plus function

% Add the two intervals
% aIsum = aI1 + aI2;


%%
% figure;
clf;hold on;axis equal
aI1.plot('k-');
aI1.plotGaussMap(.1,'k');
aI2.plot('k-');
aI2.plotGaussMap(.1,'k');

aI3Arc.plot('r')
aI3Edge.plot('b')
aI3ArcBox.plot('r--')
aI3EdgeBox.plot('b--')

aI3ArcMid = arcIn.Center + arcIn.Radius .* exp(1j*arcIn.ArcAngle.Midpoint);
text(aI3ArcBox.Real.Midpoint,aI3ArcBox.Imag.Midpoint,string(1:length(aI3Arc))');
text(real(aI3Edge.Midpoint),imag(aI3Edge.Midpoint),string(1:length(aI3Edge)));