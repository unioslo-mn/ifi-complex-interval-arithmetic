clear
% close all

%% Define a specific and a random polyarc 

% Specific polyarc
center1 = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
radius1 = [ 1 ; .1 ; -0.7 ; 0 ];
arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
arcs1 = ciat.Arc(center1,radius1,arcAngles);
aI1 = ciat.PolyarcularInterval( arcs1 );

%%
% Polyarc from a random polar interval
pI = ciat.PolarInterval(rand,rand,rand,rand);
paI = ciat.PolyarcularInterval(pI);

%%

% Add the two intervals
% aIsum = aI1 + paI;

%%
% figure;
clf;hold on;axis equal
aI1.plot('k-');
aI1.plotGaussMap(.1,'g');
aI1.Vertices.plotGaussMap(.1,'b')
paI.plot('k-');
paI.plotGaussMap(.1,'k');