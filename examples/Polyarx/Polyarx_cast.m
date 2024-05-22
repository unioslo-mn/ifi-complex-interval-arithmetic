clear
% close all

%%

% Rectangle to Polyarx
rI = ciat.RectangularInterval(rand,rand,rand,rand);
rxI = ciat.PolyarxInterval(rI);

% Circle to Polyarx
cI = ciat.CircularInterval(complex(rand,rand),rand);
cxI = ciat.PolyarxInterval(cI);

% Polar to Polyarx
pI = ciat.PolarInterval(rand,rand,rand*2*pi,rand*2*pi);
pxI = ciat.PolyarxInterval(pI);

% Polyarcular to Polyarx
center1 = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
radius1 = [ 1 ; .1 ; -0.7 ; 0 ];
arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
arcs = ciat.Arc(center1,radius1,arcAngles);
aI = ciat.PolyarcularInterval( arcs );
axI = ciat.PolyarxInterval(aI);

%% 

% Polyarx to Rectangle
pxrI = ciat.RectangularInterval(pxI);

% Polyarx to Circle
rxcI = ciat.CircularInterval(rxI);

% Polyarx to Polar
cxpI = ciat.PolarInterval(cxI);

%% Plot
% figure;clf
subplot(2,2,1);cla;hold on;axis equal
rI.plot('c-')
rxI.plot('k--')
rxcI.plot('b-')
subplot(2,2,2);cla;hold on;axis equal
cI.plot('c-')
cxI.plot('k--')
cxpI.plot('b-')
subplot(2,2,3);cla;hold on;axis equal
pI.plot('c-')
pxI.plot('k--')
pxrI.plot('b-')
subplot(2,2,4);cla;hold on;axis equal
aI.plot('c-')
axI.plot('k--')
% pxrI.plot('b-')