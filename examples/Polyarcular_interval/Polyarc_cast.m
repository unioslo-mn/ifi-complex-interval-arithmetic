clear
% close all

%% Cast a polyarcular interval to other interval types

% % Specify polyarc
center = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
radius = [ 1 ; .1 ; -0.7 ; 0 ];
arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
arcs = ciat.Arc(center,radius,arcAngles);
aI = ciat.PolyarcularInterval( arcs );

% Cast 
arI = ciat.RectangularInterval(aI);
apI = ciat.PolarInterval(aI);
aoI = ciat.CircularInterval(aI);

%% Cast bounding intervals back to polyarc
araI = ciat.PolyarcularInterval(arI);
apaI = ciat.PolyarcularInterval(apI);
aoaI = ciat.PolyarcularInterval(aoI);


%% Plot
% figure;
cla;hold on;axis equal
aI.plot('k','linewidth',2);
arI.plot('b-','linewidth',2);
apI.plot('r-','linewidth',2);
aoI.plot('c-','linewidth',2);
araI.plot('y:','linewidth',2);
apaI.plot('y:','linewidth',2);
aoaI.plot('y:','linewidth',2);