clear
% close all

%% Define a specific and a random polyarc 

% % Specify polyarc
center1 = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
radius1 = [ 1 ; .1 ; -0.7 ; 0 ];
arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
arcs = ciat.Arc(center1,radius1,arcAngles);
aI = ciat.PolyarcularInterval( arcs );

%% Sample polyarcular interval and calculat polygonal area

aIsamp = sample(aI,1e2);
aIsampArea = polyarea(real(aIsamp{:}),imag(aIsamp{:}));


%%

aI.Area - aIsampArea


% Plot
% figure;
clf;hold on;axis equal
aI.plot('g-')
% aIpoly.plot('r:')
plot(real(aIsamp{:}),imag(aIsamp{:}),'b--')
