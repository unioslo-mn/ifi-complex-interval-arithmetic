clear
% close all

%% Define a specific and a random polyarc 

% % Specify polyarc
center1 = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
radius1 = [ 1 ; .1 ; -0.7 ; 0 ];
arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
arcs1 = ciat.Arc(center1,radius1,arcAngles);
aI1 = ciat.PolyarcularInterval( arcs1 );

%% Define random points inside the polyarc rectangle
N = 1e3;
points = complex( rand(N,1) * aI1.Real.Width + aI1.Real.Infimum , ...
                  rand(N,1) * aI1.Imag.Width + aI1.Imag.Infimum);

pointIn(N,1) = false;
for n = 1:N
    pointIn(n) = aI1.isin(points(n));
end

%% Plot

% figure;
clf;hold on;axis equal
aI1.plot('k-');
plot(real(points(pointIn)),imag(points(pointIn)),'b.')
plot(real(points(~pointIn)),imag(points(~pointIn)),'r.')
