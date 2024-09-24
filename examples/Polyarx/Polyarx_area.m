clear
% close all

%% Define a specific and a random polyarc 

% % % Specify polyarc
% center1 = [ 3+0.5i ; 3.1 - 0.5i ; 4+2.5i ; 2+1.5i];
% radius1 = [ 1 ; .1 ; -0.7 ; 0 ];
% arcAngles = [ [-.7,-.5] ; [-1,-.3] ; [.3,.5] ; [0 0] ] *pi;
% arcAngles = ciat.RealInterval(arcAngles(:,1),arcAngles(:,2));
% arcs = ciat.Arc(center1,radius1,arcAngles);
% aI = ciat.PolyarcularInterval( arcs );
% xI = ciat.PolyarxInterval(aI);

% Specify random polyarc
N = 3;
absMin = rand(N,1);
absMax = absMin + rand(N,1)/10;
angMin = 2*pi*rand(N,1);
angMax = angMin + 2*pi*rand(N,1)/20;
pI = ciat.PolarInterval(absMin,absMax,angMin,angMax);
paI = ciat.PolyarcularInterval(pI);
pxI = ciat.PolyarxInterval(pI);
aI = sum(paI);
xI = sum(pxI);


%% Sample polyarcular interval and calculat polygonal area

aIsamp = sample(aI,1e2);
aIsampArea = polyarea(real(aIsamp{:}),imag(aIsamp{:}));

xIsamp = sample(xI,1e2);
xIsampArea = polyarea(real(xIsamp{:}),imag(xIsamp{:}));


%%

sprintf(['Area of the polyarc: %0.3f (approximation: %0.3f) \n',...
        'Area of the polyarx: %0.3f (approximation: %0.3f) \n',...
        'Area differences respectively: %0.3e and %0.3e\n'],...
        aI.Area, aIsampArea,...
        xI.Area, xIsampArea,...
        aI.Area - aIsampArea , xI.Area - xIsampArea)


% Plot
% figure;
clf;hold on;axis equal
aI.plot('r-','linewidth',2);  
xI.plot('b-','linewidth',2);
plot(real(aIsamp{:}),imag(aIsamp{:}),'ro');
plot(real(xIsamp{:}),imag(xIsamp{:}),'bo');