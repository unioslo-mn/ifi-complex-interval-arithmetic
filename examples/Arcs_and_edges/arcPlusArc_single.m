clear
% close all

%%  Generate a random arc

% Create arc array
centers = 3*complex(rand(2,1),rand(2,1));
radii = 2*rand(2,1)-1;
angInf = 2*pi*rand(2,1);
angSup = 2*pi*rand(2,1);
arcs = ciat.Arc(centers, radii, ciat.RealInterval(angInf,angSup));
  

% Add arcs
arcPlusArc = arcs(1) + arcs(2);

%% 

% Normalize arcs
arcsNorm = arcs + (-arcs.Center);
arcPlusArcNorm = arcPlusArc + (-arcPlusArc.Center);

% figure;clf
% Plot arcs
subplot(1,2,1);cla;hold on;axis equal;title('Original')
arcs.plot('k','lineWidth',2);
arcPlusArc.plot('b','lineWidth',2);
% Plot normalized arcs
subplot(1,2,2);cla;hold on;axis equal;title('Normalized')
arcsNorm.plot('k','lineWidth',2)
arcsNorm.plotMap(0,.1,'k');
arcPlusArcNorm.plot('b','lineWidth',2)
arcPlusArcNorm.plotMap(0,.1,'b');
fimplicit(@(x,y) x.^2 + y.^2 - radii(1).^2,'k:')
fimplicit(@(x,y) x.^2 + y.^2 - radii(2).^2,'k:')
fimplicit(@(x,y) x.^2 + y.^2 - sum(radii).^2,'b:')


