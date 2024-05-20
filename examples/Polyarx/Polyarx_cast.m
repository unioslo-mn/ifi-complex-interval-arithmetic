clear
% close all

%%

% Rectangle
rI = ciat.RectangularInterval(rand,rand,rand,rand);
rxI = ciat.PolyarxInterval(rI);

% Circle
cI = ciat.CircularInterval(complex(rand,rand),rand);
cxI = ciat.PolyarxInterval(cI);

% Polar
pI = ciat.PolarInterval(rand,rand,rand*2*pi,rand*2*pi);
pxI = ciat.PolyarxInterval(pI);


%% Plot
% figure;clf
subplot(2,2,1);cla;hold on;axis equal
rI.plot('c-')
rxI.plot('k--')
subplot(2,2,2);cla;hold on;axis equal
cI.plot('c-')
cxI.plot('k--')
subplot(2,2,3);cla;hold on;axis equal
pI.plot('c-')
pxI.plot('k--')