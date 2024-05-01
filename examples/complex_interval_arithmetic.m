%% Single complex intervals (defined by double type)

clear;
close all

figure(1);hold on;axis equal
plot(0,0,'k+')
xlabel('Real')
ylabel('Imag')
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultContourLineWidth',2)
set(0,'DefaultFunctionContourLineWidth',2)
set(0,'defaultimplicitfunctionlineLineWidth',2)
title('Example results of complex interval arithmetic operations')

% Rectangular interval
rI = ciat.RectangularInterval(1,2,2,4);
rI.plot('b-');

% Circular interval
cI = ciat.CircularInterval(-2+3i,1);
cI.plot('b-');

% Polar interval
pI = ciat.PolarInterval(1,2,-2,-3);
pI.plot('b-');

% Polygonal interval
gI = ciat.PolygonalInterval([1-1i, 3-2i, 2-3i, 1-2i]);
gI.plot('b-');

% Polyarcular interval
% aI = ciat.PolyarcularInterval([ciat.Arc(3+0.5i,1,ciat.RealInterval([-.7,-.4]*pi)),...
%                                ciat.Arc(4+2.5i,-0.7,ciat.RealInterval([.3,.5]*pi)),...
%                                ciat.Arc(2+1.5i,0,ciat.RealInterval(0))]);
% aI.plot('b-');
%% Cast complex interval types

% Rectangular to circular
rcI = ciat.CircularInterval(rI);
rcI.plot('r--');

% Rectangular to polar
rpI = ciat.PolarInterval(rI);
rpI.plot('r--');

% Rectangular to polygonal
rgI = ciat.PolygonalInterval(rI);
rgI.plot('k:','LineWidth',2);

% % Rectangular to polyarcular
% raI = ciat.PolyarcularInterval(rI);
% raI.plot('k:','LineWidth',2);

% Circular to rectangular
crI = ciat.RectangularInterval(cI);
crI.plot('r--');

% Circular to polar
cpI = ciat.PolarInterval(cI);
cpI.plot('r--');

% Circular to polygonal
cgI = ciat.PolygonalInterval(cI);
cgI.plot('r--');

% % Circular to polyarcular
% caI = ciat.PolyarcularInterval(cI);
% caI.plot('k:','LineWidth',2);

% Polar to rectangular
prI = ciat.RectangularInterval(pI);
prI.plot('r--');

% Polar to circular
pcI = ciat.CircularInterval(pI);
pcI.plot('r--');

% Polar to polygonal
pgI = ciat.PolygonalInterval(pI);
pgI.plot('r--');

% % Polar to polyarcular
% paI = ciat.PolyarcularInterval(pI);
% paI.plot('k:','LineWidth',2);

% Polygonal to rectangular
grI = ciat.RectangularInterval(gI);
grI.plot('r--');

% Polygonal to circular
gcI = ciat.CircularInterval(gI);
gcI.plot('r--');

% Polygonal to polar
gpI = ciat.PolarInterval(gI);
gpI.plot('r--');

% % Polygonal to polyarcular
% gaI = ciat.PolyarcularInterval(gI);
% gaI.plot('k:','LineWidth',2);

% % Polyarcular to rectangular
% arI = ciat.RectangularInterval(aI);
% arI.plot('r--');

% % Polyarcular to polar
% apI = ciat.PolarInterval(aI);
% apI.plot('r--');

%% Transform complex intervals

% Negative intervals of all types
plot(-rI,'c-');
plot(-cI,'c-');
plot(-pI,'c-');
plot(-gI,'c-');

%% Combine complex intervals

% Add rectangular intervals
rI_plus = rI + crI;
rI_plus.plot('g-');

% Multiply rectangular intervals
rI_times = rI * crI;
rI_times.plot('g-');

% Add circular intervals
cI_plus = cI + rcI;
cI_plus.plot('g-');

% Multiply rectangular intervals
cI_times = cI * rcI;
cI_times.plot('g-');

% Multiply polar intervals
pI_times = cpI * rpI;
pI_times.plot('g-');

% Add polygonal intervals
gI_plus = rgI + cgI;
gI_plus.plot('g-');

% Add polygonal intervals
gI_plus = rgI + cgI;
gI_plus.plot('g-')

% Multiply polygonal intervals
gI_times = rgI * cgI;
gI_times.plot('g-');

% % Add polyarcular intervals
% aI_plus = ciat.PolyarcularInterval(-pI) + aI;
% aI_plus.plot('g-');

%% Unite or intersect complex intervals

% Union of rectangular intervals
cI_union = union([-rI,prI]);
cI_union.plot('r--');

% Intersection of rectangular intervals
gI_union = intersection([ciat.RectangularInterval(-cI),grI]);
gI_union.plot('k--');

% Union of circular intervals
cI_union = union([-cI,gcI]);
cI_union.plot('r--');

% Union of polar intervals
cI_union = union([cpI,pI_times]);
cI_union.plot('r--');

% Union of polygonal intervals
gI_union = union([ciat.PolygonalInterval(pI_times),gI_times]);
gI_union.plot('r--');

% Intersection of polygonal intervals
gI_union = intersection([ciat.PolygonalInterval(pI_times),gI_times]);
gI_union.plot('k--');