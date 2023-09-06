clear all;
close all

% Sum and product of two rectangular intervals

figure()

r1 =  ciat.RectangularInterval(1, 1.5, 1, 1.5);
r1 = r1.setProbaGrid("normal");

r2 = ciat.RectangularInterval(2, 2.5, 2, 2.5);
r2 = r2.setProbaGrid("uniform");

rsum = r1 + r2;
rprod = r1 .* r2;

r1.plot;
r2.plot;
rsum.plot
rprod.plot
axis equal


% Sum and product of two polygonal intervals

figure()

r1 =  ciat.RectangularInterval(1, 1.5, 1, 1.5);
r1 = r1.setProbaGrid("normal", 'nx', 250, 'ny', 250);
p1 = ciat.PolygonalInterval(r1);

r2 = ciat.PolarInterval(2, 2.5, 2, 2.5);
r2 = r2.setProbaGrid("polarnormal", 'nx', 250, 'ny', 250);
p2 = ciat.PolygonalInterval(r2);

pprod = p1 .* p2;

p1.plot;
p2.plot;
pprod.plot
axis equal
