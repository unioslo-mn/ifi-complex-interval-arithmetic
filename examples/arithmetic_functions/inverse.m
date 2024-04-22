close all

%% Inverse function on a single random rectangle

x1=-1.5;x2=pi/3;x3=-0.2;x4=0.2;
%x1=-pi/3;x2=1.5;x3=-0.2;x4=0.2;
% Random points in the unit square
x1 = rand(1)-1; x2 = 2*rand(1)-1; x3 = 2*rand(1)-1; x4 = 2*rand(1)-1;
x1 = x1*2;
x2 = x2*2;
x3 = x3*2;
x4 = x4*2;
%x1=-1.5;x2=-1;x3=-1;x4=2;
r1 = ciat.RectangularInterval(x1, x2, x3, x4);
tr1 = recip(r1);

% x1b=-4;x2b=3;x3b=3;x4b=4;
% r1b = ciat.RectangularInterval(x1b, x2b, x3b, x4b);

% Define the number of points in each dimension
numPoints = 100;

% Generate grid of complex points in the rectangular region
x = linspace(x1, x2, numPoints);
y = linspace(x3, x4, numPoints);
[X, Y] = meshgrid(x, y);
points = X(:) + 1i*Y(:);

% Generate a few lines of point for visualization
% Horizontal lines
line1 = linspace(x1, x2, numPoints) + 1j*(x4+x3)/2;
line2 = line1 - 1j*(x4-x3)/2;
line3 = line1 + 1j*(x4-x3)/2;

% Vertical lines
line4 = 1i*linspace(x3, x4, numPoints) + (x2+x1)/2;
line5 = line4 - 1*(x2-x1)/2;
line6 = line4 + 1*(x2-x1)/2;


% Their images
tline1 = 1./(line1);
tline2 = 1./(line2);
tline3 = 1./(line3);
tline4 = 1./(line4);
tline5 = 1./(line5);
tline6 = 1./(line6);

% xb = linspace(x1b, x2b, numPoints);
% yb = linspace(x3b, x4b, numPoints);
% [X, Y] = meshgrid(xb, yb);
% pointsb = X(:) + 1i*Y(:);

sqrtuaredPoints = 1./(points);

% Initialize figure
figure;clf

% Plot the original points
subplot(1, 2, 1);
scatter(real(points), imag(points), 5, 'b', 'filled');
hold on;

% Plot lines
scatter(real(line1), imag(line1), 15, 'r', 'filled');
scatter(real(line2), imag(line2), 15, 'c', 'filled');
scatter(real(line3), imag(line3), 15, 'g', 'filled');
scatter(real(line4), imag(line4), 15, 'm', 'filled');
scatter(real(line5), imag(line5), 15, 'y', 'filled');
scatter(real(line6), imag(line6), 15, 'k', 'filled');


% scatter(real(pointsb), imag(pointsb), 5, 'r', 'filled');
r1.plot;
% r1b.plot;
title('Rectangular interval points');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot the sqrtuared points
subplot(1, 2, 2);
scatter(real(sqrtuaredPoints), imag(sqrtuaredPoints), 5, 'b', 'filled');
hold on;

% Plot images of lines
scatter(real(tline1), imag(tline1), 15, 'r', 'filled');
scatter(real(tline2), imag(tline2), 15, 'c', 'filled');
scatter(real(tline3), imag(tline3), 15, 'g', 'filled');
scatter(real(tline4), imag(tline4), 15, 'm', 'filled');
scatter(real(tline5), imag(tline5), 15, 'y', 'filled');
scatter(real(tline6), imag(tline6), 15, 'k', 'filled');

% Plot two vertical lines at tr1.Real extremes
realinf = tr1.Real.Infimum;
realsup = tr1.Real.Supremum;
imaginf = tr1.Imag.Infimum;
imagsup = tr1.Imag.Supremum;
plot([realinf realinf], [min(imag(1./(points))) max(imag(1./(points)))], 'k', 'LineWidth', 5);
plot([realsup realsup], [min(imag(1./(points))) max(imag(1./(points)))], 'k', 'LineWidth', 5);
plot([min(real(1./(points))) max(real(1./(points)))], [imaginf imaginf], 'k', 'LineWidth', 5);
plot([min(real(1./(points))) max(real(1./(points)))], [imagsup imagsup], 'k', 'LineWidth', 5);

title('Mapping of the inverse function');
xlabel('Real^2');
ylabel('Imaginary^2');
axis equal;

%% Inverse function on four predefined rectangular intervals

% Four intervals for a 2x2 matrix
x1=-1;x2=-3;x3=-2;x4=-1;
r1 = ciat.RectangularInterval(x1, x2, x3, x4);

x1b=-4;x2b=-3;x3b=3;x4b=4;
r1b = ciat.RectangularInterval(x1b, x2b, x3b, x4b);

x1c=2;x2c=3;x3c=-2;x4c=-1;
r1c = ciat.RectangularInterval(x1c, x2c, x3c, x4c);

x1d=1;x2d=2;x3d=-2;x4d=-1;
r1d = ciat.RectangularInterval(x1d, x2d, x3d, x4d);

% Define the number of points in each dimension
numPoints = 100;

% Generate grid of complex points in the rectangular region
x = linspace(x1, x2, numPoints);
y = linspace(x3, x4, numPoints);
[X, Y] = meshgrid(x, y);
points = X(:) + 1i*Y(:);

xb = linspace(x1b, x2b, numPoints);
yb = linspace(x3b, x4b, numPoints);
[X, Y] = meshgrid(xb, yb);
pointsb = X(:) + 1i*Y(:);

xc = linspace(x1c, x2c, numPoints);
yc = linspace(x3c, x4c, numPoints);
[X, Y] = meshgrid(xc, yc);
pointsc = X(:) + 1i*Y(:);

xd = linspace(x1d, x2d, numPoints);
yd = linspace(x3d, x4d, numPoints);
[X, Y] = meshgrid(xd, yd);
pointsd = X(:) + 1i*Y(:);


% Compute the inverse of the matrix (here 2x2 matrix)
det = (points.*pointsd - pointsb.*pointsc);
inv =   pointsd  ./ det;
invb = -pointsb ./ det;
invc = -pointsc ./ det;
invd =  points  ./ det;

% Compute example for a single point for verification
idx = 1;
A = [points(idx) pointsb(idx); pointsc(idx) pointsd(idx)];
invA = A^(-1);
Aint = [r1 r1b; r1c r1d];

% Initialize figure
figure;clf

% Plot the original points
subplot(1, 2, 1);
hold on;
scatter(real(points), imag(points), 5, 'b', 'filled');
scatter(real(pointsb), imag(pointsb), 5, 'r', 'filled');
scatter(real(pointsc), imag(pointsc), 5, 'g', 'filled');
scatter(real(pointsd), imag(pointsd), 5, 'y', 'filled');

% Plot example of a single point
scatter(real(A), imag(A), 25, 'k', 'filled');

r1.plot("LineWidth", 3, "Color", "b");
r1b.plot("LineWidth", 3, "Color", "r");
r1c.plot("LineWidth", 3, "Color", "g");
r1d.plot("LineWidth", 3, "Color", "y");

title('Rectangular interval points');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot the coefficients of the inverse matrix
subplot(1, 2, 2);
hold on;
scatter(real(inv), imag(inv),   5, 'b', 'filled');
scatter(real(invb), imag(invb), 5, 'r', 'filled');
scatter(real(invc), imag(invc), 5, 'g', 'filled');
scatter(real(invd), imag(invd), 5, 'y', 'filled');

% Plot example of a single point
scatter(real(invA), imag(invA), 25, 'k', 'filled');

title('Mapping of the inverse function');
xlabel('Real');
ylabel('Imaginary');
axis equal;
