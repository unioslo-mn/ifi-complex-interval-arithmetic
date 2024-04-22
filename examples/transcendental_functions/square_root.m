% Imaging of the rectangular square root function : smallest rectangle enclosing
% the image of a rectangle by the square root function

close all

x1=-1;x2=-3;x3=-1;x4=2;
r1 = ciat.RectangularInterval(x1, x2, x3, x4);
r1_sqrt = sqrt(r1);

x1b=-4;x2b=3;x3b=3;x4b=4;
r1b = ciat.RectangularInterval(x1b, x2b, x3b, x4b);
r1_sqrtb = sqrt(r1b);

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

% sqrtuare the points
sqrtuaredPoints = sqrt(points);
sqrtuaredPointsb = sqrt(pointsb);

% Plot the original points
subplot(1, 2, 1);
scatter(real(points), imag(points), 5, 'b', 'filled');
hold on;
scatter(real(pointsb), imag(pointsb), 5, 'r', 'filled');
r1.plot;
r1b.plot;
title('Rectangular interval points');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot the sqrtuared points
subplot(1, 2, 2);
scatter(real(sqrtuaredPoints), imag(sqrtuaredPoints), 5, 'b', 'filled');
hold on;
scatter(real(sqrtuaredPointsb), imag(sqrtuaredPointsb), 5, 'r', 'filled');
r1_sqrt.plot;
r1_sqrtb.plot;
title('Mapping of the square-root function');
xlabel('Real^2');
ylabel('Imaginary^2');
axis equal;
