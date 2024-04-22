close all

%% Square function on two random rectangular intervals
% Imaging of the rectangular square function : smallest rectangle enclosing
% the image of a rectangle by the square function
% Compared to the multiplication of the interval by itself (with dependency)

x1=-1;x2=-3;x3=-1;x4=2;
r1 = ciat.RectangularInterval(x1, x2, x3, x4);
r1_sq = r1^2;
r1_t = r1*r1;

x1b=-4;x2b=5;x3b=3;x4b=4;
r1b = ciat.RectangularInterval(x1b, x2b, x3b, x4b);
r1_sqb = r1b^2;
r1_tb = r1b*r1b;

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

% Square the points
squaredPoints = points.^2;
squaredPointsb = pointsb.^2;

% Plot the original points
subplot(1, 2, 1);
scatter(real(points), imag(points), 5, 'b', 'filled','DisplayName','A');
hold on;
scatter(real(pointsb), imag(pointsb), 5, 'r', 'filled','DisplayName','B');
r1.plot('DisplayName','\partialA');
r1b.plot('DisplayName','\partialB');
title('Rectangular interval points');
xlabel('Real');
ylabel('Imaginary');
axis equal;
legend

% Plot the squared points
subplot(1, 2, 2);
scatter(real(squaredPoints), imag(squaredPoints), 5, 'b', 'filled','DisplayName','A^2');
hold on;
scatter(real(squaredPointsb), imag(squaredPointsb), 5, 'r', 'filled','DisplayName','B^2');
r1_sq.plot('DisplayName','(\partialA)^2');
r1_sqb.plot('DisplayName','(\partialB)^2');
r1_t.plot('DisplayName','(\partialA)\times\partialA)');
r1_tb.plot('DisplayName','(\partialB)\times\partialB)');
title('Mapping of the square function');
xlabel('Real^2');
ylabel('Imaginary^2');
axis equal;
legend
