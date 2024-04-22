close all

%% Tangen function on a random polar interval
% Imaging of the rectangular tan function : smallest rectangle enclosing
% the image of a rectangle by the tan function

% Define a polar sector
% r1=1; r2=2;
% theta1=pi/4; theta2=pi/2;
% Random polar sector
r1 = 1 + 2*rand;
r2 = r1 + 2*rand;
theta1 = pi*rand;
theta2 = theta1 + pi*rand;
p1 = ciat.PolarInterval(r1, r2, theta1, theta2);

% Define the number of points in each dimension
numPoints = 100;

% Generate grid of complex points in the polar region
r = linspace(r1, r2, numPoints);
theta = linspace(theta1, theta2, numPoints);
points = exp(1j*theta).'*r;

tanPoints = tan(points);

% Initialize figure
figure;clf

% Plot the original points
subplot(1, 2, 1);
scatter(real(points), imag(points), 5, 'filled');
hold on;

title('Polar interval points');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot the sqrtuared points
subplot(1, 2, 2);
scatter(real(tanPoints), imag(tanPoints), 5, 'filled');
hold on;

title('Mapping of the tangent function');
xlabel('Real');
ylabel('Imaginary');
axis equal;

%% Tangen function on a random rectangular interval

% Imaging of the rectangular tan function : smallest rectangle enclosing
% the image of a rectangle by the tan function

% Random points in the unit square
x1 = rand(1)-1; x2 = 2*rand(1)-1; x3 = 2*rand(1)-1; x4 = 2*rand(1)-1;
x1 = x1*2;
x2 = x2*2;
x3 = x3*2;
x4 = x4*2;


r1 = ciat.RectangularInterval(x1, x2, x3, x4);
tr1 = tan(r1);

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
tline1 = tan(line1);
tline2 = tan(line2);
tline3 = tan(line3);
tline4 = tan(line4);
tline5 = tan(line5);
tline6 = tan(line6);

tanPoints = tan(points);

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


r1.plot;
% r1b.plot;
title('Rectangular interval points');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot the sqrtuared points
subplot(1, 2, 2);
scatter(real(tanPoints), imag(tanPoints), 5, 'b', 'filled');
hold on;

% Plot images of lines
scatter(real(tline1), imag(tline1), 15, 'r', 'filled');
scatter(real(tline2), imag(tline2), 15, 'c', 'filled');
scatter(real(tline3), imag(tline3), 15, 'g', 'filled');
scatter(real(tline4), imag(tline4), 15, 'm', 'filled');
scatter(real(tline5), imag(tline5), 15, 'y', 'filled');
scatter(real(tline6), imag(tline6), 15, 'k', 'filled');

% Plot two vertical lines at tr1.Real extremes
tr1.plot('k', 'LineWidth', 1);

% scatter(real(tanPointsb), imag(tanPointsb), 5, 'r', 'filled');
title('Mapping of the tangent function');
xlabel('Real');
ylabel('Imaginary');
axis equal;
