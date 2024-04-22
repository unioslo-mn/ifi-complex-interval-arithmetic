close all

%% Logarithm function on a random rectangular interval

% Random points in the unit square
x1 = rand(1)-1; x2 = 2*rand(1)-1; x3 = 2*rand(1)-1; x4 = 2*rand(1)-1;
x1 = x1*2;
x2 = x2*2;
x3 = x3*2;
x4 = x4*2;

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
tline1 = log(line1);
tline2 = log(line2);
tline3 = log(line3);
tline4 = log(line4);
tline5 = log(line5);
tline6 = log(line6);

% sqrtuare the points
sqrtuaredPoints = log(points);

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

title('Mapping of the logarithm function');
xlabel('Real');
ylabel('Imaginary');
axis equal;
