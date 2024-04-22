% Test to see if it is better for the interval to compute the function as x²/y² or (x/y)²

x1=1;x2=2;x3=1;x4=2;
r1 = ciat.RectangularInterval(x1, x2, x3, x4);

x1b=2;x2b=3;x3b=-1;x4b=0;
r1b = ciat.RectangularInterval(x1b, x2b, x3b, x4b);

% Corners of the squares
c1 = x1 + 1j * x3; % Bottom left
c2 = x1 + 1j * x4; % Top left
c3 = x2 + 1j * x3; % Bottom right
c4 = x2 + 1j * x4; % Top right

% Corners of the squares
c1b = x1b + 1j * x3b; % Bottom left
c2b = x1b + 1j * x4b; % Top left
c3b = x2b + 1j * x3b; % Bottom right
c4b = x2b + 1j * x4b; % Top right

res1 = (r1./r1b).^2;
res2 = r1.^2./r1b.^2;

% Define the number of points in each dimension
numPoints = 20;

% Generate grid of complex points in the rectangular region
points = linspace(x1,x2,numPoints) + 1j * linspace(x3,x4,numPoints).';
pointsb = linspace(x1b,x2b,numPoints) + 1j * linspace(x3b,x4b,numPoints).';

respoints = points(:).^2 ./ (pointsb(:).^2).';



%% Plot
close all

% Plot the original points
subplot(1, 2, 1);
hold on;
scatter(real(points), imag(points), 5, 'b', 'filled', 'HandleVisibility', 'off');
scatter(real(pointsb), imag(pointsb), 5, 'r', 'filled', 'HandleVisibility', 'off');

r1.plot('k');
r1b.plot('k');


% Plot r1./r1b
divres = r1./r1b;
divres.plot('b');
divrespoints = points(:) ./ pointsb(:).';
scatter(real(divrespoints), imag(divrespoints), 2, 'b', 'filled', 'HandleVisibility', 'off');



title('Original Points');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot the sqrtuared points
subplot(1, 2, 2);
% close all
% figure()
hold on;
scatter(real(respoints), imag(respoints), 2, 'k', 'filled', 'HandleVisibility', 'off');

res1.plot('k', 'LineWidth', 2, 'DisplayName', '(r1/r1b)²');
res2.plot('k', 'LineWidth', 2, 'DisplayName', 'r1²/r1b²');
resinter = intersection([res1, res2]);
resinter.plot('r', 'LineWidth', 2, 'DisplayName', 'Intersection');

xlabel('Real');
ylabel('Imaginary');
axis equal;
legend('show');