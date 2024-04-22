% Compute and plot the inverse of a matrix of rectangular intervals
% Each coefficient of the matrix is a rectangular interval and is plotted separately

% % Four intervals for a 2x2 matrix
% x1=-1;x2=-3;x3=-2;x4=-1;
% r1 = ciat.RectangularInterval(x1, x2, x3, x4);

% x1b=-4;x2b=-3;x3b=3;x4b=4;
% r1b = ciat.RectangularInterval(x1b, x2b, x3b, x4b);

% x1c=2;x2c=3;x3c=-2;x4c=-1;
% r1c = ciat.RectangularInterval(x1c, x2c, x3c, x4c);

% x1d=1;x2d=2;x3d=-2;x4d=-1;
% r1d = ciat.RectangularInterval(x1d, x2d, x3d, x4d);

% % Define the number of points in each dimension
% numPoints = 100;

% % For each interval, sample points in the rectangle

% % Generate grid of complex points in the rectangular region
% x = linspace(x1, x2, numPoints);
% y = linspace(x3, x4, numPoints);
% [X, Y] = meshgrid(x, y);
% points = X(:) + 1i*Y(:);

% xb = linspace(x1b, x2b, numPoints);
% yb = linspace(x3b, x4b, numPoints);
% [X, Y] = meshgrid(xb, yb);
% pointsb = X(:) + 1i*Y(:);

% xc = linspace(x1c, x2c, numPoints);
% yc = linspace(x3c, x4c, numPoints);
% [X, Y] = meshgrid(xc, yc);
% pointsc = X(:) + 1i*Y(:);

% xd = linspace(x1d, x2d, numPoints);
% yd = linspace(x3d, x4d, numPoints);
% [X, Y] = meshgrid(xd, yd);
% pointsd = X(:) + 1i*Y(:);


% % Compute the inverse of the matrix (here 2x2 matrix)
% det = (points.*pointsd - pointsb.*pointsc);
% inv =   pointsd  ./ det;
% invb = -pointsb ./ det;
% invc = -pointsc ./ det;
% invd =  points  ./ det;

% % Compute example for a single point for verification
% idx = 1;
% A = [points(idx) pointsb(idx); pointsc(idx) pointsd(idx)];
% invA = A^(-1);












% Now deal with real matrices

epsilon = 0.0001;
r1 = ciat.RealInterval(1-epsilon,1+epsilon);
r1b = ciat.RealInterval(-epsilon,epsilon);
r1c = ciat.RealInterval(-epsilon,epsilon);
r1d = ciat.RealInterval(1-epsilon,1+epsilon);

Aint = [r1 r1b; r1c r1d];
Ac = Aint.Midpoint;
Delta = Aint.Width/2;
I = eye(2);
M = (I - abs(Ac^-1)*Delta)^(-1);
mu = diag(M);
Tnu = (2*diag(mu) - I)^-1;
Btdown = -M*abs(Ac^-1) + Tnu*(Ac^-1 + abs(Ac^-1));
Btup   =  M*abs(Ac^-1) + Tnu*(Ac^-1 - abs(Ac^-1));







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

title('Original Points');
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

% rinvA1.plot("LineWidth", 3, "Color", "b");
% rinvA2.plot("LineWidth", 3, "Color", "r");
% rinvA3.plot("LineWidth", 3, "Color", "g");
% rinvA4.plot("LineWidth", 3, "Color", "y");


title('Inverse Matrix Coefficients');
xlabel('Real');
ylabel('Imaginary');
axis equal;
