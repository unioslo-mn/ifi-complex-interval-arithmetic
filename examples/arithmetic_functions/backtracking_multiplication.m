%% Backtracking test
% Resulting interval of the multiplication of two squres, take one point one the border
% of the product interval, and backtrack the original factors in the squares

close all

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

rprod = r1 * r1b;

% Define the number of points in each dimension
numPoints = 20;

% Generate grid of complex points in the rectangular region
points = linspace(x1,x2,numPoints) + 1j * linspace(x3,x4,numPoints)';
pointsb = linspace(x1b,x2b,numPoints) + 1j * linspace(x3b,x4b,numPoints)';

prod = points(:) * pointsb(:).';


% Find the vertex of the product of maximum modulus
alt = zeros(4,1);
alt(1) = rprod.Real.Infimum + 1j * rprod.Imag.Infimum;
alt(2) = rprod.Real.Infimum + 1j * rprod.Imag.Supremum;
alt(3) = rprod.Real.Supremum + 1j * rprod.Imag.Infimum;
alt(4) = rprod.Real.Supremum + 1j * rprod.Imag.Supremum;
[~,ind] = max(abs(alt));
z = alt(ind);
[~,ind] = max(abs(prod(:)));
z = prod(ind);
z = (c4)*c3b;

absv = ciat.RealInterval(abs(z)/r1b.abs.Supremum, abs(z)/r1b.abs.Infimum);
angles = ciat.RealInterval(angle(z) - r1b.angle.Supremum, angle(z) - r1b.angle.Infimum);
r1_inv = ciat.PolarInterval(absv, angles);
r1_inv2 = ciat.RectangularInterval(r1_inv);
tmp = r1_inv2 * r1b;

E_inv_sample = r1_inv.sample(100);
center = (x2-x1)/2 + 1j * (x4-x3)/2;

% Compute distances from E_inv_sample to square of vertices c1, c2, c3, c4
% Loop over E_inv_sample and compute distance to the square
dist = zeros(length(E_inv_sample),1);
for i=1:length(E_inv_sample)
    % Right side
    if real(E_inv_sample(i)) > x2
        if imag(E_inv_sample(i)) > x4
            % Above right corner
            dist(i) = abs(E_inv_sample(i) - c4);
        elseif imag(E_inv_sample(i)) < x3
            % Below right corner
            dist(i) = abs(E_inv_sample(i) - c3);
        else
            % Right side
            dist(i) =abs(x2 - real(E_inv_sample(i)));
        end
    % Left side
    elseif real(E_inv_sample(i)) < x1
        if imag(E_inv_sample(i)) > x4
            % Above left corner
            dist(i) = abs(E_inv_sample(i) - c2);
        elseif imag(E_inv_sample(i)) < x3
            % Below left corner
            dist(i) = abs(E_inv_sample(i) - c1);
        else
            % Left side
            dist(i) = abs(real(E_inv_sample(i)) - x1);
        end
    else
        % Top side
        if imag(E_inv_sample(i)) > x4
            % Top side
            dist(i) = abs(imag(E_inv_sample(i)) - x4);
        % Bottom side
        elseif imag(E_inv_sample(i)) < x3
            % Bottom side
            dist(i) = abs(x3 - imag(E_inv_sample(i)));
        else
            % Inside
            dist(i) = 0;
        end
    end
end

% [~, min_idx] = min( abs( E_inv_sample - center ) );
[~, min_idx] = min( dist );
A_marked = E_inv_sample(min_idx);
E_marked = z / A_marked;

% Plot both grid of points and the product of the two
figure;clf;hold on
plot(real(points(:)),imag(points(:)),'b.','DisplayName','A')
plot(real(pointsb(:)),imag(pointsb(:)),'r.','DisplayName','B')
plot(real(prod(:)),imag(prod(:)),'g.','DisplayName','A\timesB')
plot(real(z),imag(z),'k*','DisplayName','Z')
r1.plot('b','DisplayName','\partial A');
r1b.plot('r','DisplayName','\partial B');
rprod.plot('g','DisplayName','\partial (A\timesB)');
plot(real(A_marked),imag(A_marked),'ko','DisplayName','A\capZ/B')
plot(real(E_marked),imag(E_marked),'kx','DisplayName','B\capZ/A')
plot(z./pointsb(:),'k.','DisplayName','Z/B');
r1_inv.plot('k','DisplayName','\partial(Z/B)');
axis equal
grid on
legend()