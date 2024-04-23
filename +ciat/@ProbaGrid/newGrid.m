function obj = newGrid(obj,interval, distribution_name, varargin)

    % Takes as an input an interval, a distribution name and
    % optional parameters for the distribution
    % This creates an object representing a probability density
    % function on the interval.

    % Default values for nx and ny
    default_nx = 100;
    default_ny = 100;

    % Default values for the mean and the covariance matrix
    if distribution_name == "polarnormal"
        % If the distribution is a polar normal distribution, the
        % default parameters are the mean and the covariance matrix
        % of the distribution in polar coordinates
        default_mu = [interval.Abs.Midpoint, interval.Angle.Midpoint];
        % sigma_r = width(interval.Abs) / 2 / 3;
        % sigma_theta = width(interval.Angle) / 2 / 3;
        % default_sigma = [sigma_r, 0; 0, sigma_theta].^2;
        default_sigma = ([width(interval.Abs) / 2, 0; 0, ...
                          width(interval.Angle) / 2] / 3).^2;
    else
        default_mu = [interval.Real.Midpoint, interval.Imag.Midpoint];
        default_sigma = ([width(interval.Real) / 2, 0; 0, width(interval.Imag) / 2] / 3).^2;
    end

    % Parse optional parameters using inputParser
    p = inputParser;
    addOptional(p, 'nx', default_nx, @isscalar);
    addOptional(p, 'ny', default_ny, @isscalar);
    addParameter(p, 'mu', default_mu, @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'Sigma', default_sigma, @ismatrix);
    parse(p, varargin{:});

    % Extract the parameters
    nx = p.Results.nx;
    ny = p.Results.ny;
    mu = p.Results.mu;
    sigma = p.Results.Sigma;

    % Create a grid of points according to the parameters
    x = linspace(interval.Real.Infimum, interval.Real.Supremum, nx);
    y = linspace(interval.Imag.Infimum, interval.Imag.Supremum, ny);
    [X,Y] = meshgrid(x,y);

    obj.x = x;
    obj.y = y;

    % Check if the grid is real (1D) or complex (2D)
    if ny == 1
        % The grid is real
        obj.nx = nx;
        obj.ny = 1;
        dx = x(2) - x(1);
        dy = 0;
    else
        % The grid is complex
        obj.nx = nx;
        obj.ny = ny;
        dx = x(2) - x(1);
        dy = y(2) - y(1);
    end

    obj.dx = dx;
    obj.dy = dy;

    % obj.nx = nx;
    % obj.ny = ny;

    switch distribution_name
        case "uniform"
            obj.Pdf = ones(length(x),length(y))./interval.Area;

        case "normal"
            obj.Pdf = mvnpdf([X(:),Y(:)], mu, sigma);
            obj.Pdf = reshape(obj.Pdf, size(X));

        case "polarnormal"
            % Gaussian distribution along radial and angular axis
            % The parameters are the mean and the standard deviation
            % of the radial and angular distributions

            R = sqrt(X.^2 + Y.^2);
            Theta = atan2(Y, X);

            % In this case, the parameters must be given in polar coordinates
            mu_r = mu(1);
            mu_theta = mu(2);
            sigma_r = sqrt(sigma(1, 1)); % Independence is assumed
            sigma_theta = sqrt(sigma(2, 2)); % Independence is assumed

            % This method does not work because the angle distribution
            % is not considered as a circular distribution
            % Aswell, the radial distribution should only have
            % positive support
            % Instead, we use the von Mises distribution as an approximation
            % of the normal distribution on the circle, and the folded
            % normal distribution for the radial distribution


            % If sigma_theta is too small, the von Mises distribution
            % is not a good approximation of the normal distribution
            % on the circle. In this case, we use the normal
            % distribution on the circle instead.

            if sigma_theta < 5e-2
                pdf_theta = normpdf(Theta, mu_theta, sigma_theta);
            else
                % Compute the pdf
                % We use the von Mises distribution to compute the pdf, as an approximation of the wrapped normal distribution
                pdf_theta = exp(cos(Theta - mu_theta) / sigma_theta^2) / (2 * pi * besseli(0, 1 / sigma_theta^2));
            end

            % For the radial distribution, we use a folded normal distribution
            pdf_r = normpdf(R, mu_r, sigma_r) + normpdf(R, -mu_r, sigma_r);

            % Compute the pdf
            obj.Pdf = pdf_r .* pdf_theta ./ R;

        otherwise
            error("Unknown pdf name")
    end

    % Restrict the pdf to the interval and normalize it
    obj.Pdf = obj.Pdf .* interval.ininterval(X + 1i * Y);
    %obj.Pdf = obj.Pdf./sum(obj.Pdf(:)).* dx.*dy;
    % obj = obj.normalize();

    % Check if the pdf contains a NaN value
    if any(isnan(obj.Pdf(:)))
        error("The pdf contains a NaN value")
    end
end