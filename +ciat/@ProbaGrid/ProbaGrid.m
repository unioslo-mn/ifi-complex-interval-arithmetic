classdef ProbaGrid

    properties
        Pdf
        x
        y
        nx
        ny
        dx
        dy
    end


    methods
        function obj = ProbaGrid(interval, distribution_name, varargin)
            % Takes as an input an interval, a distribution name and
            % optional parameters for the distribution
            % This creates an object representing a probability density
            % function on the interval.

            if nargin == 0
                return
            end

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
                default_sigma = ([width(interval.Abs) / 2, 0; 0, width(interval.Angle) / 2] / 3).^2;
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

            dx = x(2) - x(1);
            dy = y(2) - y(1);

            obj.dx = dx;
            obj.dy = dy;

            obj.nx = nx;
            obj.ny = ny;

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

        % Normalization
        function obj = normalize(obj)
            % Normalize the pdf
            obj.Pdf = obj.Pdf./(sum(obj.Pdf(:)) .* obj.dx .* obj.dy);
        end

        % Addition
        function obj = plus(obj, other)
            % Compute the distribution of the sum of the two random
            % variables. This is done by convolving the two pdfs.

            new_nx = max(obj.nx, other.nx);
            new_ny = max(obj.ny, other.ny);

            % Interpolate the grid with the smallest step size
            [pdf1, pdf2, new_dx, new_dy] = match_step_size(obj, other);

            % Convolve the two pdfs
            pdf = conv2(pdf1, pdf2, "full") .* new_dx .* new_dy;

            % Create a new ProbaGrid object
            new_x = linspace(obj.x(1) + other.x(1), obj.x(end) + other.x(end), length(pdf(1,:)));
            new_y = linspace(obj.y(1) + other.y(1), obj.y(end) + other.y(end), length(pdf(:,1)));
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
            % obj = obj.normalize();  

            % Keep the same grid size
            new_x = linspace(obj.x(1), obj.x(end), new_nx);
            new_y = linspace(obj.y(1), obj.y(end), new_ny);
            obj = obj.adjust(new_x, new_y);
        end

        % Log
        function obj = log(obj)
            % Compute the distribution of the log of the random variable.

            % First compute the new grid
            log_points = log(obj.x + 1i * obj.y.');
            new_x = linspace(min(real(log_points(:))), max(real(log_points(:))), obj.nx);
            new_y = linspace(min(imag(log_points(:))), max(imag(log_points(:))), obj.ny);

            exp_points = exp(new_x + 1i * new_y.');

            pdf = interp2(obj.x, obj.y, obj.Pdf, real(exp_points), imag(exp_points), "linear", 0) .* abs(exp_points);

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
            % obj = obj.normalize();
        end

        % Exp
        function obj = exp(obj)
            % % Compute the distribution of the exp of the random variable.

            % First compute the new grid
            exp_points = exp(obj.x + 1i * obj.y.');
            new_x = linspace(min(real(exp_points(:))), max(real(exp_points(:))), obj.nx);
            new_y = linspace(min(imag(exp_points(:))), max(imag(exp_points(:))), obj.ny);
            [new_X, new_Y] = meshgrid(new_x, new_y);

            % Note that this would need to be summed for each multiple of 2*pi for the actual result
            % However, for the case of the product (which is the only one we need), the result only
            % needs three terms, as the other are zero.
            pdf = (interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)),         imag(log(new_X + 1i * new_Y)), "linear", 0) + ...
                   interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)+2*1j*pi), imag(log(new_X + 1i * new_Y)+2*1j*pi), "linear", 0) + ...
                   interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)-2*1j*pi), imag(log(new_X + 1i * new_Y)-2*1j*pi), "linear", 0)) ...
                    ./ abs(new_X + 1i * new_Y);

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
            % obj = obj.normalize();
        end

        % Multiplication
        function obj = times(obj, other)
            % Compute the distribution of the product of the two random
            % variables.

            new_nx = max(obj.nx, other.nx);
            new_ny = max(obj.ny, other.ny);

            obj = exp(log(obj) + log(other));
            % obj = obj.normalize();

            % Keep the same grid size
            new_x = linspace(obj.x(1), obj.x(end), new_nx);
            new_y = linspace(obj.y(1), obj.y(end), new_ny);
            obj = obj.adjust(new_x, new_y);
        end

        function obj = mtimes(obj, other)
            % Compute the distribution of the product of the two random
            % variables.

            obj = exp(log(obj) + log(other));
            % obj = obj.normalize();
        end


        % 
        function obj = adjust(obj, new_x, new_y)
            % Adjust the grid to the new values of x and y
            % The new values of x and y must be a subset of the old ones
            % The pdf is interpolated to the new grid

            % % Check that the new values of x and y are a subset of the old ones
            % if ~all(ismember(new_x, obj.x)) || ~all(ismember(new_y, obj.y))
            %     error("The new values of x and y must be a subset of the old ones")
            % end

            % Interpolate the pdf to the new grid
            [new_X, new_Y] = meshgrid(new_x, new_y);
            pdf = interp2(obj.x, obj.y, obj.Pdf, new_X, new_Y, "linear", 0);

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
        end

        % Fit to interval
        function obj = fitToInterval(obj, interval)
            % Adjust the grid to the new values of x and y

            new_x = linspace(interval.Real.Infimum, interval.Real.Supremum, obj.nx);
            new_y = linspace(interval.Imag.Infimum, interval.Imag.Supremum, obj.ny);

            obj = obj.adjust(new_x, new_y);
        end

        % Plot
        function plot(obj, varargin)
            % Plot heatmap of the pdf on the grid
            imagesc(obj.x, obj.y, obj.Pdf./max(obj.Pdf(:)), varargin{:})
            % contour(obj.x, obj.y, obj.Pdf./max(obj.Pdf(:)), 10)
            set(gca,'YDir','normal')
            colormap turbo
            axis equal
            colorbar
        end
        
        % Plot areas
        function plot_areas(obj, interval, n_areas, varargin)
            % Separate the pdf into n_areas areas, of the same shape
            % as the original interval but smaller and centered at the
            % baricenter of the original interval
            % For each area, plot the probability to be on that zone
            % as a heatmap
            % Thus the sum of the probabilities of all the areas is 1,
            % and the graphs contains n_areas colored zones.

            % Compute the center of the interval
            center = interval.Real.Midpoint + 1i * interval.Imag.Midpoint;

            interval_centered = interval - center;

            [X,Y] = meshgrid(obj.x, obj.y);

            probas = zeros(n_areas, 1);
            probas_grid = zeros(size(obj.Pdf));
            old_mask = zeros(size(obj.Pdf));

            t_values = linspace(0, 1, n_areas+1);
            hold on
            for i = 1:n_areas
                % Compute the new interval
                new_interval = interval_centered * t_values(i+1) + center;

                % Compute the mask of the new interval
                mask = new_interval.ininterval(X + 1i * Y);
                mask = mask & ~old_mask;

                % Save the mask for the next iteration to avoid overlapping
                old_mask = mask | old_mask;

                % Compute the probability of being in the new interval
                proba_area = sum(obj.Pdf(mask), "all") * obj.dx * obj.dy;

                probas(i) = proba_area;

                % Fill the probability grid with the computed probability
                probas_grid(mask) = proba_area;
            end

            % Plot the heatmap
            imagesc(obj.x, obj.y, probas_grid, 'AlphaData', interval.ininterval(X + 1i * Y), varargin{:})
            set(gca,'YDir','normal')
            % colormap turbo
            axis equal
            colorbar

            % Plot all the intervals
            for i = 1:n_areas
                % Compute the new interval
                new_interval = interval_centered * t_values(i+1) + center;

                % Plot the new interval
                % new_interval.plot("Color", [0 0 0], "LineWidth", 2, 'HandleVisibility','off');
                
                % Plot the interval in red if it's the first with cumulative probability greater than 0.5
                if sum(probas(1:i)) > 0.5 && sum(probas(1:i-1)) < 0.5
                    new_interval.plot("Color", [1 0 0], "LineWidth", 4, 'DisplayName', "50% probability");
                end

                % Same in blue for 0.9
                if sum(probas(1:i)) > 0.9 && sum(probas(1:i-1)) < 0.9
                    new_interval.plot("Color", [0 0 1], "LineWidth", 4, 'DisplayName', "90% probability");
                end

                % Same in green for 0.99
                if sum(probas(1:i)) > 0.99 && sum(probas(1:i-1)) < 0.99
                    new_interval.plot("Color", [0 1 0], "LineWidth", 4, 'DisplayName', "99% probability");
                end 
            end

            
        end
 
    end

    
    % Protected methods
    methods (Access = protected)
        function pdf = interp_pdf(obj, dx, dy)
            % Interpolate the pdf on a grid with step sizes dx and dy

            % Create the new grid
            new_x = obj.x(1):dx:obj.x(end);
            new_y = obj.y(1):dy:obj.y(end);
            [X,Y] = meshgrid(new_x,new_y);

            % Interpolate the pdf
            pdf = interp2(obj.x, obj.y, obj.Pdf, X, Y, "linear", 0);
        end
    end

    % Static methods
    methods (Static)
        function obj = from_pdf(pdf, x, y)
            % Create a ProbaGrid object from a pdf and a grid of points
            % pdf must be a matrix of the same size as x and y
            % x and y must be vectors of the same size

            obj = ciat.ProbaGrid;
            obj.Pdf = pdf;
            obj.x = x;
            obj.y = y;
            obj.nx = length(x);
            obj.ny = length(y);
            obj.dx = x(2) - x(1);
            obj.dy = y(2) - y(1);
        end

    end
end

% Local functions
% Function that takes two ProbaGrid objects and returns their two pdfs,
% where the one with the smallest step size is interpolated on the grid of
% the other
function [pdf1, pdf2, new_dx, new_dy] = match_step_size(obj1, obj2)
    new_dx = max(obj1.dx, obj2.dx);
    new_dy = max(obj1.dy, obj2.dy);

    pdf1 = interp_pdf(obj1, new_dx, new_dy);
    pdf2 = interp_pdf(obj2, new_dx, new_dy);
end
