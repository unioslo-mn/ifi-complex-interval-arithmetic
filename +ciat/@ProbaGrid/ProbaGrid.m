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

        % Normalization
        function obj = normalize(obj)
            % Normalize the pdf
            obj.Pdf = obj.Pdf./(sum(obj.Pdf(:)) .* obj.dx .* obj.dy);
        end

        function obj = mtimes(obj, other)
            % Compute the distribution of the product of the two random
            % variables.

            obj = exp(log(obj) + log(other));
            % obj = obj.normalize();
        end

        % Real part
        function obj = real(obj)
            % Compute the distribution of the real part of the random
            % variable.

            % First compute the new grid
            new_x = linspace(min(obj.x(:)), max(obj.x(:)), obj.nx);
            new_y = 0;

            pdf = sum(obj.Pdf, 1) .* obj.dy;

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
            % obj = obj.normalize();
        end

        % Imaginary part
        function obj = imag(obj)
            % Compute the distribution of the imaginary part of the random
            % variable.

            % First compute the new grid
            new_x = linspace(min(obj.y(:)), max(obj.y(:)), obj.ny);
            new_y = 0;

            pdf = sum(obj.Pdf, 2) .* obj.dx;

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
            % obj = obj.normalize();
        end

        % Adjust grid
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
            [M,N] = size(obj);
            
            for m = 1:M
                for n = 1:N
                    new_x = linspace(interval(m,n).Real.Infimum, ...
                                     interval(m,n).Real.Supremum, obj(m,n).nx);
                    new_y = linspace(interval(m,n).Imag.Infimum, ...
                                     interval(m,n).Imag.Supremum, obj(m,n).ny);
                    obj(m,n) = obj(m,n).adjust(new_x, new_y);
                end
            end
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

        %% Function headers
        obj = plus(obj, other)
        r = sum(obj,varargin)
        obj = log(obj)
        obj = exp(obj)
        obj = times(obj, other)
        obj = abs2(obj)
        plot_areas(obj, interval, n_areas, varargin)
        
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

            % Check if the grid is complex
            if obj.ny == 1
                % obj.is_complex = false;
                obj.dx = x(2) - x(1);
                obj.dy = 0;
            else
                % obj.is_complex = true;
                obj.dx = x(2) - x(1);
                obj.dy = y(2) - y(1);
            end
        end

    end
end

