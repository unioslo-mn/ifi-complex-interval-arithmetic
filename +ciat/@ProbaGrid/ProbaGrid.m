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

            % Parse optional parameters using inputParser
            p = inputParser;
            addOptional(p, 'nx', default_nx, @isscalar);
            addOptional(p, 'ny', default_ny, @isscalar);
            addParameter(p, 'mu', [interval.Real.Midpoint, interval.Imag.Midpoint], @(x) isnumeric(x) && numel(x) == 2);
            addParameter(p, 'Sigma', [width(interval.Real) / 2, 0; 0, width(interval.Imag) / 2] / 3, @ismatrix);
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
                    % Restrict the pdf to the interval
                case "normal"


                    obj.Pdf = mvnpdf([X(:),Y(:)], mu, sigma);
                    obj.Pdf = reshape(obj.Pdf, size(X));

                    % Restrict the pdf to the interval

                    % obj.Pdf = obj.Pdf./sum(obj.Pdf(:).*dx.*dy);
                otherwise
                    error("Unknown pdf name")
            end


        end


        
        % Addition
        function obj = plus(obj, other)
            % Compute the distribution of the sum of the two random
            % variables. This is done by convolving the two pdfs.

            % Interpolate the grid with the smallest step size
            [pdf1, pdf2, new_dx, new_dy] = match_step_size(obj, other);

            % Convolve the two pdfs
            pdf = conv2(pdf1, pdf2, "full") .* new_dx .* new_dy;

            % Create a new ProbaGrid object
            % new_x = obj.x(1) + other.x(1):new_dx:obj.x(end) + other.x(end);
            % new_y = obj.y(1) + other.y(1):new_dy:obj.y(end) + other.y(end);
            new_x = linspace(obj.x(1) + other.x(1), obj.x(end) + other.x(end), length(pdf(1,:)));
            new_y = linspace(obj.y(1) + other.y(1), obj.y(end) + other.y(end), length(pdf(:,1)));
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
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

            % LITTLE CHEAT
            % Normalize the pdf
            new_dx = new_x(2) - new_x(1);
            new_dy = new_y(2) - new_y(1);
            pdf = pdf ./ sum(pdf(:) .* new_dx .* new_dy);

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
        end

        % Exp
        function obj = exp(obj)
            % % Compute the distribution of the exp of the random variable.

            % First compute the new grid
            exp_points = exp(obj.x + 1i * obj.y.');
            new_x = linspace(min(real(exp_points(:))), max(real(exp_points(:))), obj.nx);
            new_y = linspace(min(imag(exp_points(:))), max(imag(exp_points(:))), obj.ny);
            [new_X, new_Y] = meshgrid(new_x, new_y);

            % log_points = log(new_x + 1i * new_y.');

            % pdf =( interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points), "linear", 0) + ...
            %     interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points) + 2*pi, "linear", 0) + ...
            %     interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points) - 2*pi, "linear", 0) + ...
            %     interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points) + 4*pi, "linear", 0) + ...
            %     interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points) - 4*pi, "linear", 0) + ...
            %     interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points) + 6*pi, "linear", 0) + ...
            %     interp2(obj.x, obj.y, obj.Pdf, real(log_points), imag(log_points) - 6*pi, "linear", 0) ...
            %     ) ./ abs(exp_points);




            % ---------------------------------------------------


            % $$f_{exp(log X + log Y)}(z) = f_{log X + log Y}(ln(z)) / z$$

            % x = new_x;
            % y = new_y;
            % [X,Y] = meshgrid(x,y);

            % pdfsumloghandle = @(x, y) interp2(obj.x, obj.y, obj.Pdf, x, y, "linear", 0);

            % pdfprodhandle = @(x, y) (pdfsumloghandle(real(log(x + 1i * y)), imag(log(x + 1i * y))) + ...
            %                          pdfsumloghandle(real(log(x + 1i * y)+2*1j*pi), imag(log(x + 1i * y)+2*1j*pi)) + ...
            %                          pdfsumloghandle(real(log(x + 1i * y)-2*1j*pi), imag(log(x + 1i * y)-2*1j*pi))) ...
            %                          ./ abs(x + 1i * y);
            % pdfprod = pdfprodhandle(X, Y);

            % Unfold the computation of the pdf
            pdf = (interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)),         imag(log(new_X + 1i * new_Y)), "linear", 0) + ...
                       interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)+2*1j*pi), imag(log(new_X + 1i * new_Y)+2*1j*pi), "linear", 0) + ...
                       interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)-2*1j*pi), imag(log(new_X + 1i * new_Y)-2*1j*pi), "linear", 0)) ...
                    ./ abs(new_X + 1i * new_Y);


            % LITTLE CHEAT
            % Normalize the pdf
            new_dx = new_x(2) - new_x(1);
            new_dy = new_y(2) - new_y(1);
            pdf = pdf ./ sum(pdf(:) .* new_dx .* new_dy);

            % Create a new ProbaGrid object
            obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);

        end

        % Multiplication
        function obj = mtimes(obj, other)
            % Compute the distribution of the product of the two random
            % variables.

            obj = exp(log(obj) + log(other));
        end

        % Plot
        function plot(obj, varargin)
            % Plot heatmap of the pdf on the grid
            imagesc(obj.x, obj.y, obj.Pdf, varargin{:})
            set(gca,'YDir','normal')
            colormap turbo
            axis equal
            colorbar
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
