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
            
            if nargin ~= 0
                obj = newGrid(obj,interval, distribution_name, varargin{:});
            end

        end

        % Normalization
        function obj = normalize(obj)
            % Normalize the pdf
            if obj.dy > 0
                obj.Pdf = obj.Pdf ./ (sum(obj.Pdf(:)) .* obj.dx .* obj.dy);
            else
                % obj.Pdf = obj.Pdf ./ (sum(obj.Pdf) .* obj.dx);
                obj.Pdf = obj.Pdf ./ sum(obj.Pdf);
            end
        end

        function obj = mtimes(obj, other)
            % Compute the distribution of the product of the two random
            % variables.

            obj = exp(log(obj) + log(other));
            obj = obj.normalize();
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
            obj = obj.normalize();
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
            obj = obj.normalize();
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
        obj = newGrid(obj,interval, distribution_name, varargin)
        plotAreas(obj, interval, n_areas, varargin)
        
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

