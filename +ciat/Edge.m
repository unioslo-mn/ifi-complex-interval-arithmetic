classdef Edge

	properties
        Endpoints       % Endpoints of the edge as complex numbers
    end	

	properties (Dependent)
        Vector          % Endpoint difference as a complex vector
        Slope           % Slope of the line the edge fits on
        ZeroCrossing    % Real axis crossing point of the line the edge fits on
		Length			% Length of the edge curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 crossing vertical edge (zero is zero crossing)
		CurveParameter	% Parameter interval of the normalized curve as real interval
	end

	methods
		%% Constructor
        function obj = Edge(varargin)
            for varIdx = 1:length(varargin)
                mustBeA(varargin{varIdx},'double')
            end

            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 1
                    obj.Endpoints = [varargin{1},varargin{1}];
                case 2
                    obj.Endpoints = [varargin{1},varargin{2}];
            end
		end

		%% Defining properties

        % Vector
        function value = get.Vector(obj)
            value = diff(obj.Endpoints);
        end

        % Slope
		function value = get.Slope(obj)
            vd = obj.Vector;
            value = imag(vd)/real(vd);
        end

        % Zero crossing
		function value = get.ZeroCrossing(obj)
            v1 = obj.Endpoints(1);
            value = real(v1) - imag(v1) / obj.Slope;
        end

        % Length
		function value = get.Length(obj)
            value = abs(obj.Vector);
        end

        % GaussMap
		function value = get.GaussMap(obj)
            value = ciat.RealInterval(ciat.wrapPi(angle(obj.Vector)-pi/2));
        end

        % Log-GaussMap
		function value = get.LogGaussMap(obj)
            value = obj.GaussMap - ciat.RealInterval(angle(obj.Endpoints));
        end

        % Normalization factor
		function value = get.NormFactor(obj)
            if obj.ZeroCrossing == 0
                value = 0;
            else
                % Extract line parameters
                v1 = obj.Endpoints(1);
                va = obj.Slope;        
                
                % Find scale-rotate factor
                vr = exp(-1i*(atan(va)+pi/2));   % rotation factor
                vs = 1/real(vr * v1);       % scale factor
                value = vr * vs;                        % complex factor

            end
        end

        %% Other methods

        % Plot
        function h = plot(obj, varargin)
        % Plot polygonal intervals 
        %
        % This function plots a set of polygonal intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.Arc class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   
        % _________________________________________________________________________
            tf = ishold;
            hold on
            h = [];
            for n = 1:length(obj(:))
                points = obj(n).Endpoints;
                h = [h;plot(real(points), imag(points), varargin{:})];    
            end
            if tf == false 
                hold off
            end
        end

    end


end