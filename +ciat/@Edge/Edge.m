classdef Edge

	properties
        Endpoints       % Endpoints of the edge as complex numbers
    end	

	properties (Dependent)
        Midpoint        % Midpoint of the arc as complex numbers
        Vector          % Endpoint difference as a complex vector
        Slope           % Slope of the line the edge fits on
        ZeroCrossing    % Real axis crossing point of the line the edge fits on
		Length			% Length of the edge curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 crossing vertical edge (zero is zero crossing)
		CurveParameter	% Parameter interval of the normalized curve as real interval
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
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

		%% Dependent properties

        % Midpoint
        function value = get.Midpoint(obj)
            value = sum(obj.Endpoints)/2;
        end

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
            % value = ciat.RealInterval(ciat.wrapToPi(angle(obj.Vector)-pi/2));
            value = ciat.wrapToPi(angle(obj.Vector)-pi/2);
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

        function value = get.CurveParameter(obj)
            value = ciat.RealInterval(imag(obj.Endpoints * obj.NormFactor));
        end

        % Real
        function value = get.Real(obj)
            value = ciat.RealInterval(min(real(obj.Endpoints)),...
                                      max(real(obj.Endpoints)));
        end
        function value = real(obj)
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
        end
        
        % Imag
        function value = get.Imag(obj)
            value = ciat.RealInterval(min(imag(obj.Endpoints)),...
                                      max(imag(obj.Endpoints)));
        end
        function value = imag(obj)
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            if obj.CurveParameter.isin(0)
                minAbs = 1/abs(obj.NormFactor);
            else
                minAbs = min(abs(obj.Endpoints));
            end
            maxAbs = max(abs(obj.Endpoints));
            value = ciat.RealInterval(minAbs,maxAbs);
        end
        function value = abs(obj)
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            value = ciat.RealInterval(min(angle(obj.Endpoints)),...
                                      max(angle(obj.Endpoints)));
        end
        function value = angle(obj)
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
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
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                points = obj(n).Endpoints;
                if obj(n).Length ~= 0
                    h = [h;plot(real(points), imag(points), varargin{:})];    
                else
                    h = [h;plot(real(points), imag(points), varargin{:},...
                        'Marker','.')];    
                end
            end
            if tf == false 
                hold off
            end
        end

        % Plot
        function h = plotMap(obj,logMap,arrowSize,varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];

            % Plot normal vectors
            for n = 1:length(obj(:))
                % Extact variables
                edge = obj(n);
                if logMap == 0
                    map = edge.GaussMap;
                else
                    map = edge.LogGaussMap;
                end

                % Set arrow positions
                p(1) = edge.Endpoints(1);
                p(2) = edge.Midpoint;
                p(3) = edge.Endpoints(2);

                % Set vector lengths
                if length(map) == 1
                    a = zeros(1,3);
                else
                    a = zeros(1,6);
                    p = repmat(p,1,2);
                end

                % Set arrow angle and length
                for m = 1:length(map)
                    a(1+3*(m-1)) = arrowSize * exp(1i*map(m).Infimum);
                    a(2+3*(m-1)) = arrowSize * exp(1i*map(m).Midpoint);
                    a(3+3*(m-1)) = arrowSize * exp(1i*map(m).Supremum);
                end
                
                % Plot arrows
                h=[h; ...
                    quiver(real(p),imag(p),real(a),imag(a),...
                           varargin{:},'AutoScale','off')];
            end

            if tf == false 
                hold off
            end
        end
        function h = plotGaussMap(obj, arrowSize, varargin)
            h = obj.plotMap(0,arrowSize,varargin{:});
        end
        function h = plotLogGaussMap(obj, arrowSize, varargin)
            h = obj.plotMap(1,arrowSize,varargin{:});
        end

    end


end