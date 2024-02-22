classdef Arc

	properties
        Center          % Arc center as a complex number
        Radius          % Arc radius as a real number (negative means concave arc, zero means vertex)
	end	

	properties (Dependent)
    	Angles			% Arc angle as a real interval (always counter-clockwise)
		Endpoints		% Endpoints of the arc as complex numbers
        Midpoint        % Midpoint of the arc as complex numbers
		Length			% Length of the arc curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 centered arc (zero is zero centered)
		CurveParameter	% Parameter interval of the normalized curve as real interval
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
    end

    properties (Access = private)
       ArcAngles        % Storage property for the arc angle
    end

	methods
		%% Constructor
        function obj = Arc(varargin)
            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 3
                    center = varargin{1};
                    radius = varargin{2};
                    angles = varargin{3};
                    mustBeA(center,'double')
                    mustBeA(radius,'double')
                    mustBeA(angles,{'double','ciat.RealInterval'})
                    obj.Center = center;
                    obj.Radius = radius;
                    if strcmp(class(angles),'double')
                        obj.Angles = ciat.RealInterval(angles);
                    else
                        obj.Angles = angles;
                    end
            end
		end

        %% Defining properties
                   
        % Set angles (store in the hidden property ArcAngles after wrapping)
        function obj = set.Angles(obj,angles)
            if angles.Width >= 2*pi
                obj.ArcAngles = ciat.RealInterval(-pi,pi);
            else
                obj.ArcAngles = angles - ...
                                2*pi*floor((angles.Infimum+pi)/(2*pi));
            end
        end

        % Get points (retrieve Angles from hidden property ArcAngles)
        function value = get.Angles(obj)
            value = obj.ArcAngles; 
        end

		%% Dependent properties
        
        % Endpoints
		function value = get.Endpoints(obj)
            value = obj.Center + obj.Radius * exp(1i*obj.Angles.Bounds);
            if obj.Radius < 0
                value = flip(value);
            end
        end

        % Midpoint
        function value = get.Midpoint(obj)
            value = obj.Center + obj.Radius * exp(1i*obj.Angles.Midpoint);
        end

        % Length
        function value = get.Length(obj)
            value = 2*pi*obj.Radius * obj.Angles.Width/(2*pi);
        end

        % Gauss map angle interval
        function value = get.GaussMap(obj)
            value = obj.Angles;
            if ~isempty(value) && value.isin(pi) && all(pi~=value.Bounds)
                value = [ciat.RealInterval(wrapToPi(value.Infimum),pi);...
                         ciat.RealInterval(-pi,wrapToPi(value.Supremum))];
            end
        end

        % Log-Gauss map angle interval
        function value = get.LogGaussMap(obj)
            value = obj.GaussMap - ciat.RealInterval(angle(obj.Endpoints));
        end

        % Normalization factor
        function value = get.NormFactor(obj)
            if obj.Center == 0
                value = 0;
            else
                value = 1/obj.Center;
            end
        end

        % Curve parameter
        function value = get.CurveParameter(obj)
            value = obj.Angles + angle(obj.NormFactor);
        end

        % Real
        function value = get.Real(obj)
            if obj.Radius ~= 0
                conc = obj.Radius < 0;
                if sum( obj.Angles.isin(([-1,1,3]+conc)*pi) )
                    realInf = real(obj.Center) - abs(obj.Radius);
                else
                    realInf = min(real(obj.Endpoints));
                end
                if sum(obj.Angles.isin(([0,2]+conc)*pi))
                    realSup = real(obj.Center) + abs(obj.Radius);
                else
                    realSup = max(real(obj.Endpoints));
                end
                value = ciat.RealInterval(realInf,realSup);
            else
                value = ciat.RealInterval(real(obj.Center));
            end
        end
        function value = real(obj)
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
        end
        
        % Imag
        function value = get.Imag(obj)
            if obj.Radius ~= 0
                conc = obj.Radius < 0;
                if sum(obj.Angles.isin(([-0.5,1.5]+conc)*pi))
                    imagInf = imag(obj.Center) - abs(obj.Radius);
                else
                    imagInf = min(imag(obj.Endpoints));
                end
                if sum(discretize(([0.5,2.5]+conc)*pi,obj.Angles.Bounds)==1)
                    imagSup = imag(obj.Center) + abs(obj.Radius);
                else
                    imagSup = max(imag(obj.Endpoints));
                end
                value = ciat.RealInterval(imagInf,imagSup);
            else
                value = ciat.RealInterval(imag(obj.Center));
            end
        end
        function value = imag(obj)
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            if obj.Radius ~= 0
                conc = obj.Radius < 0;
                if sum(obj.Angles.isin(angle(obj.Center)+([-1,1,3]+conc)*pi))
                    absInf = abs(obj.Center) - abs(obj.Radius);
                else
                    absInf = min(abs(obj.Endpoints));
                end
                if sum(obj.Angles.isin(angle(obj.Center)+([0,2]+conc)*pi))
                    absSup = abs(obj.Center) + abs(obj.Radius);
                else
                    absSup = max(abs(obj.Endpoints));
                end
                value = ciat.RealInterval(absInf,absSup);
            else
                value = ciat.RealInterval(abs(obj.Center));
            end
        end
        function value = abs(obj)
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            if obj.Radius ~= 0
                conc = obj.Radius < 0;
                if sum(obj.Angles.isin(angle(obj.Center)+([-0.5,1.5]+conc)*pi))
                    angleInf = angle(obj.Center) - ...
                        asin(abs(obj.Radius/obj.Center));
                else
                    angleInf = min(angle(obj.Endpoints));
                end
                if sum(obj.Angles.isin(angle(obj.Center)+([0.5,2.5]+conc)*pi))
                    angleSup = angle(obj.Center) + ...
                        asin(abs(obj.Radius/obj.Center));
                else
                    angleSup = max(angle(obj.Endpoints));
                end
                value = ciat.RealInterval(angleInf,angleSup);      
            else
                value = ciat.RealInterval(angle(obj.Center));
            end
        end
        function value = angle(obj)
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
        end

		%% Other methods

        % Plot
        function h = plot(obj, varargin)
        % Plot arcs
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
                % Extract parameters
                cr = real(obj(n).Center);
                ci = imag(obj(n).Center);
                r = obj(n).Radius;
                % Plot arc by angle quadrants
                q = [-pi,-pi/2 ; -pi/2,0 ; 0,pi/2 ; pi/2,pi];
                if r~=0
                    for idx = 1:4
                        qi = intersection([obj(n).Angles, ...
                                            ciat.RealInterval(q(idx,:))]);
                        if ~isempty(qi)
                            xBound = sort(cr + r*cos(qi.Bounds) );
                            yBound = sort(ci + r*sin(qi.Bounds) );
                            h = [h; ...
                                fimplicit(@(x,y) (x-cr).^2 + (y-ci).^2 - r^2 , ...
                                            [xBound yBound],varargin{:})];    
                        end
                    end
                else
                    h = [h; plot(cr,ci,varargin{:},'Marker','.')];
                end
            end
            if tf == false 
                hold off
            end
        end

        % Plot Gauss maps
        function h = plotMap(obj, logMap, arrowSize, varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];

            % Plot normal vectors
            for n = 1:length(obj(:))
                % Extract variables
                arc = obj(n);
                if logMap == 0
                    map = arc.GaussMap;
                else
                    map = arc.LogGaussMap;
                end

                % Set arrow positions
                p = zeros(1,3);
                p(1) = arc.Endpoints(1);
                p(2) = arc.Midpoint;
                p(3) = arc.Endpoints(2);

                % Set vector lengths
                if length(map) == 1
                    a = zeros(1,3);
                else
                    a = zeros(1,6);
                    p = reshape(repmat(p,1,2),1,[]);
                end

                % Set arrow angle and length
                for m = 1:length(map)
                    a(1+3*(m-1)) = arrowSize * exp(1i*map(m).Infimum);
                    a(2+3*(m-1)) = arrowSize * exp(1i*map(m).Midpoint);
                    a(3+3*(m-1)) = arrowSize * exp(1i*map(m).Supremum);
                end

                % Flip vectors if arc is concave 
                if obj(n).Radius < 0
                    a = flip(a);
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