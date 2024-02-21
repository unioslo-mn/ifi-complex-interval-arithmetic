classdef Arc

	properties
        Center          % Arc center as a complex number
        Radius          % Arc radius as a real number (negative means concave arc, zero means vertex)
    	Angle			% Arc angle as a real interval (always counter-clockwise)
	end	

	properties (Dependent)
		Endpoints		% Endpoints of the arc as complex numbers
		Length			% Length of the arc curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 centered arc (zero is zero centered)
		CurveParameter	% Parameter interval of the normalized curve as real interval
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
                    angle = varargin{3};
                    mustBeA(center,'double')
                    mustBeA(radius,'double')
                    mustBeA(angle,'ciat.RealInterval')
                    obj.Center = center;
                    obj.Radius = radius;
                    obj.Angle = ciat.RealInterval(angle.Bounds);
            end
		end

		%% Defining properties
        
        % Endpoints
		function value = get.Endpoints(obj)
            value = obj.Center + obj.Radius * exp(1i*obj.Angle.Bounds);
        end

        % Length
        function value = get.Length(obj)
            value = 2*pi*obj.Radius * obj.Angle.Width/(2*pi);
        end

        % Gauss map angle interval
        function value = get.GaussMap(obj)
            value = obj.Angle;
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
            value = obj.Angle + angle(obj.NormFactor);
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
                % Extract parameters
                cr = real(obj(n).Center);
                ci = imag(obj(n).Center);
                r = obj(n).Radius;
                % Plot arc by angle quadrants
                q = [-pi,-pi/2 ; -pi/2,0 ; 0,pi/2 ; pi/2,pi];
                for idx = 1:4
                    qi = intersection([obj(n).Angle, ...
                                        ciat.RealInterval(q(idx,:))]);
                    if ~isempty(qi)
                        xBound = sort(cr + r*cos(qi.Bounds) );
                        yBound = sort(ci + r*sin(qi.Bounds) );
                        h = [h; ...
                            fimplicit(@(x,y) (x-cr).^2 + (y-ci).^2 - r^2 , ...
                                        [xBound yBound],varargin{:})];    
                    end
                end
            end
            if tf == false 
                hold off
            end
        end

	end


end