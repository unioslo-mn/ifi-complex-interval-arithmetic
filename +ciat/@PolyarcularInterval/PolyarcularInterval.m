classdef PolyarcularInterval

	properties
        Tolerance;      % Maximum distance from the boundary of the represented interval
    end

    properties (Dependent)
        Arcs;        	% Center points of the arcs defining the polygonal interval boundary
        ArcCount;     	% Number of vertex points
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
        Area;           % Area of the polygonal interval
    end

    properties (Access = private)
       Boundary         % Storage property for the boundary arc segments
    end

	methods
		%% Constructor
        function obj = PolyarcularInterval(inObj,inObj2,optional)

        end
	end


	 %% Static methods
     methods (Static)
        % Function headers
        outObj = segmentInverse(obj)
        outObj = segmentProduct(obj1, obj2)
    end
end