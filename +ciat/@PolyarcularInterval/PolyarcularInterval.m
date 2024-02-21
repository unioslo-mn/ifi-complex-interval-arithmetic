classdef PolyarcularInterval

	properties
        Tolerance;      % Maximum distance from the boundary of the represented interval
    end

    properties (Dependent)
        Arcs;        	% Arcs defining the polygonal interval boundary
        ArcCount;     	% Number of vertex points
        Edges;          % Implicit edges
        Vertices;       % Implicit vertices 
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
        function obj = PolyarcularInterval(inObj)
            arguments
                inObj                (:,:)   = []
            end    

            switch class(inObj)
                case 'double'
                    if isempty(inObj)
                        % This is for initializing an array of objects
                    else
                        obj = ciat.PolyarcularInterval.cast(inObj);
                    end
                case 'ciat.Arc'
                    % This is the default way of defining polygonal
                    % intervals the points are assumed to belong to a
                    % single interval no matter how many dimensions
                    obj.Arcs = inObj(:);
                case 'cell'
                    % This is how multiple polyarcs can be defined using
                    % cells of ciat.Arc arrays
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    for n = 1:M*N
                        obj(n).Arcs = inObj{n};
                    end
                otherwise
                    % Input object will be casted
                    tol = obj.Tolerance;
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    for n = 1:M*N
                        obj(n) = ciat.PolyarcularInterval.cast(inObj(n));
                    end
            end % switch
        end % function
        
        %% Defining properties
                   
        % Set points (store in the hidden property Boundary after sorting)
        function obj = set.Arcs(obj,arcs)
            N = length(arcs);
            if N > 1
                for n=1:N
                    if arcs(n).Radius == 0
                        % set current, previous and next arc
                        arcCurr = arcs(n);
                        if n > 1
                            arcPrev = arcs(n-1);
                        else
                            arcPrev = arcs(N);
                        end
                        if n < N
                            arcNext = arcs(n+1);
                        else
                            arcNext = arcs(1);
                        end
                        
                        angMin = angle(arcCurr.Center - ...
                                       arcPrev.Endpoints(2)) - pi/2;
                        angMax = angle(arcNext.Endpoints(1) - ...
                                       arcCurr.Center) - pi/2;
                        if angMin > angMax
                            angMax = angMax + 2*pi;
                        end
                        arcs(n).Angle = ciat.RealInterval(angMin,angMax);
                    end
                end
            else
                arcs.Angle = ciat.RealInterval(-pi,pi);
            end
            obj.Boundary = arcs;
        end 
    
        % Get points (retrieve from hidden property Boundary)
        function value = get.Arcs(obj)
            value = obj.Boundary;
        end

        %% Dependent properties
                        
        % Get arc count
        function value = get.ArcCount(obj)
            value = length(obj.Arcs);
        end

        % Get implicit edges
        function value = get.Edges(obj)
            N = obj.ArcCount;
            value(N,1) = ciat.Edge;
            for n=1:N-1
                value(n).Endpoints = [obj.Arcs(n).Endpoints(2);...
                                      obj.Arcs(n+1).Endpoints(1)];
            end
            value(N).Endpoints = [obj.Arcs(N).Endpoints(2);...
                                      obj.Arcs(1).Endpoints(1)];
        end

        %% Other functions

        % Plot
        function h = plot(obj, varargin)
        % Plot polygonal intervals 
        %
        % This function plots a set of polyarcular intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolyarcularInterval class
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
                for m = 1:length(obj(n).Arcs)
                    h = [h;obj(n).Arcs(m).plot];    
                end
            end
            if tf == false 
                hold off
            end
        end

    end % methods

	 %% Static methods
     methods (Static)
        % Function headers
        outObj = segmentInverse(obj)
        outObj = segmentProduct(obj1, obj2)
    end
end