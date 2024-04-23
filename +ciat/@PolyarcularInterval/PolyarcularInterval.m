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
                    % Fix the angles property of zero radius arcs
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
                        arcs(n).Angles = ciat.RealInterval(angMin,angMax);
                    end
                end
            else
                arcs.Angles = ciat.RealInterval(-pi,pi);
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

        % Get implicit vertices
        function value = get.Vertices(obj)
            % Filter out vertices with zero radius
            N = obj.ArcCount;
            M = sum([obj.Arcs.Radius]~=0);
            if N > 1 && M > 0
                value(2*M,1) = ciat.Arc;
                
                % Initialize vertex-pair index
                m = 1;
                for n=1:N
                    arc = obj.Arcs(n);
                    if arc.Radius ~= 0
                        % Extract neighbouring edges
                        edgeNext = obj.Edges(n);
                        if n==1
                            edgePrev = obj.Edges(N);
                        else
                            edgePrev = obj.Edges(n-1);                   
                        end
        
                        % Calculate angle intervals
                        ang1Min = edgePrev.GaussMap.Midpoint;
                        ang2Max = edgeNext.GaussMap.Midpoint;
                        if arc.Radius >= 0
                            ang1Max = arc.Angles.Infimum;
                            ang2Min = arc.Angles.Supremum;
                        else
                            ang1Max = arc.Angles.Supremum;
                            ang2Min = arc.Angles.Infimum;
                        end

                        % Adjust order of angle bounds
                        if ang1Min > ang1Max
                            ang1Max = ang1Max + 2*pi;
                        end
                        if ang2Min > ang2Max
                            ang2Max = ang2Max + 2*pi;
                        end
                        
                        % Set vertex parameters
                        value(2*m-1).Center = arc.Endpoints(1);
                        value(2*m-1).Radius = 0;
                        value(2*m-1).Angles = ciat.RealInterval(ang1Min,ang1Max);
        
                        value(2*m).Center = arc.Endpoints(2);
                        value(2*m).Radius = 0;
                        value(2*m).Angles = ciat.RealInterval(ang2Min,ang2Max);
                        
                        % Increase vertex-pair index
                        m = m + 1;
                    end
                end
            else
                value = [];
            end
        end

        % Real
        function value = get.Real(obj)
            value = ciat.RealInterval(min(inf(real(obj.Arcs))),...
                                      max(sup(real(obj.Arcs))));
        end
        function value = real(obj)
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
        end
        
        % Imag
        function value = get.Imag(obj)
            value = ciat.RealInterval(min(inf(imag(obj.Arcs))),...
                                      max(sup(imag(obj.Arcs))));
        end
        function value = imag(obj)
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            value = ciat.RealInterval(min(min(inf(abs(obj.Arcs))),...
                                            min(inf(abs(obj.Edges))) ),...
                                      max(sup(abs(obj.Arcs))));
        end
        function value = abs(obj)
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            value = ciat.RealInterval(min(inf(angle(obj.Arcs))),...
                                      max(sup(angle(obj.Arcs))));         
        end
        function value = angle(obj)
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
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
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                h = [h;obj(n).Arcs.plot(varargin{:})];    
                h = [h;obj(n).Edges.plot(varargin{:})];    
            end
            if tf == false 
                hold off
            end
        end

        % Plot Gauss maps
        function h = plotGaussMap(obj, arrowSize, varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];

            % Plot normal vectors
            % for n = 1:obj.ArcCount
                h = [h; obj.Arcs.plotGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Edges.plotGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Vertices.plotGaussMap(arrowSize,varargin{:})];
            % end

            if tf == false 
                hold off
            end
        end
        function h = plotLogGaussMap(obj, arrowSize, varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];

            % Plot normal vectors
            for n = 1:obj.ArcCount
                h = [h; obj.Arcs(n).plotLogGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Edges(n).plotLogGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Vertices(n).plotLogGaussMap(arrowSize,varargin{:})];
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
        outObj = cast(inObj,options)
    end
end