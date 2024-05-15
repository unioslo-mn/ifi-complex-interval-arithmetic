classdef PolyarcularInterval

	properties (Dependent)
        Arcs;           % Defining arcs of the polygonal interval boundary
        Edges;          % Implicit edges
        Vertices;       % Implicit vertices 
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
        Area;           % Area of the polygonal interval
    end

    properties (Access = private)
       ArcStorage         % Storage property for the boundary arc segments
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
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    for n = 1:M*N
                        obj(n) = ciat.PolyarcularInterval.cast(inObj(n));
                    end
            end % switch
        end % function
        
        %% Defining properties
                   
        % Set points (store in the hidden property ArcStorage after sorting)
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
                                       arcPrev.Endpoint) - pi/2;
                        angMax = angle(arcNext.Startpoint - ...
                                       arcCurr.Center) - pi/2;
                        if angMin > angMax
                            angMax = angMax + 2*pi;
                        end
                        arcs(n).ArcAngle = ciat.RealInterval(angMin,angMax);
                    end
                end
            else
                arcs.ArcAngle = ciat.RealInterval(-pi,pi);
            end
            obj.ArcStorage = arcs;
        end 
    
        % Get points (retrieve from hidden property ArcStorage)
        function value = get.Arcs(obj)
            value = obj.ArcStorage(obj.ArcStorage.Length ~= 0);
        end

        %% Dependent properties

                       
        % Get implicit edges
        function value = get.Edges(obj)
            % N = length(obj.ArcStorage);
            % value(N,1) = ciat.Edge;
            % for n=1:N-1
            %     value(n).Startpoint =  obj.Arcs(n).Endpoint;
            %     value(n).Endpoint = obj.Arcs(n+1).Startpoint;
            % end
            % value(N).Startpoint =  obj.Arcs(N).Endpoint;
            % value(N).Endpoint = obj.Arcs(1).Startpoint;

            value = ciat.Edge(obj.ArcStorage.Endpoint , ...
                              circshift(obj.ArcStorage.Startpoint,-1));
            value = value(value.Length>10*eps);
        end

        % Get implicit vertices
        function value = get.Vertices(obj)

             % Extract elements
            arcStart = obj.Arcs.Startpoint;
            arcEnd = obj.Arcs.Endpoint;
            arcRadius = obj.Arcs.Radius;
            arcGauss = obj.Arcs.GaussMap;
            arcStartAngle = (arcRadius>0) .* arcGauss.Infimum + ...
                            (arcRadius<0) .* arcGauss.Supremum;
            arcEndAngle = (arcRadius>0) .* arcGauss.Supremum + ...
                          (arcRadius<0) .* arcGauss.Infimum;
            edgeStart = obj.Edges.Startpoint;
            edgeEnd = obj.Edges.Endpoint;
            edgeAngle = obj.Edges.GaussMap.Infimum;

            % Define vertices at the intersection of arc pairs
            [M,N] = find(abs(arcEnd - arcStart.')<10*eps);
            angStart = arcEndAngle(M);
            angEnd = arcStartAngle(N);
            angDif = wrapToPi(angEnd - angStart);
            angInf = (angDif>=0) .* angStart + (angDif<0) .* angEnd;
            angSup = angInf + sign(angDif) .* angDif;
            adjArcVertices = ciat.Arc(arcEnd(M),  zeros(length(M),1),...
                                       ciat.RealInterval(angInf,angSup) );


            % Define vertices at the intersection of edge pairs
            [M,N] = find(abs(edgeEnd - edgeStart.')<10*eps);
            angStart = edgeAngle(M);
            angEnd = edgeAngle(N);
            angDif = wrapToPi(angEnd - angStart);
            angInf = (angDif>=0) .* angStart + (angDif<0) .* angEnd;
            angSup = angInf + sign(angDif) .* angDif;
            % angSup = angSup + (angInf>0 & angSup<0)*2*pi;
            adjEdgeVertices = ciat.Arc(edgeEnd(M), zeros(length(M),1),...
                                       ciat.RealInterval(angInf,angSup) );

            % Define vertices at the intersection of arc-edge pairs
            [M,N] = find(abs(arcEnd - edgeStart.')<10*eps);
            angStart = arcEndAngle(M);
            angEnd = edgeAngle(N);
            angDif = wrapToPi(angEnd - angStart);
            angInf = (angDif>=0) .* angStart + (angDif<0) .* angEnd;
            angSup = angInf + sign(angDif) .* angDif;
            % angSup = angSup + (angInf>0 & angSup<0)*2*pi;
            adjArcEdgeVertices = ciat.Arc(arcEnd(M),  zeros(length(M),1),...
                                       ciat.RealInterval(angInf,angSup) );


            % Define vertices at the intersection of edge-arc pairs
            [M,N] = find(abs(edgeEnd - arcStart.')<10*eps);
            angStart = edgeAngle(M);
            angEnd = arcStartAngle(N);
            angDif = wrapToPi(angEnd - angStart);
            angInf = (angDif>=0) .* angStart + (angDif<0) .* angEnd;
            angSup = angInf + sign(angDif) .* angDif;
            % angSup = angSup + (angInf>0 & angSup<0)*2*pi;
            adjEdgeArcVertices = ciat.Arc(edgeEnd(M),  zeros(length(M),1),...
                                       ciat.RealInterval(angInf,angSup) );

            value = [adjArcVertices ; ...
                    adjEdgeVertices ; ...
                    adjArcEdgeVertices ; ...
                    adjEdgeArcVertices];

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

        % Area
        function value = get.Area(obj)
            [M,N] = size(obj);
            value = zeros(M,N);
            % Calculate the polygon area for each polyarcular interval
            for m = 1:M
                for n = 1:N
                    if isempty(obj.ArcStorage)
                        value(m,n) = nan;
                    else
                        points = [obj(m,n).ArcStorage.Startpoint , ...
                                  obj(m,n).ArcStorage.Endpoint].';
                        polygonArea = polyarea(real(points(:)), ...
                                               imag(points(:)));
                        arcArea = sum(obj(m,n).Arcs.Area);
                        value(m,n) = polygonArea + arcArea;
                    end
                end
            end
        end

        %% Other functions

        % Sample
        function value = sample(obj, nPoints)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    % Sample arcs and edges
                    arcs = obj(m,n).ArcStorage;
                    edges = ciat.Edge(arcs.Endpoint , ...
                                      circshift(arcs.Startpoint,-1));
                    arcPoints = sample(arcs,nPoints);
                    edgePoints = sample(edges,nPoints);

                    % Interleave arc and edge samples
                    allPoints = [arcPoints.';edgePoints.'];
                    allPoints = allPoints(:);
                    value{m,n} = [allPoints{:}];
                end
            end
        end

        % % IsNaN
        % function r = isnan(obj)
        %     r = isnan(obj.Area);
        % end 

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

        %% Function headers
        r = sum(obj,varargin)

    end % methods

	 %% Static methods
     methods (Static)
        % Function headers
        outObj = segmentInverse(obj)
        outObj = segmentProduct(obj1, obj2)
        outObj = cast(inObj,options)
        [arcOut,edgeOut] = splitSegments(arcIn,edgeIn)
        [arcOut,edgeOut] = trimSegments(arcIn,edgeIn)
        seg = orderSegments(obj)
     end


  
end