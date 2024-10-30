classdef PolyarcularInterval < matlab.mixin.indexing.RedefinesParen
	properties (Dependent)
        DefArcs         % Defining arcs of the polygonal interval boundary
        Arcs;           % Implicit arcs 
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
        function obj = PolyarcularInterval(inObj,inObj2,optional)
            arguments
                inObj                (:,:)   = []
                inObj2               (:,:)   = []
                optional.convex      (1,1)   {mustBeNumericOrLogical} = false
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
                    obj(1) = ciat.PolyarcularInterval;
                    obj.Arcs = {inObj(:)};
                    
                case 'cell'
                    % This is how multiple polyarcs can be defined using
                    % cells of ciat.Arc arrays
                    [M,N] = size(inObj);
                    obj(M,N) = ciat.PolyarcularInterval;
                    obj.Arcs = inObj;
                otherwise
                    % Input object will be casted
                    if isempty(inObj2)
                        [M,N] = size(inObj);
                        obj(M,N) = obj;
                        for n = 1:M*N
                            obj(n) = ciat.PolyarcularInterval.cast(inObj(n));
                        end
                    else
                        % Give the option to give only one second object or
                        % the same number than the first object
                        [M,N] = size(inObj);
                        [N2] = length(inObj2(:));
                        assert(N2 == M*N || N2==1)
                        obj(M,N) = obj;
                        for n = 1:M*N
                            obj(n) = ciat.PolyarcularInterval.cast(inObj(n),...
                                                     inObj2(min(n,N2)));
                        end
                    end
            end % switch
            if optional.convex
                obj = obj.convexify;
            end
        end % function
        
        %% Defining properties
                   
        % Set points (store in the hidden property ArcStorage after sorting)
        function obj = set.Arcs(obj,arcs)
            [M,N] = size(arcs);
            for m = 1:M
                for n = 1:N
                    obj.ArcStorage(m,n) = {ciat.PolyarcularInterval.setArcs( ...
                                                        arcs{m,n})};
                end
            end
        end 
    
         % Get arcs (retrieve from hidden property ArcStorage)
        function value = get.DefArcs(obj)
            value = obj.ArcStorage;
        end

        %% Dependent properties

        % Get implicit arcs
        function value = get.Arcs(obj)
            [M,N] = size(obj);
            % value = cell(M,N);
            value = obj.ArcStorage;
            for m = 1:M
                for n = 1:N
                    % arcs = obj.ArcStorage{m,n};
                    % value{m,n} = arcs(arcs.Length ~= 0);
                    value{m,n} = value{m,n}(value{m,n}.Length~=0);
                end
            end

            if M*N==1
                value = value{:};
            end
        end               

        % Get implicit edges
        function value = get.Edges(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    arcs = obj.ArcStorage{m,n};
                    edges = ciat.Edge(arcs.Endpoint , ...
                                      circshift(arcs.Startpoint,-1));
                    value{m,n} = edges(edges.Length>10*eps);
                end
            end

            if M*N==1
                value = value{:};
            end
        end

        % Get implicit vertices
        function value = get.Vertices(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            arcs = obj.Arcs;
            edges = obj.Edges;
            for m = 1:M
                for n = 1:N
                    [value{m,n},~] = obj(m,n).getVertices;
                end
            end

            if M*N==1
                value = value{:};
            end
        end

        % Real
        function value = get.Real(obj)
            [M,N] = size(obj);
            minReal = zeros(M,N);
            maxReal = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    minReal(m,n) = min(inf(real(obj.ArcStorage{m,n})));
                    maxReal(m,n) = max(sup(real(obj.ArcStorage{m,n})));
                end
            end
            value = ciat.RealInterval( minReal,maxReal );
        end
        function value = real(obj)
            value = obj.Real;
        end
        
        % Imag
        function value = get.Imag(obj)
            [M,N] = size(obj);
            minImag = zeros(M,N);
            maxImag = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    minImag(m,n) = min(inf(imag(obj.ArcStorage{m,n})));
                    maxImag(m,n) = max(sup(imag(obj.ArcStorage{m,n})));
                end
            end
            value = ciat.RealInterval( minImag,maxImag );
        end
        function value = imag(obj)
            value = obj.Imag;
        end
        
        % Abs
        function value = get.Abs(obj)
            [M,N] = size(obj);
            minAbs = zeros(M,N);
            maxAbs = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    if ~isempty(obj.Edges{m,n})
                        minAbs(m,n) = min(min(inf(abs(obj.ArcStorage{m,n}))),...
                                          min(inf(abs(obj.Edges{m,n}))) );
                    else
                        minAbs(m,n) = min(inf(abs(obj.ArcStorage{m,n})));
                    end
                    maxAbs(m,n) = max(sup(abs(obj.ArcStorage{m,n})));
                    
                end
            end
            value = ciat.RealInterval( minAbs,maxAbs );

            % Check for intervals that contain the zero
            pointIn = obj.isin(0);
            if any(pointIn,'all')
                value(pointIn).Infimum = 0;
            end
        end
        function value = abs(obj)
            value = obj.Abs;
        end
        
        % Angle
        function value = get.Angle(obj)
            [M,N] = size(obj);
            minAng = zeros(M,N);
            maxAng = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    arcs = obj.ArcStorage{m,n};
                    if obj(m,n).isin(0)
                        minAng(m,n) = -pi;
                        maxAng(m,n) = pi;
                    elseif obj(m,n).Imag.isin(0) && ...
                           obj(m,n).Real.Supremum < 0
                        minAng(m,n) = min(wrapToPi(inf(angle(arcs)+pi)))+pi;
                        maxAng(m,n) = max(wrapToPi(sup(angle(arcs)+pi)))+pi;
                    else
                        minAng(m,n) = min(inf(angle(arcs)));
                        maxAng(m,n) = max(sup(angle(arcs)));
                    end
                end
            end
            value = ciat.RealInterval( minAng,maxAng );
        end
        function value = angle(obj)
            value = obj.Angle;
        end

        % Area
        function value = get.Area(obj)
            [M,N] = size(obj);
            value = zeros(M,N);
            % Calculate the polygon area for each polyarcular interval
            for m = 1:M
                for n = 1:N
                    if isempty(obj.ArcStorage{m,n})
                        value(m,n) = nan;
                    else
                        points = [obj.ArcStorage{m,n}.Startpoint , ...
                                  obj.ArcStorage{m,n}.Endpoint].';
                        polygonArea = polyarea(real(points(:)), ...
                                               imag(points(:)));
                        arcArea = sum(obj(m,n).Arcs.Area);
                        value(m,n) = polygonArea + arcArea;
                    end
                end
            end
        end

        %% Other functions

        % Unary negative operator
        function r = uminus(obj)
            [M,N] = size(obj);
            r(M,N) = ciat.PolyarcularInterval;
            for m = 1:M
                for n = 1:N
                    r(m,n) = ciat.PolyarcularInterval(-obj(m,n).DefArcs{:});
                end
            end
        end

        % Sample
        function points = sample(obj, nPoints)
            [M,N] = size(obj);
            points = cell(M,N);
            for m = 1:M
                for n = 1:N
                    % Sample arcs and edges
                    arcs = obj.ArcStorage{m,n};
                    edges = ciat.Edge(arcs.Endpoint , ...
                                      circshift(arcs.Startpoint,-1));
                    arcPoints = sample(arcs,nPoints);
                    edgePoints = sample(edges,nPoints);

                    % Interleave arc and edge samples
                    allPoints = [arcPoints.';edgePoints.'];
                    allPoints = allPoints(:);
                    points{m,n} = [allPoints{:}].';
                end
            end

            if M*N == 1
                points = points{:};
            end
        end

        
        % Convex
        function r = isconvex(obj)
            [M,N] = size(obj);
            r(M,N) = false;
            for m = 1:M
                for n = 1:N
                    % Check if all vertices are convex
                    [~,vertexConvexity] = obj(m,n).getVertices;
                    r(m,n) = all(vertexConvexity,'all') && ...
                             all(obj.Arcs{m,n}.Radius >= 0);
                end
            end
        end

        % Convexify
        function outObj = convexify(inObj)
            [M,N] = size(inObj);
            outObj(M,N) = ciat.PolyarcularInterval;
            for m = 1:M
                for n = 1:N
                    arcs = inObj.ArcStorage{m,n};

                    % Create convex hull
                    arcSmp = arcs(arcs.Radius>0).sample(10);
                    arcSmp = [arcSmp{:}];
                    vertexPoly = [arcs.Startpoint ; arcs.Endpoint ; arcSmp(:)];
                    idx = convhull(real(vertexPoly),imag(vertexPoly));
                    edges = ciat.Edge(vertexPoly(idx),...
                                      circshift(vertexPoly(idx),-1));
                    arcs = arcs(arcs.Radius>0);
                    [arcs,edges] = ciat.PolyarcularInterval.splitSegments(arcs,edges);
                    arcs = arcs(arcs.Length > 100*eps);
                    edges = edges(edges.Length > 100*eps);
                    arcs = ciat.PolyarcularInterval.trimSegments(arcs,edges,1);

                    % Generate output object
                    outObj(m,n) = ciat.PolyarcularInterval(arcs);

                    % Join segments
                    % outObj(m,n) = outObj(m,n).joinSegments;
                end
            end
        end

        % Joining segments
        function r = joinSegments(obj)
            arcs = obj.ArcStorage{:};
            K = length(arcs);
            k = 1;
            while k < K
                if arcs(k).Radius == arcs(k+1).Radius && ...
                   arcs(k).Center == arcs(k+1).Center
                    angCup = cup(arcs(k).ArcAngle , arcs(k+1).ArcAngle);
                    if ~isnan(angCup)
                        arcs(k).ArcAngle = angCup;
                        arcs = arcs(setdiff(1:end,k+1));
                        K = K-1;
                    else
                        k = k+1;
                    end
                elseif arcs(k).ArcAngle.Width < 10*eps
                    arcs = arcs(setdiff(1:end,k));
                    K = K-1;
                else
                     k = k+1;
                end
            end
            r = ciat.PolyarcularInterval(arcs);
        end

        % Union
        function r = union(obj) %Needs debugging
            arcs = [obj.Arcs{:}];
            edges = [obj.Edges{:}];
            [arcs,edges] = ciat.PolyarcularInterval.splitSegments(arcs(:),edges(:));
            arcDef = ciat.PolyarcularInterval.trimSegments(arcs,edges,1);
            r = ciat.PolyarcularInterval(arcDef);
            r = joinSegments(r);
        end

        % IsNaN
        function r = isnan(obj)
            r = isnan(obj.Area);
        end 

        % Transpose
        function r = transpose(obj)
            [M,N] = size(obj);
            r = reshape(obj,N,M);
        end

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
            for n = 1:length(obj(:))
                h = [h; obj.Arcs{n}.plotGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Edges{n}.plotGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Vertices{n}.plotGaussMap(arrowSize,varargin{:})];
            end

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
            for n = 1:length(obj(:))
                h = [h; obj.Arcs{n}.plotLogGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Edges{n}.plotLogGaussMap(arrowSize,varargin{:})];
                h = [h; obj.Vertices{n}.plotLogGaussMap(arrowSize,varargin{:})];
            end

            if tf == false 
                hold off
            end
        end

        %% Function headers
        r = sum(obj,varargin)
        r = isin(obj,x)
        [vertices,convexity] = getVertices(obj)
        r = quickSum(obj)

    end % methods

	 %% Static methods
     methods (Static)
        % Function headers
        arcOut = setArcs(arcIn)
        outObj = segmentInverse(obj)
        outObj = segmentProduct(obj1, obj2)
        outObj = cast(inObj,options)
        [arcOut,edgeOut] = splitSegments(arcIn,edgeIn)
        [arcOut,edgeOut] = trimSegments(arcIn,edgeIn,recurCount)
        seg = orderSegments(obj)
        r = plusConvex(obj1,obj2)
        r = plusConcave(obj1,obj2)
     end


     %% Vectorization
     methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.ArcStorage = obj.ArcStorage.(indexOp(1));
            if isscalar(indexOp)
                varargout{1} = obj;
                return;
            end
            [varargout{1:nargout}] = obj.(indexOp(2:end));
        end

        function obj = parenAssign(obj,indexOp,varargin)
            % disp('parenAssign')
            % Ensure object instance is the first argument of call.
            if isempty(obj)
                % This part is for initializing an array of objects
                % such as doing obj(5,2) = ciat.RectangularInterval
                % Might not be the place or the way to do it

                % Instanciate object with zero values of correct size.
                obj = ciat.PolyarcularInterval;
                obj.ArcStorage = cell([indexOp.Indices{:}]);

                % obj = varargin{1};
                varargin{1} = obj.(indexOp);
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.ArcStorage.(indexOp) = rhs.ArcStorage;
                return;
            end
            [obj.(indexOp(2:end))] = varargin{:};
        end

        function n = parenListLength(obj,indexOp,ctx)
            % disp('parenListLength')
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj,indexOp(3:end),ctx);
        end

        function obj = parenDelete(obj,indexOp)
            % disp('parenDelete')
            obj.ArcStorage.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            % disp('cat')
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                newArgs{ix} = varargin{ix}.ArcStorage;
            end
            out = ciat.PolyarcularInterval(cat(dim,newArgs{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.ArcStorage,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            disp('empty')
            obj = ciat.PolyarcularInterval;
        end
    end

    methods
        function obj = reshape(obj,varargin)
            obj.ArcStorage = reshape(obj.ArcStorage,varargin{:});
        end
    end
end