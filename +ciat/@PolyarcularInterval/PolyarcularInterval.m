classdef PolyarcularInterval < matlab.mixin.indexing.RedefinesParen

% Polyarcular interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of complex intervals represented by polyarcular
% regions in the complex plane defined by an ordered series of arcs.
% The object allows performing arithmetic operations on and between them. 
% The object automatically calculates properties of the interval used for 
% casting to other representation types and allows the calculation with 
% arrays and matrices of intervals.
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

	properties (Dependent)
        DefArcs         % Defining arcs of the polygonal interval boundary
        Arcs;           % Implicit arcs 
        Edges;          % Implicit edges
        Segments;       % Boundary segments in counter-clockwise order
        Boundary;       % Boundary segments and vertices in counter-clockwise order
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
        %POLYARCULARINTERVAL Construct an instance of this class
        %
        % This function generates one or more polyarcular intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin an array of arcs as a single argument of
        % ciat.Arc type, which results a single polyarcular interval
        % no matter how the arcs are structured (array or matrix.
        % In order to generate multiple polygonal intervals from arcs
        % a cell array has to be given as a single argument with each cell
        % containing a set of arcs of the ciat.Arc type.
        % It is also possible to give only a single input argument, of one
        % of the other complex interval types in which case it will be 
        % converted to polygonal intervals, representing the smallest 
        % enclosing interval. 
        %__________________________________________________________________________
        % USAGE        
        %   ciat.PolyarxInterval(arx)
        %   ciat.PolyarxInterval(obj)
        %   ciat.PolyarxInterval
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        % _________________________________________________________________________
        % OPTIONS
        %   arx  : struct containing the properties of each convex arc
        % _________________________________________________________________________
        % EXAMPLES
        % _________________________________________________________________________
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

            [M,N] = size(obj);
            if M*N==1
                value = value{:};
            end
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

        % Get all boundary segments
        function value = get.Segments(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    arcs = obj.ArcStorage{m,n};
                    edges = ciat.Edge(arcs.Endpoint , ...
                                      circshift(arcs.Startpoint,-1));
                    % arcs = arcs(arcs.Length > 100*eps);
                    % edges = edges(edges.Length > 100*eps);
                    arcs = mat2cell(arcs,ones(length(arcs),1));
                    edges = mat2cell(edges,ones(length(edges),1));
                    segments = [arcs';edges'];
                    value{m,n} = segments(:);
                end
            end

            if M*N==1
                value = value{:};
            end
        end

        % Get all boundary segments
        function value = get.Boundary(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    % Extract arcs
                    arcs = obj.ArcStorage{m,n};

                    % Extract edges
                    edges = ciat.Edge(arcs.Endpoint , ...
                                      circshift(arcs.Startpoint,-1));

                    % Extract vertices after the arc
                    ang1 = arcs.GaussMap.Supremum .* (arcs.Radius > 0) + ...
                           arcs.GaussMap.Infimum .* (arcs.Radius < 0);
                    ang2 = edges.GaussMap.mid;
                    mask = ang1 > 0 & ang2 < 0;
                    ang1(mask) = ang1(mask) - 2*pi;
                    vertA = ciat.Arc(arcs.Endpoint,zeros(length(arcs),1),...
                                     ciat.RealInterval( ang1,ang2));

                    % Extract vertices after the edge
                    ang1 = edges.GaussMap.mid;
                    ang2 = arcs.GaussMap.Infimum .* (arcs.Radius > 0) + ...
                           arcs.GaussMap.Supremum .* (arcs.Radius < 0);
                    ang2 = circshift(ang2,-1);
                    mask = ang1 > 0 & ang2 < 0;
                    ang1(mask) = ang1(mask) - 2*pi;
                    vertE = ciat.Arc(edges.Endpoint,zeros(length(arcs),1),...
                                     ciat.RealInterval( ang1,ang2));

                    % Create ordered cell array
                    arcs = mat2cell(arcs,ones(length(arcs),1));
                    vertA = mat2cell(vertA,ones(length(vertA),1));
                    edges = mat2cell(edges,ones(length(edges),1));
                    vertE = mat2cell(vertE,ones(length(vertE),1));
                    segments = [arcs.';vertA.';edges.';vertE.'];
                    value{m,n} = segments(:);
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
                    if ~isempty(obj(m,n).Edges)
                        minAbs(m,n) = min(min(inf(abs(obj.ArcStorage{m,n}))),...
                                          min(inf(abs(obj(m,n).Edges))) );
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
                    isOnBoundary = any([obj(m,n).Arcs.ison(0);...
                                         obj(m,n).Edges.ison(0)]);
                    if obj(m,n).isin(0) && ~isOnBoundary
                        minAng(m,n) = -pi;
                        maxAng(m,n) = pi;
                    elseif obj(m,n).Imag.isin(0) && ...
                           obj(m,n).Real.Supremum < 0
                        minAng(m,n) = min(ciat.wrapToPi(inf(angle(arcs)+pi)))+pi;
                        maxAng(m,n) = max(ciat.wrapToPi(sup(angle(arcs)+pi)))+pi;
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
                        arcArea = abs(sum(obj(m,n).Arcs.Area));
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
                    r(m,n) = ciat.PolyarcularInterval(-obj(m,n).DefArcs);
                end
            end
        end

        % Unary reciprocal operator
        function r = recip(obj)
            [M,N] = size(obj);
            r(M,N) = ciat.PolyarcularInterval;
            for m = 1:M
                for n = 1:N
                    inSeg = obj(m,n).Segments;
                    L = length(inSeg);
                    outSeg(L,1) = ciat.Arc;
                    for l = 1:L
                        recipSeg = recip(inSeg{l});
                        if isa(recipSeg,'ciat.Arc')
                            outSeg(l) = recipSeg;
                        end
                    end
                    mask = ~isnan(outSeg) & (outSeg.Radius ~= 0);
                    r(m,n) = ciat.PolyarcularInterval(outSeg(mask));
                end
            end
        end

        % Binary minus operation
        function r = minus(obj1,obj2)
            r = obj1 + (-obj2);
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
                    if isa(allPoints,'cell')
                        points{m,n} = [allPoints{:}].';
                    end
                    % points{m,n} = [allPoints].';
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
                    if ~isempty(obj.Arcs)
                        r(m,n) = all(vertexConvexity,'all') && ...
                                 all(obj(m,n).Arcs.Radius >= 0);
                    else
                        r(m,n) = all(vertexConvexity,'all');
                    end
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
                    vertexPoly = [arcs.Startpoint ; arcs.Endpoint ; arcSmp(:)];
                    idx = convhull(real(vertexPoly),imag(vertexPoly));
                    edges = ciat.Edge(vertexPoly(idx),...
                                      circshift(vertexPoly(idx),-1));
                    arcs = arcs(arcs.Radius>0);
                    [arcs,edges] = ciat.PolyarcularInterval.splitSegments(arcs,edges);
                    arcs = arcs(arcs.Length > 100*eps);
                    edges = edges(edges.Length > 100*eps);
                    arcs = ciat.PolyarcularInterval.trimSegments(arcs,edges,'attempts',5);

                    % Generate output object
                    outObj(m,n) = ciat.PolyarcularInterval(arcs);

                    % Merge arcs
                    % outObj(m,n).mergeArcs;
                end
            end
        end

        % Joining segments
        function r = mergeArcs(obj)
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
        function r = union(obj) % Needs check for empty union
            N = length(obj(:));
            if N == 1
                r = obj;
                return
            end

            r = obj(1);
            for n=2:N
                arcIn = [r.Arcs ; obj(n).Arcs];
                edgeIn = [r.Edges ; obj(n).Edges];
                [arcSp,edgeSp] = ciat.PolyarcularInterval.splitSegments( ...
                                                             arcIn,edgeIn);
                if length(arcSp)+length(edgeSp)>length(arcIn)+length(edgeIn)
                    arcDef = ciat.PolyarcularInterval.trimSegments( ...
                                                             arcSp,edgeSp);
                    r = ciat.PolyarcularInterval(arcDef);
                    r.mergeArcs;
                else
                    warning('The union of intervals is disconnected, returning NaN')
                    r = ciat.PolyarcularInterval;
                end
            end
        end
        % Alias for the union function
        function r = cup(obj,varargin)
            r = union(obj,varargin{:});
        end

        % Intersection
        function r = intersection(obj) % Needs check for empty union
            N = length(obj(:));
            if N == 1
                r = obj;
                return
            end

            r = obj(1);
            for n=2:N
                arcIn = [r.Arcs ; obj(n).Arcs];
                edgeIn = [r.Edges ; obj(n).Edges];
                [arcSp,edgeSp] = ciat.PolyarcularInterval.splitSegments( ...
                                                            arcIn,edgeIn);
                if length(arcSp)+length(edgeSp)>length(arcIn)+length(edgeIn)
                    arcDef = ciat.PolyarcularInterval.trimSegments( ...
                                            arcSp,edgeSp,'inner',true);
                    r = ciat.PolyarcularInterval(arcDef);
                    r.mergeArcs;
                else
                    warning('The intersection of intervals is empty, returning NaN')
                    r = ciat.PolyarcularInterval;
                end
            end
        end
        % Alias for the union function
        function r = cap(obj,varargin)
            r = intersection(obj,varargin{:});
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
                h = [h; obj(n).Arcs.plotGaussMap(arrowSize,varargin{:})];
                h = [h; obj(n).Edges.plotGaussMap(arrowSize,varargin{:})];
                h = [h; obj(n).Vertices.plotGaussMap(arrowSize,varargin{:})];
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
                h = [h; obj(n).Arcs.plotLogGaussMap(arrowSize,varargin{:})];
                h = [h; obj(n).Edges.plotLogGaussMap(arrowSize,varargin{:})];
                h = [h; obj(n).Vertices.plotLogGaussMap(arrowSize,varargin{:})];
            end

            if tf == false 
                hold off
            end
        end

        %% Function headers
        r = plus(obj1,obj2)
        r = sum(obj,varargin)
        r = times(obj1,obj2)
        r = mtimes(obj1,obj2)
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
        [arcOut,edgeOut] = trimSegments(arcIn,edgeIn,optional)
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