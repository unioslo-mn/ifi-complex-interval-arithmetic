classdef PolyarxInterval < matlab.mixin.indexing.RedefinesParen

% Convex polyarc interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of complex intervals represented by convex
% polyarcular regions in the complex plane defined by an ordered set of
% convex arcs.
% The object allows performing arithmetic operations on and between them. 
% In particular, it is optimized for performing addition of polar
% intervals, in the context of beampattern interval analysis.
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



    properties (Access = private)
       Boundary         % Storage property for the vertex points
    end

    properties (Dependent)
        Arx;            % Arx vector defining the polyarx boundary
        Arcs;           % Arc segments of the boundar
        Edges;          % Edge segments of the boundary
        Vertices;       % Vertex segments of the boundary
        ArxCount;       % Number of vertex points
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
        Area;           % Area of the polygonal interval
    end
    
    methods
        
        %% Constructor
        function obj = PolyarxInterval(inObj,inObj2)
        %POLYARXINTERVAL Construct an instance of this class
        %
        % This function generates one or more convex polyarcular intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin a struct containing the properties of the convex 
        % arcs, which results a single polygonal interval.
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
            end    
          
            switch class(inObj)
                case 'double'
                    if isempty(inObj)
                        % This is for initializing an array of objects
                    else
                        obj(1) = ciat.PolyarxInterval;
                        obj.Arx = {inObj};
                    end
                case 'cell'
                    % This is how multiple polygons can be defined using
                    % cells of double arrays
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    obj.Arx = inObj;
                otherwise
                    % Input object will be casted
                    if isempty(inObj2)
                        [M,N] = size(inObj);
                        obj(M,N) = ciat.PolyarxInterval;
                        for n = 1:M*N
                            obj(n) = ciat.PolyarxInterval.cast(inObj(n));
                        end
                    else
                        % Give the option to give only one second object or
                        % the same number than the first object
                        [M,N] = size(inObj);
                        [N2] = length(inObj2(:));
                        assert(N2 == M*N || N2==1)
                        obj(M,N) = ciat.PolyarxInterval;
                        for n = 1:M*N
                            obj(n) = ciat.PolyarxInterval.cast(inObj(n),...
                                                     inObj2(min(n,N2)));
                        end
                    end
            end
        end
                
        %% Defining properties
               
        % Set points (store in the hidden property Boundary after sorting)
        function obj = set.Arx(obj,arx)
            [M,N] = size(obj);
            if M*N == 1
                obj.Boundary = {ciat.PolyarxInterval.sortArx(arx)};
            else
                for m = 1:M
                    for n = 1:N
                        obj.Boundary(m,n) = {ciat.PolyarxInterval.sortArx(arx{m,n})};
                    end
                end
            end
        end 
        
        % Get points (retrieve from hidden property Boundary)
        function value = get.Arx(obj)
            value = obj.Boundary;
            [M,N] = size(obj);
            if M*N == 1
                value = value{:};
            end
        end

        % Get arcs
        function value = get.Arcs(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    arx = obj(m,n).Arx;
                    if size(arx,1)>1
                        arxPrev = circshift(arx,1,1);
                        arxPrev(1,4) = -pi;
                        mask = arx(:,3)~=0;
                        if any(mask)
                            center = complex(arx(mask,1) , arx(mask,2));
                            radius = arx(mask,3);
                            angInf = arxPrev(mask,4);
                            angSup = arx(mask,4);
                            value{m,n} = ciat.Arc(center,radius,...
                                       ciat.RealInterval(angInf,angSup));
                        else
                            value{m,n} = ciat.Arc;
                        end
                    else
                        value{m,n} = ciat.Arc(complex(arx(1),arx(2)),arx(3),...
                                            ciat.RealInterval(-pi,pi));
                    end
                end
            end
            if M*N == 1
                value = value{:};
            end
        end

        % Get edges
        function value = get.Edges(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    arx = obj(m,n).Arx;
                    arxPrev = circshift(arx,1,1);
                    startpoint = complex(arxPrev(:,1),arxPrev(:,2)) + ...
                                 arxPrev(:,3) .* exp(1j*arxPrev(:,4));
                    endpoint =  complex(arx(:,1),arx(:,2)) + ...
                                 arx(:,3) .* exp(1j*arxPrev(:,4));
                    value{m,n} = ciat.Edge(startpoint,endpoint);
                    value{m,n} = value{m,n}(value{m,n}.Length > 10*eps);
                end
            end
            if M*N == 1
                value = value{:};
            end
        end

        % Get vertices
        function value = get.Vertices(obj)
            [M,N] = size(obj);
            value = cell(M,N);
            for m = 1:M
                for n = 1:N
                    arx = obj(m,n).Arx;
                    arxPrev = circshift(arx,1,1);
                    arxPrev(1,4) = -pi;
                    mask = arx(:,3) == 0;
                    center = complex(arx(mask,1),arx(mask,2));
                    angInf = arxPrev(mask,4);
                    angSup = arx(mask,4);
                    value{m,n} = ciat.Arc(center,zeros(size(center)),...
                                          angInf,angSup);
                end
            end
            if M*N == 1
                value = value{:};
            end
        end
        
        %% Dependent properties
                        
        % Get point count
        function value = get.ArxCount(obj)
            if isa(obj,'cell')
                value = cellfun(@length,obj.Arx);
            else
                value = length(obj.Arx);
            end
        end


        % Real
        function value = get.Real(obj)
            [M,N] = size(obj);
            minReal = zeros(M,N);
            maxReal = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    arcs = [obj(m,n).Arcs;obj(m,n).Vertices];
                    minReal(m,n) = min(inf(real(arcs)));
                    maxReal(m,n) = max(sup(real(arcs)));
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
                    arcs = [obj(m,n).Arcs;obj(m,n).Vertices];
                    minImag(m,n) = min(inf(imag(arcs)));
                    maxImag(m,n) = max(sup(imag(arcs)));
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

            if M*N == 1
                arcs = [obj.Arcs; obj.Vertices];
                edges = obj.Edges;
                if ~isempty(edges)
                    minAbs = min([ arcs.Abs.Infimum; ...
                                   edges.Abs.Infimum ] );   
                else
                    minAbs = min(inf(abs(arcs)));
                end
                maxAbs = max(sup(abs(arcs)));
            else
                vert = obj.Vertices;
                arcs = obj.Arcs;
                edges = obj.Edges;
                minAbs = zeros(M,N);
                maxAbs = zeros(M,N);
                for m = 1:M
                    for n = 1:N
                        minAbs(m,n) = min([ vert{m,n}.Abs.Infimum;...
                                            arcs{m,n}.Abs.Infimum; ...
                                            edges{m,n}.Abs.Infimum ] );   
                        maxAbs(m,n) = max([ vert{m,n}.Abs.Infimum;...
                                            arcs{m,n}.Abs.Infimum] ); 
                        
                    end
                end
            end

            % Generate interval
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
                    arcs = [obj(m,n).Arcs;obj(m,n).Vertices];
                    if obj(m,n).isin(0)
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
                    arcs = obj(m,n).Arcs; 
                    vertices = obj(m,n).Vertices; 
                    if isempty(arcs)
                        % value(m,n) = nan;
                        points = vertices.Center;
                        idx = convhull(real(points),imag(points));
                        points = points(idx);
                        value(m,n) = polyarea(real(points(:)), ...
                                               imag(points(:)));
                    else
                        points = [vertices.Center ; ...
                                  arcs.Startpoint ; ...
                                  arcs.Endpoint];
                        idx = convhull(real(points),imag(points));
                        points = points(idx);
                        polygonArea = polyarea(real(points(:)), ...
                                               imag(points(:)));
                        arcArea = sum(arcs.Area);
                        value(m,n) = polygonArea + arcArea;
                    end
                end
            end
        end
        function value = area(obj)
        % Area of polygonal intervals
        %
        % This function returns the area valuea of a set of 
        % polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolyarxInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = abs(ciat.PolyarxInterval([0,1,1i]));
        % _________________________________________________________________________
            value = obj.Area;
        end

        
        
        %% Other functions

        % Sample
        function points = sample(obj, nPoints)
            [M,N] = size(obj);
            points = cell(M,N);
            for m = 1:M
                for n = 1:N
                    % Sample arcs and edges
                    arcs = [obj(m,n).Arcs; obj(m,n).Vertices]; 
                    edges = obj(m,n).Edges;
                    arcPoints = sample(arcs,nPoints);
                    edgePoints = sample(edges,nPoints);

                    % Interleave arc and edge samples
                    allPoints = [[arcPoints{:}].';[edgePoints{:}].'];
                    % allPoints = allPoints(:);
                    % value{m,n} = [allPoints{:}];
                    idx = convhull(real(allPoints),imag(allPoints));
                    points{m,n} = allPoints(idx);
                end
            end

            if M*N == 1
                points = points{:};
            end
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
        % This function plots a set of polygonal intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolyarxInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   h = plot(ciat.PolyarxInterval([0,1,1i]));
        % _________________________________________________________________________
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                if ~isempty(obj(n).Arcs)
                    h = [h;obj(n).Arcs.plot(varargin{:})];    
                end
                if ~isempty(obj(n).Edges)
                    h = [h;obj(n).Edges.plot(varargin{:})];    
                end
            end
            if tf == false 
                hold off
            end
        end
        
        %% Function headers
        r = plus(obj1,obj2)
        r = quickPlus(arx1,arx2)
        r = sum(obj,varargin)
        r = times(obj1,obj2)
        r = mtimes(obj1,obj2)
        r = timesDouble(xObj,d)
        r = isin(obj,x)
        points = backtrack(obj,trackAngle)
                
    end
    
    %% Static methods
    methods (Static)
        % Function headers
        outObj = cast(inObj,options)
        arx = sortArx(arx)
    end

    %% Vectorization
    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Boundary = obj.Boundary.(indexOp(1));
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
                obj = ciat.PolyarxInterval;
                obj.Boundary = cell([indexOp.Indices{:}]);

                % obj = varargin{1};
                varargin{1} = obj.(indexOp);
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.Boundary.(indexOp) = rhs.Boundary;
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
            obj.Boundary.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            % disp('cat')
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                newArgs{ix} = varargin{ix}.Boundary;
            end
            out = ciat.PolyarxInterval(cat(dim,newArgs{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Boundary,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            disp('empty')
            obj = ciat.PolyarxInterval;
        end
    end

    methods
        function obj = reshape(obj,varargin)
            obj.Boundary = reshape(obj.Boundary,varargin{:});
        end
    end
end

