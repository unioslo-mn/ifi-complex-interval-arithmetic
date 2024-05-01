classdef PolygonalInterval < matlab.mixin.indexing.RedefinesParen

% Polygonal interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of complex intervals represented by polygonal
% regions in the complex plane defined by an ordered series of points
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



    properties
        Tolerance;      % Maximum distance from the boundary of the represented interval (except where it is concave)
    end
    
    properties (Access = private)
       Boundary         % Storage property for the vertex points
    end

    properties (Dependent)
        Points;         % Vertex points of the polygonal interval boundary
        PointCount;     % Number of vertex points
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
        Area;           % Area of the polygonal interval
    end
    
    methods
        
        %% Constructor
        function obj = PolygonalInterval(inObj,inObj2,optional)
        %POLYGONALINTERVAL Construct an instance of this class
        %
        % This function generates one or more polygonal intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin an array of points as a single argument of
        % complex double type, which results a single polygonal interval
        % no matter how the points are structured (array or matrix.
        % In order to generate multiple polygonal intervals from points
        % a cell array has to be given as a single argument with each cell
        % containing a set of complex values.
        % It is also possible to give only a single input argument, of one
        % of the other complex interval types in which case it will be 
        % converted to polygonal intervals, representing the smallest 
        % enclosing interval. 
        %__________________________________________________________________________
        % USAGE        
        %   ciat.PolygonalInterval(center,radius)
        %   ciat.PolygonalInterval(obj)
        %   ciat.PolygonalInterval
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        % _________________________________________________________________________
        % OPTIONS
        %   points  : vertex points as complex double type array
        %   cells   : cell array of vertex points 
        %   obj       : object of double or another complex interval type
        %               (see cast function for details)
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = ciat.PolygonalInterval(1,2,3,4)
        %   polyInt = ciat.PolygonalInterval([1,2],[2,3],[3,4],[4,5])
        %   polyInt = ciat.PolygonalInterval([1+1i,2+2i])
        %   polyInt = ciat.PolygonalInterval(ciat.PolarInterval(0,1,2,3))
        %   polyInt(5,2) = ciat.PolygonalInterval
        % _________________________________________________________________________
            arguments
                inObj                (:,:)   = []
                inObj2               (:,:)   = []
                optional.tolerance   (1,1)   {mustBeNumeric}     = 1e-6
            end    

            % obj.Tolerance = optional.tolerance;
            
            switch class(inObj)
                case 'double'
                    if isempty(inObj)
                        % This is for initializing an array of objects
                    else
                        % This is the default way of defining polygonal
                        % intervals the points are assumed to belong to a
                        % single interval no matter how many dimensions
                        obj.Points = {inObj(:)};
                    end
                case 'cell'
                    % This is how multiple polygons can be defined using
                    % cells of double arrays
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    obj.Points = inObj;
                otherwise
                    % Input object will be casted
                    tol = obj.Tolerance;
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    if isempty(inObj2)
                        for n = 1:M*N
                            obj(n) = ciat.PolygonalInterval.cast(inObj(n),...
                                                    'tolerance',tol);
                        end
                    else
                        % Give the option to give only one second object or
                        % the same number than the first object
                        [N2] = length(inObj2(:));
                        assert(N2 == M*N || N2==1)
                        for n = 1:M*N
                            obj(n) = ciat.PolygonalInterval.cast(inObj(n),...
                                                     inObj2(min(n,N2)),...
                                                    'tolerance',tol);
                        end
                    end
            end
        end
                
        %% Defining properties
               
        % Set points (store in the hidden property Boundary after sorting)
        function obj = set.Points(obj,points)
            [M,N] = size(obj);
            for m = 1:M
                for n = 1:N
                    obj.Boundary(m,n) = {ciat.PolygonalInterval.sortPoints( ...
                                          points{m,n},obj.Tolerance(m,n))};
                end
            end
        end 
        
        % Get points (retrieve from hidden property Boundary)
        function value = get.Points(obj)
            value = obj.Boundary;
        end
        
        %% Dependent properties
                        
        % Get point count
        function value = get.PointCount(obj)
            value = cellfun(@length,obj.Points);
        end
        
        % Real
        function value = get.Real(obj)
            [M,N] = size(obj);
            minReal = zeros(M,N);
            maxReal = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    minReal(m,n) = min(real(obj.Points{m,n}));
                    maxReal(m,n) = min(real(obj.Points{m,n}));
                end
            end

            value = ciat.RealInterval( minReal,maxReal );
        end
        function value = real(obj)
        % Real value of polygonal intervals
        %
        % This function creates the real interval representing the 
        % real values of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = real(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = real(ciat.PolygonalInterval([0,1,1i]));
        % _________________________________________________________________________
            value = obj.Real;
        end
        
        % Imag
        function value = get.Imag(obj)
            [M,N] = size(obj);
            minImag = zeros(M,N);
            maxImag = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    minImag(m,n) = min(imag(obj.Points{m,n}));
                    maxImag(m,n) = min(imag(obj.Points{m,n}));
                end
            end

            value = ciat.RealInterval( minImag,maxImag );
        end
        function value = imag(obj)
        % Imaginary value of polygonal intervals
        %
        % This function creates the real interval representing the 
        % imaginary values of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = imag(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = imag(ciat.PolygonalInterval([0,1,1i]));
        % _________________________________________________________________________
            value = obj.Imag;
        end
        
        % Abs
        function value = get.Abs(obj)
            [M,N] = size(obj);
            minAbs = zeros(M,N);
            maxAbs = zeros(M,N);
            pointIn = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    minAbs(m,n) = min(abs(obj.Points{m,n}));
                    maxAbs(m,n) = min(abs(obj.Points{m,n}));
                    pointIn(m,n) = obj.PointCount(m,n) >= 3 && ...
                                  inpolygon(0,0,real(obj.Points{m,n}), ...
                                                imag(obj.Points{m,n}));
                end
            end

            value = ciat.RealInterval( minAbs,maxAbs );
            value(pointIn) = 0;
        end
        function value = abs(obj)
        % Absolute value of polygonal intervals
        %
        % This function creates the real interval representing the 
        % absolute value of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = abs(ciat.PolygonalInterval([0,1,1i]));
        % _________________________________________________________________________
            value = obj.Abs;
        end
        
        % Angle
        function value = get.Angle(obj)
            [M,N] = size(obj);
            minAng = zeros(M,N);
            maxAng = zeros(M,N);
            pointIn = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    minAng(m,n) = min(angle(obj.Points{m,n}));
                    maxAng(m,n) = min(angle(obj.Points{m,n}));
                    pointIn(m,n) = obj.PointCount >= 3 && ...
                                  inpolygon(0,0,real(obj.Points{m,n}), ...
                                                imag(obj.Points{m,n}));
                end
            end

            value = ciat.RealInterval( minAng,maxAng );
            value(pointIn).Infimum = 0;
            value(pointIn).Supremum = 2*pi;
        end
        function value = angle(obj)
        % Angle of polygonal intervals
        %
        % This function creates the real interval representing the 
        % angle of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = abs(ciat.PolygonalInterval([0,1,1i]));
        % _________________________________________________________________________
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
        end
        
        % Area
        function value = get.Area(obj)
            if obj.PointCount == 0
                value = nan;
            else
                value = polyarea(real(obj.Points{:}),imag(obj.Points{:}));
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
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = abs(ciat.PolygonalInterval([0,1,1i]));
        % _________________________________________________________________________
            [M,N] = size(obj);
            value = reshape([obj.Area],M,N);
        end

        
        
        %% Other functions

        % In interval
        function r = ininterval(obj,points)
        % Check if points are in polygonal intervals
        %
        % This function checks if a set of points are in a set of
        % polygonal intervals
        % _________________________________________________________________________
        % USAGE
        %   r = ininterval(obj,points)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj       : array of objects from the ciat.PolygonalInterval class
        %   points    : array of points to be checked
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = ciat.PolygonalInterval([0,1,1i]);
        %   r = ininterval(polyInt,[0.5,0.5i])
        % _________________________________________________________________________
            r = inpolygon(real(points),imag(points),...
                          real(obj.Points),imag(obj.Points));
        end
        
               
        % Subtraction (minus)
        function r = minus(obj1,obj2)
        % Subtraction of polygonal intervals (- operator)
        %
        % This function creates the polygonal interval representing the 
        % difference of two sets of polygonal intervals (see MATLAB minus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = obj1 - obj2
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = ciat.PolygonalInterval([0,1,1i]) - ...
        %             ciat.PolygonalInterval([0,-1,-1i]);
        % _________________________________________________________________________
             r = obj1 + (-obj2);
        end
        
        % Union
        function r = union(obj)
        % Union of polygonal intervals
        %
        % This function creates the polygonal interval representing the 
        % union of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = union(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = union([ciat.PolygonalInterval([0,1,1i]), ...
        %                    ciat.PolygonalInterval([0,-1,-1i])]);
        % _________________________________________________________________________
            r = ciat.PolygonalInterval(cat(1,obj.Points));
        end
        
        % Intersection
        function r = intersection(obj)
        % Intersection of polygonal intervals
        %
        % This function creates the polygonal interval representing the 
        % intersection of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = intersection(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = intersection([ciat.PolygonalInterval([0,1,1i]), ...
        %                    ciat.PolygonalInterval([0,-1,-1i])]);
        % _________________________________________________________________________
            r = obj(1);
            for n = 2:length(obj(:))
                poly1 = polyshape(real(r.Points),imag(r.Points));
                poly2 = polyshape(real(obj(n).Points),imag(obj(n).Points));
                poly3 = intersect(poly1,poly2);
                r.Points = complex(poly3.Vertices(:,1),poly3.Vertices(:,2));
            end
        end
        
        % Negative (uminus)
        function r = uminus(obj)
        % Negative of polygonal intervals (- operator)
        %
        % This function creates the polygonal interval representing the 
        % negative of a set of polygonal intervals (see MATLAB uminus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = -obj
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = -ciat.PolygonalInterval([0,1,1i]);
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                r(n).Points = -r(n).Points;
            end
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
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   h = plot(ciat.PolygonalInterval([0,1,1i]));
        % _________________________________________________________________________
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                points = cat(1,obj(n).Points, obj(n).Points(1));
                h = [h;plot(real(points), imag(points), varargin{:})];    
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
        points = backtrack(obj,trackAngle)
                
    end
    
    %% Static methods
    methods (Static)
        % Function headers
        outObj = cast(inObj,options)
        points = sortPoints(points,tolerance)
    end

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Tolerance = obj.Tolerance.(indexOp(1));
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
                obj = ciat.PolygonalInterval;
                obj.Tolerance = 1e-6 * ones([indexOp.Indices{:}]);
                obj.Boundary = cell([indexOp.Indices{:}]);

                % obj = varargin{1};
                varargin{1} = obj.(indexOp);
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.Tolerance.(indexOp) = rhs.Tolerance;
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
            obj.Tolerance.(indexOp) = [];
            obj.Boundary.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            % disp('cat')
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                newArgs{ix} = varargin{ix}.Tolerance;
                newArgs{ix} = varargin{ix}.Boundary;
            end
            out = ciat.PolygonalInterval(cat(dim,newArgs{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Tolerance,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            disp('empty')
            obj = ciat.PolygonalInterval;
        end
    end
end

