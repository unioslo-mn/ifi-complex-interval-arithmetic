classdef PolygonalInterval

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
    
    properties (Dependent)
        Points;         % Vertex points of the polygonal interval boundary
        PointCount;     % Number of vertex points
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
        Area;           % Area of the polygonal interval
    end
    
    properties (Access = private)
       Boundary         % Storage property for the vertex points
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
            obj.Tolerance = optional.tolerance;
            
            switch class(inObj)
                case 'double'
                    if isempty(inObj)
                        % This is for initializing an array of objects
                    else
                        % This is the default way of defining polygonal
                        % intervals the points are assumed to belong to a
                        % single interval no matter how many dimensions
                        obj.Points = inObj(:);
                    end
                case 'cell'
                    % This is how multiple polygons can be defined using
                    % cells of double arrays
                    [M,N] = size(inObj);
                    obj(M,N) = obj;
                    for n = 1:M*N
                        obj(n).Points = inObj{n};
                    end
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
            obj.Boundary = points;
            obj.Boundary = obj.sortPoints;
        end 
        
        % Get points (retrieve from hidden property Boundary)
        function value = get.Points(obj)
            value = obj.Boundary;
        end
        
        %% Dependent properties
                        
        % Get point count
        function value = get.PointCount(obj)
            value = length(obj.Points);
        end
        
        % Real
        function value = get.Real(obj)
            value = ciat.RealInterval( min(real(obj.Points)),...
                                       max(real(obj.Points)) );
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
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
        end
        
        % Imag
        function value = get.Imag(obj)
            value = ciat.RealInterval( min(imag(obj.Points)),...
                                       max(imag(obj.Points)) );
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
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            value = ciat.RealInterval( min(abs(obj.Points)),...
                                       max(abs(obj.Points)) );
            if obj.PointCount >= 3 && ...
                        inpolygon(0,0,real(obj.Points),imag(obj.Points))
                value.Infimum = 0;
            end
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
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            value = ciat.RealInterval( min(angle(obj.Points)),...
                                       max(angle(obj.Points)) );  
            if obj.PointCount >= 3 && ...
                        inpolygon(0,0,real(obj.Points),imag(obj.Points))
                value.Infimum = 0;
                value.Supremum = 2*pi;
            end
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
           value = polyarea(real(obj.Points),imag(obj.Points));
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
        
        % Sum
        function r = sum(obj)
        % Sum of polygonal intervals
        %
        % This function creates the polygonal interval representing the 
        % sum of a set of polygonal intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = sum(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolygonalInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polyInt = sum([ciat.PolygonalInterval([0,1,1i]), ...
        %                    ciat.PolygonalInterval([0,-1,-1i])]);
        % _________________________________________________________________________
            r = obj(1);
            for n = 2:length(obj(:))
                r = r + obj(n);
            end
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
        r = times(obj1,obj2)
        r = mtimes(obj1,obj2)
        points = sortPoints(obj)
        points = backtrack(obj,trackAngle)
                
    end
    
    %% Static methods
    methods (Static)
        % Function headers
        outObj = cast(inObj,options)

        function angleOut = wrap2Pi(angleIn)
            angleOut = rem(2*pi+angleIn, 2*pi);
        end

    end
end

