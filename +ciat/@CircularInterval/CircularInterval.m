classdef CircularInterval < matlab.mixin.indexing.RedefinesParen
    
% Circular interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of complex intervals represented by circular
% regions in the complex plane defined by a center and radius. The object 
% allows performing arithmetic operations on and between them. 
% The object automatically calculates properties of the 
% interval used for casting to other representation types and allows the 
% calculation with arrays and matrices of intervals.
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
        Center; % Center point defining the interval
        Radius; % Radius defining the interval
    end
    
    properties (Dependent)
        Real;   % Projection of the circular interval to the real axis
        Imag;   % Projection of the circular interval to the imaginary axis
        Abs;    % Projection of the circular interval to the absolute value axis
        Angle;  % Projection of the circular interval to the angle axis
        Area;   % Area of the circular interval
    end
    
    methods
        %% Constructor
        function obj = CircularInterval(varargin)
        %CIRCULARINTERVAL Construct an instance of this class
        %
        % This function generates one or more circular intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin a center and a radius value as two seperate
        % arguments of double type (it can be an array of matrix). 
        % It is also possible to give only a single input argument,
        % in which case it will be converted to circular intervals, 
        % representing the smallest enclosing interval. In case of 
        % double type input this results in degenerate intervals.
        %__________________________________________________________________________
        % USAGE        
        %   ciat.CircularInterval(center,radius)
        %   ciat.CircularInterval(obj)
        %   ciat.CircularInterval
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        % _________________________________________________________________________
        % OPTIONS
        %   center    : center point as complex double type array
        %   radius    : radius as double type array
        %   obj       : object of double or another complex interval type
        %               (see cast function for details)
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = ciat.CircularInterval(1+1i,1)
        %   circInt = ciat.CircularInterval([1+1i,2+2i],[1,2])
        %   circInt = ciat.CircularInterval([1+1i,2+2i])
        %   circInt = ciat.CircularInterval(ciat.PolarInterval(0,1,2,3))
        %   circInt(5,2) = ciat.CircularInterval
        % _________________________________________________________________________
        
            switch length(varargin)
                    case 0
                        % This is for initializing an array of objects
                    case 1
                        % Single input argument is casted
                        [M,N] = size(varargin{1});
                        obj(M,N) = obj;
                        for n = 1:M*N
                            obj(n) = ciat.CircularInterval.cast(...
                                                        varargin{1}(n));
                        end
                    case 2
                        % Two input arguments is the default way of defining
                        % center and radius for the circular intervals
                        mustBeA(varargin{1},'double')
                        mustBeA(varargin{2},'double')

                        % Get sizes
                        [M1,N1] = size(varargin{1});
                        [M2,N2] = size(varargin{2});
                        assert(M1 == M2 && N1 == N2)
                        M = M1;
                        N = N1;

                        % Create object array
                        obj(M,N) = obj;
                        obj.Center = varargin{1};
                        obj.Radius = varargin{2};
                    otherwise
                        error('Too many input arguments.')
            end
        end
        
        %% Defining properties
        
        % Center
        function value = center(obj)
        % Center point of circular intervals
        %
        % This function returns the center point values of a set of 
        % circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = center(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = center(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
        
            value = obj.Center;
        end
        
        % Radius
        function value = radius(obj)
        % Radius value of circular intervals
        %
        % This function returns the radius values of a set of 
        % circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = radius(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = radius(ciat.CircularInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Radius;
        end
        
        %% Dependent properties
        
        % Real
        function value = get.Real(obj)
            centerReal = real(obj.Center);
            value = ciat.RealInterval(centerReal - obj.Radius,...
                                      centerReal + obj.Radius);
        end
        function value = real(obj)
        % Real value of circular intervals
        %
        % This function creates the real intervals representing the 
        % real values of a set of circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = real(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = real(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
            value = obj.Real;
        end
        
        % Imag
        function value = get.Imag(obj)
            centerImag = imag(obj.Center);
            value = ciat.RealInterval(centerImag - obj.Radius,...
                                      centerImag + obj.Radius);
        end
        function value = imag(obj)
        % Real value of circular intervals
        %
        % This function creates the real intervals representing the 
        % imaginary values of a set of circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = imag(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = imag(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
            value = obj.Imag;
        end
        
        % Abs
        function value = get.Abs(obj)
            value = max(0, ciat.RealInterval(abs(obj.Center) - obj.Radius,...
                                      abs(obj.Center) + obj.Radius));
        end
        function value = abs(obj)
        % Absolute value of circular intervals
        %
        % This function creates the real intervals representing the 
        % absolute values of a set of circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = abs(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
            value = obj.Abs;
        end
        
        % Angle
        function value = get.Angle(obj)
            if abs(obj.Center) <= obj.Radius
                value = ciat.RealInterval(0 , 2*pi);
            else
                centerAngle = angle(obj.Center);
                circleAngle = asin(obj.Radius / abs(obj.Center)); 
                value = ciat.RealInterval(centerAngle - circleAngle,...
                                          centerAngle + circleAngle);
            end         
        end
        function value = angle(obj)
        % Angle value of circular intervals
        %
        % This function creates the real intervals representing the 
        % angle values of a set of circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = angle(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = angle(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
            value = obj.Angle;
        end
        
        % Area
        function value = get.Area(obj)
           value = pi * obj.Radius.^2;
        end
        function value = area(obj)
        % Area value of circular intervals
        %
        % This function returns the area values of a set of 
        % circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = area(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = area(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
            value = obj.Area;
        end
        
        %% Other methods

        function r = eq(obj1,obj2)
            % Equality of circular intervals
            %
            % This function checks if two circular intervals are equal
            % _________________________________________________________________________
            % USAGE
            %   r = eq(obj1,obj2)
            % _________________________________________________________________________
            % NECESSARY ARGUMENTS
            %   obj1      : array of objects from the ciat.CircularInterval class
            %   obj2      : array of objects from the ciat.CircularInterval class
            % _________________________________________________________________________
            % OPTIONS
            % _________________________________________________________________________
            % EXAMPLES
            %   r = eq(ciat.CircularInterval(0,1,2,3),ciat.CircularInterval(0,1,2,3));
            % _________________________________________________________________________
                    [M1,N1] = size(obj1);
                    [M2,N2] = size(obj2);
                    assert(M1 == M2 && N1 == N2, "Arrays have incompatible sizes for this operation.")
                    r = all(obj1.Center == obj2.Center, "all") && ...
                        all(obj1.Radius == obj2.Radius, "all");
        end

        function r = ne(obj1,obj2)
            % Inequality of circular intervals
            %
            % This function checks if two circular intervals are not equal
            % _________________________________________________________________________
            % USAGE
            %   r = ne(obj1,obj2)
            % _________________________________________________________________________
            % NECESSARY ARGUMENTS
            %   obj1      : array of objects from the ciat.CircularInterval class
            %   obj2      : array of objects from the ciat.CircularInterval class
            % _________________________________________________________________________
            % OPTIONS
            % _________________________________________________________________________
            % EXAMPLES
            %   r = ne(ciat.CircularInterval(0,1,2,3),ciat.CircularInterval(0,1,2,3));
            % _________________________________________________________________________
            [M1,N1] = size(obj1);
            [M2,N2] = size(obj2);
            assert(M1 == M2 && N1 == N2)
            r = any(obj1.Center ~= obj2.Center, "all") || ...
                any(obj1.Radius ~= obj2.Radius, "all");
        end

        function r = transpose(obj)
            r = obj;
            r.Center = r.Center.';
            r.Radius = r.Radius.';
        end

        function r = ctranspose(obj)
            r = obj;
            r.Center = r.Center';
            r.Radius = r.Radius.';
        end

        % In interval
        function r = ininterval(obj, points)
        % Check if points are in circular intervals
        %
        % This function checks if a set of points are in a set of
        % circular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = ininterval(obj, points)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj       : array of objects from the ciat.CircularInterval class
        %   points    : array of points in the complex plane
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = ininterval(ciat.CircularInterval(0,1),[0.5+0.5i,1+1i]);
        % _________________________________________________________________________
            r = abs(points - obj.Center) <= obj.Radius;
        end
        
        % Sum
        function r = sum(obj,varargin)
        % Sum of circular intervals
        %
        % This function creates the circular interval representing the 
        % sum of a set of circular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = sum(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = sum([ciat.CircularInterval(0,1), ...
        %                    ciat.CircularInterval(2,3,4,5)]);
        % _________________________________________________________________________
            r = ciat.CircularInterval(sum(obj.Center,varargin{:}), ...
                                      sum(obj.Radius,varargin{:}));
        end
                
        % Negative (uminus)
        function r = uminus(obj)
        % Negative of circular intervals (- operator)
        %
        % This function creates the circular interval representing the 
        % negative of a set of circular intervals (see MATLAB uminus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = -obj
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = -ciat.CircularInterval(0,1);
        % _________________________________________________________________________
            r = ciat.CircularInterval(-obj.Center, obj.Radius);
        end  
        
        % Subtraction (minus)
        function r = minus(obj1,obj2)
        % Subtraction of circular intervals (- operator)
        %
        % This function creates the circular interval representing the 
        % difference of two sets of circular intervals (see MATLAB minus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = obj1 - obj2
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   circInt = ciat.CircularInterval(0,1) - ...
        %             ciat.CircularInterval(2,3);
        % _________________________________________________________________________
            r = obj1 + (-obj2);
        end
                
        % Plot
        function h = plot(obj, varargin)
        % Plot circular intervals 
        %
        % This function plots a set of circular intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.CircularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   h = plot(ciat.CircularInterval(0,1));
        % _________________________________________________________________________
            tf = ishold; 
            if tf == false 
                clf
            end
            hold on
            [M,N] = size(obj); 
            h = [];
            for n = 1:length(obj(:))
                p = obj(n).Center + obj(n).Radius * ...
                                    (exp(1j*(linspace(0,2*pi,360))));
                h = [h;plot(real(p), imag(p), varargin{:})];
            end
            h = reshape(h,M,N);
            if tf == false 
                hold off
            end
        end
        
        %% Function headers
        r = plus(obj1,obj2)
        r = mtimes(obj1,obj2)
        r = times(obj1,obj2)
        r = union(obj)        
    end
    
    methods (Static)
       outObj = cast(inObj)
    end

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Center = obj.Center.(indexOp(1));
            obj.Radius = obj.Radius.(indexOp(1));
            if isscalar(indexOp)
                varargout{1} = obj;
                return;
            end
            [varargout{1:nargout}] = obj.(indexOp(2:end));
        end

        function obj = parenAssign(obj,indexOp,varargin)
            % Ensure object instance is the first argument of call.
            if isempty(obj)
                % This part is for initializing an array of objects
                % such as doing obj(5,2) = ciat.RectangularInterval
                % Might not be the place or the way to do it

                % Instanciate object with zero values of correct size.
                obj = ciat.CircularInterval;
                obj.Center = zeros([indexOp.Indices{:}]);
                obj.Radius = zeros([indexOp.Indices{:}]);

                % obj = varargin{1};
                varargin{1} = obj.(indexOp);
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.Center.(indexOp) = rhs.Center;
                obj.Radius.(indexOp) = rhs.Radius;
                return;
            end
            [obj.(indexOp(2:end))] = varargin{:};
        end

        function n = parenListLength(obj,indexOp,ctx)
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj,indexOp(3:end),ctx);
        end

        function obj = parenDelete(obj,indexOp)
            obj.Center.(indexOp) = [];
            obj.Radius.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.CircularInterval')
                    newArgs{ix} = varargin{ix}.Center;
                    newArgs2{ix} = varargin{ix}.Radius;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ciat.CircularInterval(cat(dim,newArgs{:}), ...
                                           cat(dim,newArgs2{:}));
        end

        function varargout = size(obj,varargin)
            [varargout{1:nargout}] = size(obj.Center,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            disp('empty')
            obj = ciat.CircularInterval;
        end
    end
    
    methods
        function obj = reshape(obj,varargin)
            obj.Center = reshape(obj.Center,varargin{:});
            obj.Radius = reshape(obj.Radius,varargin{:});
        end
    end
end

