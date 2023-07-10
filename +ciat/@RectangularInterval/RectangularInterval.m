classdef RectangularInterval < matlab.mixin.indexing.RedefinesParen
    
% Rectangular interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of complex intervals represented by rectangular
% regions in the complex plane defined by an interval along the real and 
% one along the imaginary axis (two real intervals).
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
        Real;    % Projection of the rectangular interval to the real axis
        Imag;   % Projection of the rectangular interval to the imaginary axis
    end
    
     properties (Dependent)
        Abs;    % Projection of the rectangular interval to the absolute value axis
        Angle;  % Projection of the rectangular interval to the angle axis
        Area;   % Area of the rectangular interval
    end

        
    
    methods
        % Class creator
        function obj = RectangularInterval(varargin)
        %RECTANGULARINTERVAL Construct an instance of this class
        %
        % This function generates one or more rectangular intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin a real and an imaginary interval as two seperate
        % arguments of real interval type (it can be an array or matrix). 
        % Another default method is giving the real infimum, real supremum,
        % imaginary infimum and imaginary supremum values as four seperate
        % arguments of double type (can be an array or matrix).
        % It is also possible to give only a single input argument,
        % in which case it will be converted to rectangular intervals, 
        % representing the smallest enclosing interval. In case of 
        % double type input this results in degenerate intervals.
        %__________________________________________________________________________
        % USAGE        
        %   ciat.RectangularInterval(center,radius)
        %   ciat.RectangularInterval(obj)
        %   ciat.RectangularInterval
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        % _________________________________________________________________________
        % OPTIONS
        %   real    : real value as real interval type array
        %   imag    : imaginary value as real interval type array
        %   realInf : real infimum value as double type array
        %   realSup : real supremum value as double type array
        %   imagInf : imaginary infimum value as double type array
        %   imagSup : imaginary supremum value as double type array
        %   obj       : object of double or another complex interval type
        %               (see cast function for details)
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = ciat.RectangularInterval(1,2,3,4)
        %   rectInt = ciat.RectangularInterval([1,2],[2,3],[3,4],[4,5])
        %   rectInt = ciat.RectangularInterval([1+1i,2+2i])
        %   rectInt = ciat.RectangularInterval(ciat.PolarInterval(0,1,2,3))
        %   rectInt(5,2) = ciat.RectangularInterval
        % _________________________________________________________________________
            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 1
                    % Single input argument is casted
                    [M,N] = size(varargin{1});
                    obj(M,N) = obj;
                    for n = 1:M*N
                        obj(n) = ciat.RectangularInterval.cast(...
                                                        varargin{1}(n));
                    end
                case 2
                    % Two input arguments of real intervals is the default
                    % way of definin rectangular intervals
                    mustBeA(varargin{1},'ciat.RealInterval')
                    mustBeA(varargin{2},'ciat.RealInterval')
                    
                    % Get sizes
                    [M1,N1] = size(varargin{1});
                    [M2,N2] = size(varargin{2});
                    assert(M1 == M2 && N1 == N2)
                    M = M1;
                    N = N1;
                    
                    % Create object array
                    obj(M,N) = obj;
                    obj.Real = varargin{1};
                    obj.Imag = varargin{2};
                case 4
                    % Four input arguments of doubles is also a default way
                    % of defining rectangular intervals 
                    mustBeA(varargin{1},'double')
                    mustBeA(varargin{2},'double')
                    mustBeA(varargin{3},'double')
                    mustBeA(varargin{4},'double')
                    
                    % Get sizes
                    [M1,N1] = size(varargin{1});
                    [M2,N2] = size(varargin{2});
                    [M3,N3] = size(varargin{3});
                    [M4,N4] = size(varargin{4});
                    assert(var([M1,M2,M3,M4]) == 0 && ...
                           var([N1,N2,N3,N4]) == 0)

                    obj.Real = ciat.RealInterval(varargin{1},varargin{2});
                    obj.Imag = ciat.RealInterval(varargin{3},varargin{4});
            end
            
        end
        
        %% Defining properties
        % Real
        function value = real(obj)
        % Real value of rectangular intervals
        %
        % This function creates the real interval representing the 
        % real values of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = real(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = real(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Real;
        end
        
        % Imag
        function value = imag(obj)
        % Imaginary value of rectangular intervals
        %
        % This function creates the real interval representing the 
        % imaginary values of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = imag(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = imag(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Imag;
        end
        
        %% Dependent properties
        
        % Abs
        function value = get.Abs(obj)
            value = sqrt( abs(obj.Real) * abs(obj.Real) +... 
                          abs(obj.Imag) * abs(obj.Imag) );
            value.Infimum(obj.contains(0)) = 0;
        end
        function value = abs(obj)
        % Absolute value of rectangular intervals
        %
        % This function creates the real interval representing the 
        % absolute value of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = abs(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Abs;
        end
        
        % Angle
        function value = get.Angle(obj)
            alt = zeros(1,4);
            alt(1) = obj.Real.Infimum  + 1j*obj.Imag.Infimum;
            alt(2) = obj.Real.Infimum  + 1j*obj.Imag.Supremum;
            alt(3) = obj.Real.Supremum + 1j*obj.Imag.Infimum;
            alt(4) = obj.Real.Supremum + 1j*obj.Imag.Supremum;
            value = ciat.RealInterval(min(angle(alt)) , max(angle(alt)));
            
            value.Infimum(obj.Real.Infimum <= 0 & ...
                          obj.Imag.Infimum <= 0 & ...
                          obj.Imag.Supremum >= 0) = 0;
            value.Supremum(obj.Real.Infimum <= 0 & ...
                           obj.Imag.Infimum <= 0 & ...
                           obj.Imag.Supremum >= 0) = 2*pi;
        end
        function value = angle(obj)
        % Angle of rectangular intervals
        %
        % This function creates the real interval representing the 
        % angle of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = abs(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Angle;
        end
        
        % Area
        function value = get.Area(obj)
           value = width(obj.Real) .* width(obj.Imag);
        end        
        function value = area(obj)
        % Area of rectangular intervals
        %
        % This function returns the area valuea of a set of 
        % rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = abs(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Area;
        end
        
        %% Other methods
        
        % Equal
        function r = eq(obj1,obj2)
        % Equality rectangular intervals
        %
        % This function checks if two sets of rectangular intervals are equal
        % _________________________________________________________________________
        % USAGE
        %   r = eq(obj1,obj2)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj1      : array of objects from the ciat.RectangularInterval class
        %   obj2      : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = eq(ciat.RectangularInterval(0,1,2,3), ...
        %          ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            [M1,N1] = size(obj1);
            [M2,N2] = size(obj2);
            assert(M1 == M2 && N1 == N2)
            r = all(obj1.Real == obj2.Real, "all") && ...
                all(obj1.Imag == obj2.Imag, "all");
        end

        % Not equal
        function r = ne(obj1, obj2)
        % Not equal rectangular intervals
        %
        % This function checks if two sets of rectangular intervals are not equal
        % _________________________________________________________________________
        % USAGE
        %   r = ne(obj1,obj2)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj1      : array of objects from the ciat.RectangularInterval class
        %   obj2      : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   r = ne(ciat.RectangularInterval(0,1,2,3), ...
        %          ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            [M1,N1] = size(obj1);
            [M2,N2] = size(obj2);
            assert(M1 == M2 && N1 == N2)
            r = any(obj1.Real ~= obj2.Real, "all") || ...
                any(obj1.Imag ~= obj2.Imag, "all");
        end

        function r = transpose(obj)
            r = obj;
            r.Real = r.Real.';
            r.Imag = r.Imag.';
        end

        function r = ctranspose(obj)
            r = obj;
            r.Real = r.Real.';
            r.Imag = -r.Imag.';
        end

        function r = contains(obj, x)
        % Contains rectangular intervals
        %
        % This function checks if a set of rectangular intervals contains
        % a given value. Returns a logical array of the same size as the
        % input array.
        % _________________________________________________________________________
        % USAGE
        %   r = contains(obj, x)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj       : array of objects from the ciat.RectangularInterval class
        %   x         : complex value
        % _________________________________________________________________________
        % EXAMPLES
        %   r = contains(ciat.RectangularInterval(0,1,2,3), 0.5);
        % _________________________________________________________________________
            r = obj.Real.contains(real(x)) & obj.Imag.contains(imag(x));
        end

        % Sum
        function r = sum(obj)
        % Sum of rectangular intervals
        %
        % This function creates the rectangular interval representing the 
        % sum of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = sum(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = sum([ciat.RectangularInterval(0,1,2,3), ...
        %                    ciat.RectangularInterval(2,3,4,5)]);
        % _________________________________________________________________________
            r = ciat.RectangularInterval(sum([obj.Real]) , ...
                                         sum([obj.Imag]));
        end
        
        % Subtraction (minus)
        function r = minus(obj1,obj2)
        % Subtraction of rectangular intervals (- operator)
        %
        % This function creates the rectangular interval representing the 
        % difference of two sets of rectangular intervals (see MATLAB minus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = obj1 - obj2
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = ciat.RectangularInterval(0,1,2,3) - ...
        %             ciat.RectangularInterval(0,1,2,3);
        % _________________________________________________________________________
            r = obj1 + (-obj2);
        end

        % Exponential
        function r = exp(obj)
        % Exponential value of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % exponential value of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = exp(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = exp(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = ciat.RectangularInterval(exp(obj.Real).*cos(obj.Imag), ...
                                         exp(obj.Real).*sin(obj.Imag));
        end

        % Sine
        function r = sin(obj)
        % Sine of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % sine of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = sin(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = sin(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = ciat.RectangularInterval(sin(obj.Real).*cosh(obj.Imag), cos(obj.Real).*sinh(obj.Imag));
        end

        % Cosine
        function r = cos(obj)
        % Cosine of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % cosine of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = cos(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = cos(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = ciat.RectangularInterval(cos(obj.Real).*cosh(obj.Imag), -sin(obj.Real).*sinh(obj.Imag));
        end

        % Cotangent
        function r = cot(obj)
        % Cotangent of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % cotangent of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = cot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = cot(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = tan(pi/2 - obj);
        end

        % Hyperbolic sine
        function r = sinh(obj)
        % Hyperbolic sine of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % hyperbolic sine of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = sinh(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = sinh(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = -1j.*sin(1j.*obj);
        end

        % Hyperbolic cosine
        function r = cosh(obj)
        % Hyperbolic cosine of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % hyperbolic cosine of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = cosh(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = cosh(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = cos(1j.*obj);
        end

        % Hyperbolic tangent
        function r = tanh(obj)
        % Hyperbolic tangent of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % hyperbolic tangent of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = tanh(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = tanh(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = -1j.*tan(1j.*obj);
        end

        % Hyperbolic cotangent
        function r = coth(obj)
        % Hyperbolic cotangent of rectangular intervals
        %
        % This function creates the rectangular interval representing the
        % hyperbolic cotangent of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE
        %   r = coth(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = coth(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            r = 1j.*cot(1j.*obj);
        end
        
        % Union
        function r = union(obj)
        % Union of rectangular intervals
        %
        % This function creates the rectangular interval representing the 
        % union of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = union(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = union([ciat.RectangularInterval(0,1,2,3), ...
        %                    ciat.RectangularInterval(2,3,4,5)]);
        % _________________________________________________________________________
            N = length(obj(:));
            assert(N>1)
            r = ciat.RectangularInterval(union([obj.Real]) , ...
                                         union([obj.Imag]));
        end
        
        % Intersection
        function r = intersection(obj)
        % Intersection of rectangular intervals
        %
        % This function creates the rectangular interval representing the 
        % intersection of a set of rectangular intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = intersection(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = intersection([ciat.RectangularInterval(0,1,2,3), ...
        %                    ciat.RectangularInterval(2,3,4,5)]);
        % _________________________________________________________________________
            N = length(obj(:));
            assert(N>1)
            r = ciat.RectangularInterval(intersection([obj.Real]) , ...
                                         intersection([obj.Imag]));
        end
        
        % Negative (uminus)
        function r = uminus(obj)
        % Negative of rectangular intervals (- operator)
        %
        % This function creates the rectangular interval representing the 
        % negative of a set of rectangular intervals (see MATLAB uminus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = -obj
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = -ciat.RectangularInterval(0,1,2,3);
        % _________________________________________________________________________
            r = obj;
            r.Real = -r.Real;
            r.Imag = -r.Imag;
        end  
        
        % Plot
        function h = plot(obj, varargin)
        % Plot rectangular intervals 
        %
        % This function plots a set of rectangular intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RectangularInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   h = plot(ciat.RectangularInterval(0,1,2,3));
        % _________________________________________________________________________
            tf = ishold; 
            hold on
            [M,N] = size(obj); 
            h = [];
            for n = 1:length(obj(:))
                p1 = obj(n).Real.Infimum + 1j*obj(n).Imag.Infimum;
                p2 = obj(n).Real.Infimum + 1j*obj(n).Imag.Supremum;
                p3 = obj(n).Real.Supremum + 1j*obj(n).Imag.Supremum;
                p4 = obj(n).Real.Supremum + 1j*obj(n).Imag.Infimum;

                h = [h;plot(real([p1,p2,p3,p4,p1]), ...
                            imag([p1,p2,p3,p4,p1]), varargin{:})];
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
        r = recip(obj)
        r = rdivide(obj1,obj2)
        r = power(obj1,obj2)
        r = sqrt(obj)
        
    end
    
    methods (Static)
       outObj = cast(inObj)
    end

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Real = obj.Real.(indexOp(1));
            obj.Imag = obj.Imag.(indexOp(1));
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
                obj = ciat.RectangularInterval;
                obj.Real = ciat.RealInterval(zeros([indexOp.Indices{:}]), zeros([indexOp.Indices{:}]));
                obj.Imag = ciat.RealInterval(zeros([indexOp.Indices{:}]), zeros([indexOp.Indices{:}]));

                % obj = varargin{1};
                varargin{1} = obj.(indexOp);
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.Real.(indexOp) = rhs.Real;
                obj.Imag.(indexOp) = rhs.Imag;
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
            obj.Real.(indexOp) = [];
            obj.Imag.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.RectangularInterval')
                    newArgs{ix} = varargin{ix}.Real;
                    newArgs2{ix} = varargin{ix}.Imag;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ciat.RectangularInterval(cat(dim,newArgs{:}), ...
                                           cat(dim,newArgs2{:}));
        end

        function varargout = size(obj,varargin)
            [varargout{1:nargout}] = size(obj.Real,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            disp('empty')
            obj = ciat.RectangularInterval;
        end
    end

    methods
        function obj = reshape(obj,varargin)
            obj.Real = reshape(obj.Real,varargin{:});
            obj.Imag = reshape(obj.Imag,varargin{:});
        end
    end
end

