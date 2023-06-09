classdef RectangularInterval
    
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
                    for m = 1:M
                        for n = 1:N
                            obj(m,n).Real = varargin{1}(m,n);
                            obj(m,n).Imag = varargin{2}(m,n);
                        end
                    end
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
                    M = M1;
                    N = N1;
                    % Create object array
                    obj(M,N) = obj;
                    for m = 1:M
                        for n = 1:N
                            obj(m,n).Real = ciat.RealInterval(...
                                                        varargin{1}(m,n),...
                                                        varargin{2}(m,n));
                            obj(m,n).Imag = ciat.RealInterval(...
                                                        varargin{3}(m,n),...
                                                        varargin{4}(m,n));
                        end
                    end
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
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        %% Dependent properties
        
        % Abs
        function value = get.Abs(obj)
            value = sqrt( abs(obj.Real) * abs(obj.Real) +... 
                          abs(obj.Imag) * abs(obj.Imag) );
            if obj.Real.Infimum <= 0 && obj.Real.Supremum >= 0 && ...
               obj.Imag.Infimum <= 0 && obj.Imag.Supremum >= 0
                value.Infimum = 0;
            end 
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
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            alt = zeros(1,4);
            alt(1) = obj.Real.Infimum  + 1j*obj.Imag.Infimum;
            alt(2) = obj.Real.Infimum  + 1j*obj.Imag.Supremum;
            alt(3) = obj.Real.Supremum + 1j*obj.Imag.Infimum;
            alt(4) = obj.Real.Supremum + 1j*obj.Imag.Supremum;
            value = ciat.RealInterval(min(angle(alt)) , max(angle(alt)));
            
            if obj.Real.Infimum <= 0 && ...
               obj.Imag.Infimum <= 0 && obj.Imag.Supremum >= 0
                value.Infimum = 0;
                value.Supremum = 2*pi;
            end         
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
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Area],M,N);
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
                r = all([obj1.Real] == [obj2.Real]) && all([obj1.Imag] == [obj2.Imag]);
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
            r = any([obj1.Real] ~= [obj2.Real]) && any([obj1.Imag] ~= [obj2.Imag]);
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
            r = obj(1);
            for n = 2:length(obj(:))
                r = r + obj(n);
            end
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
            for n = 1:length(r(:))
                r(n).Real = -r(n).Real;
                r(n).Imag = -r(n).Imag;
            end
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
        
    end
    
    methods (Static)
       outObj = cast(inObj)
    end
end

