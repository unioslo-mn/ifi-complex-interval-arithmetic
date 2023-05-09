classdef CircularInterval
    
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
                        for m = 1:M
                            for n = 1:N
                                obj(m,n).Center = varargin{1}(m,n);
                                obj(m,n).Radius = varargin{2}(m,n);
                            end
                        end
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
        
            [M,N] = size(obj);
            value = reshape([obj.Center],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Radius],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            value = ciat.RealInterval(abs(obj.Center) - obj.Radius,...
                                      abs(obj.Center) + obj.Radius);
            if value.Infimum < 0
                value.Infimum = 0;
            end 
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
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
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
            [M,N] = size(obj);
            value = reshape([obj.Area],M,N);
        end
        
        %% Other methods
        
        % Sum
        function r = sum(obj)
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
            r = obj(1);
            for n = 2:length(obj(:))
                r = r + obj(n);
            end
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
            r = obj;
            for n = 1:length(r(:))
                r(n).Center = -r(n).Center;
            end
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
end

