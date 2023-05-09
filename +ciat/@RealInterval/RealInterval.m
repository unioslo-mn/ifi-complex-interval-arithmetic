classdef RealInterval

% Real interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of real intervals defined by two real values.
% This class allows the representation of certain parameters of complex
% interval classes (such as real and imaginary intervals absolute value 
% or angle. The object allows performing arithmetic operations on and 
% between them. The object automatically calculates properties of the 
% interval used for arithmetic operations and casting and allows the 
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
        Infimum         % Lower bound of the interval
        Supremum        % Upper bound of the interval
    end
    
    properties (Dependent)
        Midpoint        % Value half-way between the infimum and supremum
        Width           % Difference of the infimum and supremum
        Bounds          % Infimum and supremum as an array
    end
    
    methods
        %% Constructor
        function obj = RealInterval(varargin)
        %REALINTERVAL Construct an instance of this class
        %
        % This function generates one or more real intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin an infimum and a supremum value as two seperate
        % arguments of real interval type (it can be an array or matrix). 
        % Another default method is giving a single argument containing
        % an array with one size of 2, containing double type values, in
        % which case the first value along the dimension of size 2 will be
        % the infimum and the second the supremum.
        %__________________________________________________________________________
        % USAGE        
        %   ciat.RealInterval(center,radius)
        %   ciat.RealInterval(obj)
        %   ciat.RealInterval
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        % _________________________________________________________________________
        % OPTIONS
        %   infimum    : infimum value of the intervals (double type array)
        %   supremum   : supremum value of the intervals (double type array)
        %   obj        : array of double type with one dimension size of 2 
        % _________________________________________________________________________
        % EXAMPLES
        %   rectInt = ciat.RealInterval(1,2)
        %   rectInt = ciat.RealInterval([1,2],[2,3])
        %   rectInt = ciat.RealInterval([1,2,3 ; 2,3,4])
        %   rectInt = ciat.RealInterval([1,2 ; 3,4 ; 5,6])
        %   rectInt(5,2) = ciat.RealInterval
        % _________________________________________________________________________
            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 1
                    mustBeNumeric(varargin{1});
                    [M,N] = size(varargin{1});
                    switch min(M,N)
                        case 1
                            obj.Infimum = min(varargin{1});
                            obj.Supremum = max(varargin{1});
                        case 2
                            M = max([M,N]);
                            infsup = reshape(varargin{1},[],2);
                            infsup = sort(infsup,2);
                            obj(M,1) = obj;
                            for m = 1:M
                                obj(m).Infimum = infsup(m,1);
                                obj(m).Supremum = infsup(m,2);
                            end
                        otherwise
                            error('Invalid input array size (>2).')
                    end
                case 2
                    mustBeNumeric(varargin{1});
                    mustBeNumeric(varargin{2});
                    assert(size(varargin{1},1) == size(varargin{2},1))
                    assert(size(varargin{1},2) == size(varargin{2},2))
                    [M,N] = size(varargin{1});
                    obj(M,N) = obj;
                   for n = 1:M*N
                        obj(n).Infimum = min(varargin{1}(n),varargin{2}(n));
                        obj(n).Supremum = max(varargin{1}(n),varargin{2}(n));
                   end
            end
        end
        
        %% Definint properties
        
        % Infimum
        function value = inf(objArr)
        % Infimum of real intervals
        %
        % This function returns the infimum values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = inf(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   val = inf(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            [M,N] = size(objArr);
            value = reshape([objArr.Infimum],M,N);
        end
        
        % Supremum
        function value = sup(objArr)
        % Supremum of real intervals
        %
        % This function returns the supremum values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = sup(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   val = sup(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            [M,N] = size(objArr);
            value = reshape([objArr.Supremum],M,N);
        end
        
        % Midpoint
        function value = mid(objArr)
        % Midpoint of real intervals
        %
        % This function returns the midpoint values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = mid(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   val = mid(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            [M,N] = size(objArr);
            value = reshape([objArr.Midpoint],M,N);
        end
        
        % Width
        function value = width(objArr)
        % Width of real intervals
        %
        % This function returns the width values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = width(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   val = width(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            [M,N] = size(objArr);
            value = reshape([objArr.Width],M,N);
        end
        
        %% Dependent properties
        
        %Midpoint
        function value = get.Midpoint(obj)
            value = (obj.Supremum + obj.Infimum)/2;
        end
        
        % Width
        function value = get.Width(obj)
            value = obj.Supremum - obj.Infimum;
        end
        
        % Bounds as a 2 element array
        function value = get.Bounds(obj)
            value = [obj.Infimum;obj.Supremum];
        end
        
        %% Other methods
        
        % Sum
        function r = sum(obj)
        % Sum of real intervals
        %
        % This function creates the real interval representing the 
        % sum of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = sum(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = sum([ciat.RealInterval(0,1), ...
        %                    ciat.RealInterval(2,3,4,5)]);
        % _________________________________________________________________________
            r = obj(1);
            for n = 2:length(obj(:))
                r = r + obj(n);
            end
        end
        
        % Negative (uminus)
        function r = uminus(obj)
        % Negative of real intervals (- operator)
        %
        % This function creates the real interval representing the 
        % negative of a set of real intervals (see MATLAB uminus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = -obj
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = -ciat.RealInterval(0,1);
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                r(n).Infimum = -r(n).Supremum;
                r(n).Supremum = -obj(n).Infimum;
            end
        end
        
        % Subtraction (minus)
        function r = minus(obj1,obj2)
        % Subtraction of real intervals (- operator)
        %
        % This function creates the real interval representing the 
        % difference of two sets of real intervals (see MATLAB minus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = obj1 - obj2
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = ciat.RealInterval(0,1) - ...
        %             ciat.RealInterval(2,3);
        % _________________________________________________________________________
            r = obj1 + (-obj2);
        end
        
        % Reciprocal
        function r = recip(obj)
        % Reciprocal of real intervals
        %
        % This function creates the real interval representing the 
        % reciprocal of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = recip(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = recip(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                if obj(n).Infimum <=0 && obj(n).Supremum>=0
                    warning('Zero division')
                    r(n).Infimum = nan();
                    r(n).Supremum = nan();
                else
                    r(n).Infimum  = 1/obj(n).Supremum;
                    r(n).Supremum = 1/obj(n).Infimum;
                end
            end
        end
        
        % Absolute value
        function r = abs(obj)
        % Absolute value of real intervals
        %
        % This function creates the real intervals representing the 
        % absolute values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = abs(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                r(n).Infimum = min(abs(obj(n).Bounds));
                r(n).Supremum = max(abs(obj(n).Bounds));
                if (obj(n).Infimum * obj(n).Supremum) < 0
                    r(n).Infimum = 0;
                end
            end
        end
        
        % Exponential
        function r = exp(obj)
        % Exponential value of real intervals
        %
        % This function creates the real intervals representing the 
        % exponential values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = exp(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = exp(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                r(n).Infimum = exp(obj(n).Infimum);
                r(n).Supremum = exp(obj(n).Supremum);
            end
        end
        
        % Logarithm
        function r = log(obj)
        % Logarithm value of real intervals
        %
        % This function creates the real intervals representing the 
        % logarithm values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = log(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = log(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                if obj(n).Supremum >= 0
                    r(n).Supremum = log(obj(n).Supremum);
                    if obj(n).Infimum >= 0
                        r(n).Infimum = log(obj(n).Infimum);
                    else
                        warning('Log of negative infimum replaced by zero.')
                        r(n).Infimum = 0;
                    end
                else
                    warning('Log of negative interval replaced by zero.')
                    r(n).Infimum = 0;
                    r(n).Supremum = 0;
                end
            end
        end
        function r = log10(obj)
        % 10-base logarithm value of real intervals
        %
        % This function creates the real intervals representing the 
        % 10-base logarithm values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = log10(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = log10(ciat.RealInterval(0,1));
        % _________________________________________________________________________
           r = obj;
            for n = 1:length(r(:))
                if obj(n).Supremum >= 0
                    r(n).Supremum = log10(obj(n).Supremum);
                    if obj(n).Infimum >= 0
                        r(n).Infimum = log10(obj(n).Infimum);
                    else
                        warning('Log of negative infimum replaced by zero.')
                        r(n).Infimum = 0;
                    end
                else
                    warning('Log of negative interval replaced by zero.')
                    r(n).Infimum = 0;
                    r(n).Supremum = 0;
                end
            end
        end
        
        % Square-root
        function r = sqrt(obj)
        % Square-root value of real intervals
        %
        % This function creates the real intervals representing the 
        % square-root values of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = sqrt(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = sqrt(ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                if obj(n).Supremum >= 0
                    r(n).Supremum = sqrt(obj(n).Supremum);
                    if obj(n).Infimum >= 0
                        r(n).Infimum = sqrt(obj(n).Infimum);
                    else
                        warning('Sqrt of negative infimum replaced by zero.')
                        r(n).Infimum = 0;
                    end
                else
                    warning('Sqrt of negative interval replaced by zero.')
                    r(n).Infimum = 0;
                    r(n).Supremum = 0;
                end
            end
        end
        
        % Union
        function r = union(obj)
        % Union of real intervals
        %
        % This function creates the real interval representing the 
        % union of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = union(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = union([ciat.RealInterval(0,1,2,3), ...
        %                    ciat.RealInterval(2,3,4,5)]);
        % _________________________________________________________________________
            N = length(obj(:));
            assert(N>1)
            r = ciat.RealInterval( min([obj.Infimum]) , ...
                                   max([obj.Supremum]) );
        end
        
        % Intersection
        function r = intersection(obj)
        % Intersection of real intervals
        %
        % This function creates the real interval representing the 
        % intersection of a set of real intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = intersection(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   realInt = intersection([ciat.RealInterval(0,1,2,3), ...
        %                    ciat.RealInterval(2,3,4,5)]);
        % _________________________________________________________________________
        	N = length(obj(:));
            assert(N>1)
            maxInf = max([obj.Infimum]);
            minSup = min([obj.Supremum]);
            if maxInf <= minSup
                r = ciat.RealInterval( maxInf , minSup );
            else
                r = ciat.RealInterval.empty;
            end
        end
        
        % Plot
        function plt = plot(obj,varargin)
        % Plot real intervals 
        %
        % This function plots a set of real intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   h = plot(ciat.RealInterval(0,1,2,3));
        % _________________________________________________________________________
            [M,N] = size(obj);
            inf = [obj.Infimum];
            sup = [obj.Supremum];
            plt = line(repmat(1:M*N,2,1),[inf;sup],varargin{:});
        end
          
        
        %% Method signatures
        r = plus(obj1,obj2)
        r = times(obj1,obj2)
        r = mtimes(obj1,obj2)
        r = rdivide(obj1,obj2)
        r = mrdivide(obj1,obj2)
        r = sin(obj)
        r = cos(obj)
                
    end
end

