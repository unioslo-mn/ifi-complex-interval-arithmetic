classdef RealInterval < matlab.mixin.indexing.RedefinesParen

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
        %   ciat.RealInterval(infimum,supremum)
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
            
            for varIdx = 1:length(varargin)
                assert(isreal(varargin{varIdx}),'input must be real valued');
            end
            
            switch length(varargin)
                case 0
                % This is for initializing an array of objects
                case 1
                mustBeNumeric(varargin{1});
                [M,N] = size(varargin{1});
                switch min(M,N)
                    case 1
                    % obj.Infimum = min(varargin{1});
                    % obj.Supremum = max(varargin{1});
                    obj.Infimum = varargin{1};
                    obj.Supremum = varargin{1};
                    case 2
                    [M,dimIdx] = max([M,N]);
                    if dimIdx == 2
                        varargin{1} = varargin{1}.';
                    end
                    infsup = reshape(varargin{1},[],2);
                    infsup = sort(infsup,2);
                    obj(M,1) = obj;
                    for m = 1:M
                        obj(m).Infimum = infsup(m,1);
                        obj(m).Supremum = infsup(m,2);
                    end
                    if M == 2
                        warning(['Ambiguous input structure',...
                        'we assume each row to be'...,
                        'a seperate interval.'])
                    end
                    otherwise
                        obj.Infimum = varargin{1};
                        obj.Supremum = varargin{1};
                        % error('Invalid input array size (>2).')
                end
                case 2
                mustBeNumeric(varargin{1});
                mustBeNumeric(varargin{2});
                assert(size(varargin{1},1) == size(varargin{2},1))
                assert(size(varargin{1},2) == size(varargin{2},2))
                
                % obj.Infimum = min(varargin{1},varargin{2});
                % obj.Supremum = max(varargin{1},varargin{2});
                obj.Infimum = varargin{1};
                obj.Supremum = varargin{2};
            end
        end
        
        %% Defining properties
        
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
            value = objArr.Infimum;
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
            value = objArr.Supremum;
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
        
        % Equality
        function r = eq(obj1,obj2)
        % Equality of real intervals
        %
        % This function returns true if two real intervals are equal
        % _________________________________________________________________________
        % USAGE
        %   r = eq(obj1,obj2)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj1      : array of objects from the ciat.RealInterval class
        %   obj2      : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = eq(ciat.RealInterval(0,1),ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = all(obj1.Infimum == obj2.Infimum, "all") && ...
                all(obj1.Supremum == obj2.Supremum, "all");
        end
        
        % Inequality
        function r = ne(obj1,obj2)
        % Inequality of real intervals
        %
        % This function returns true if two real intervals are not equal
        % _________________________________________________________________________
        % USAGE
        %   r = ne(obj1,obj2)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj1      : array of objects from the ciat.RealInterval class
        %   obj2      : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = ne(ciat.RealInterval(0,1),ciat.RealInterval(0,1));
        % _________________________________________________________________________
            r = any(obj1.Infimum ~= obj2.Infimum, "all") || ...
            any(obj1.Supremum ~= obj2.Supremum, "all");
        end

        function r = transpose(obj)
            r = obj;
            r.Infimum = r.Infimum.';
            r.Supremum = r.Supremum.';
        end

        function r = ctranspose(obj)
            r = obj.';
        end
        
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
            % r = obj(1);
            % for n = 2:length(obj(:))
            %     r = r + obj(n);
            % end
            r = ciat.RealInterval(sum([obj.Infimum]), sum([obj.Supremum]));
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
            r = ciat.RealInterval(-obj.Supremum,-obj.Infimum);
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
            r.Infimum = exp(obj.Infimum);
            r.Supremum = exp(obj.Supremum);
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

        % Minimum
        function r = min(obj1, obj2)
        % Minimum of real intervals
        %
        % This function creates the real interval representing the
        % minimum of two sets of real intervals
        % _________________________________________________________________________
        % USAGE
        %   r = min(obj1, obj2)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj1      : array of objects from the ciat.RealInterval class
        %   obj2      : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = min(ciat.RealInterval(0,1), ciat.RealInterval(2,3));
        % _________________________________________________________________________
            % Turn scalars to degenerate intervals
            if isa(obj1, 'ciat.RealInterval') == 0
                        obj1 = ciat.RealInterval(obj1, obj1);
            end
            if isa(obj2, 'ciat.RealInterval') == 0
                obj2 = ciat.RealInterval(obj2, obj2);
            end
            r = ciat.RealInterval(min(obj1.Infimum, obj2.Infimum), ...
                min(obj1.Supremum, obj2.Supremum));
        end

        % Maximum
        function r = max(obj1, obj2)
        % Maximum of real intervals
        %
        % This function creates the real interval representing the
        % maximum of two sets of real intervals
        % _________________________________________________________________________
        % USAGE
        %   r = max(obj1, obj2)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj1      : array of objects from the ciat.RealInterval class
        %   obj2      : array of objects from the ciat.RealInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = max(ciat.RealInterval(0,1), ciat.RealInterval(2,3));
        % _________________________________________________________________________
            % Turn scalars to degenerate intervals
            if isa(obj1, 'ciat.RealInterval') == 0
                        obj1 = ciat.RealInterval(obj1, obj1);
            end
            if isa(obj2, 'ciat.RealInterval') == 0
                obj2 = ciat.RealInterval(obj2, obj2);
            end
            r = ciat.RealInterval(max(obj1.Infimum, obj2.Infimum), ...
                max(obj1.Supremum, obj2.Supremum));
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
        r = power(obj1,obj2)
        
    end

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Infimum = obj.Infimum.(indexOp(1));
            obj.Supremum = obj.Supremum.(indexOp(1));
            if isscalar(indexOp)
                varargout{1} = obj;
                return;
            end
            [varargout{1:nargout}] = obj.(indexOp(2:end));
        end

        function obj = parenAssign(obj,indexOp,varargin)
            % POTENTIAL UNPEXPECTED BEHAVIOUR HERE
            % Only works for 2D arrays, not all is tested
            % Probably not all cases are covered

            % Warning, does not work for operations like
            % obj(1,1).Supremum = 1;
            % Should use
            % obj.Supremum(1,1) = 1;

            % Ensure object instance is the first argument of call.
            if isempty(obj)
                % Object must be of the correct size
                
                % If rhs is of size 1, then use indices to set size
                if isscalar(varargin{1})
                    sz = [indexOp.Indices{:}];
                    obj = ciat.RealInterval;
                    obj.Infimum = zeros(sz);
                    obj.Supremum = zeros(sz);
                else
                    obj = varargin{1};
                end
            end
            if isempty(varargin{1})
                % When rhs is empty, allocate memory for the object, size of indexOp.

                % Size to allocate
                tmp = indexOp.Indices;
                % Replace ':' with 1, to avoid errors when indexing into empty arrays.
                tmp(strcmp(':', indexOp.Indices)) = {1};
                sz = max(cellfun(@numel, tmp), cellfun(@max, tmp));
                
                obj = ciat.RealInterval;
                obj.Infimum = zeros(sz);
                obj.Supremum = zeros(sz);
                return;
            end
            if numel(indexOp) == 1
                if isscalar(indexOp(1))
                    assert(nargin==3);
                    rhs = varargin{1};
                    % If rhs is not an interval, then convert it to one.
                    if ~isa(rhs, 'ciat.RealInterval')
                        rhs = ciat.RealInterval(rhs);
                    end
                    obj.Infimum.(indexOp(1)) = rhs.Infimum;
                    obj.Supremum.(indexOp(1)) = rhs.Supremum;
                    return;
                end
            end
            tmp = obj.(indexOp(1));
            [tmp.(indexOp(2:end))] = varargin{:};
            obj.(indexOp(1)) = tmp;
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
            obj.Infimum.(indexOp) = [];
            obj.Supremum.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.RealInterval')
                    newArgs{ix} = varargin{ix}.Infimum;
                    newArgs2{ix} = varargin{ix}.Supremum;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ciat.RealInterval(cat(dim,newArgs{:}), cat(dim,newArgs2{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Infimum,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            disp('empty')
            obj = ciat.RealInterval;
        end
    end

    methods
        function obj = reshape(obj,varargin)
            obj.Infimum = reshape(obj.Infimum,varargin{:});
            obj.Supremum = reshape(obj.Supremum,varargin{:});
        end
    end
end

