classdef PolarInterval < matlab.mixin.indexing.RedefinesParen

% Polar interval class for complex interval arithmetic calculations
%
% This is a class of the Complex Interval Arithmetic Toolbox.
% It allows the definition of complex intervals represented by polar
% regions (annular sectors) in the complex plane defined by a absolute 
% value and angle interval (real intervals). The object 
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
        Abs;    % Projection of the polar interval to the absolute value axis
        Angle;  % Projection of the polar interval to the angle axis
    end
    
    properties (Dependent)
        Real;    % Projection of the polar interval to the real axis
        Imag;   % Projection of the polar interval to the imaginary axis
        Area;   % Area of the polar interval
    end
    
    methods
        %% Constructor
        function obj = PolarInterval(varargin)
        %POLARINTERVAL Construct an instance of this class
        %
        % This function generates one or more polar intervals
        % based on the optional input arguments, each having a 
        % default value. If no argument is given, the generated
        % object has empty properties, which is useful for 
        % initialization of an array of intervals.
        % There are multiple ways of defining an interval. The default
        % method is givin an absolute value and angle interval as two 
        % seperate arguments of real interval type (it can be an array or 
        % matrix). 
        % Another default method is giving the abs infimum, abs supremum,
        % angle infimum and angle supremum values as four seperate
        % arguments of double type (can be an array or matrix).
        % It is also possible to give only a single input argument,
        % in which case it will be converted to polar intervals, 
        % representing the smallest enclosing interval. In case of 
        % double type input this results in degenerate intervals.
        %__________________________________________________________________________
        % USAGE        
        %   ciat.PolarInterval(absMin,absMax,angleMin,angleMax)
        %   ciat.PolarInterval(obj)
        %   ciat.PolarInterval
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        % _________________________________________________________________________
        % OPTIONS
        %   abs     : absolute value as real interval type array
        %   angle   : angle value as real interval type array
        %   absInf  : abs infimum value as double type array
        %   absSup  : abs supremum value as double type array
        %   angleInf: angle infimum value as double type array
        %   angleSup: angle supremum value as double type array
        %   obj       : object of double or another complex interval type
        %               (see cast function for details)
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = ciat.PolarInterval(1,2,3,4)
        %   polarInt = ciat.PolarInterval([1,2],[2,3],[3,4],[4,5])
        %   polarInt = ciat.PolarInterval([1+1i,2+2i])
        %   polarInt = ciat.PolarInterval(ciat.PolarInterval(0,1,2,3))
        %   polarInt(5,2) = ciat.PolarInterval
        % _________________________________________________________________________
            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 1
                    % Single input argument is casted
                    [M,N] = size(varargin{1});
                    obj(M,N) = obj;
                    for n = 1:M*N
                        obj(n) = ciat.PolarInterval.cast(varargin{1}(n));
                    end
                case 2
                    % Two input arguments of real intervals is the default
                    % way of defining polar intervals
                    mustBeA(varargin{1},["ciat.RealInterval","double"])
                    mustBeA(varargin{2},["ciat.RealInterval","double"])

                    % Turn scalars to degenerate intervals
                    if isa(varargin{1}, 'double')
                        varargin{1} = ciat.RealInterval(varargin{1});
                    end
                    if isa(varargin{2}, 'double')
                        varargin{2} = ciat.RealInterval(varargin{2});
                    end 
                    
                    % Get sizes
                    [M1,N1] = size(varargin{1});
                    [M2,N2] = size(varargin{2});
                    assert(M1 == M2 && N1 == N2)
                    M = M1;
                    N = N1;
                    
                    % Create object array
                    obj(M,N) = obj;
                    obj.Abs = varargin{1};
                    obj.Angle = varargin{2};
                case 4
                    % Four input arguments of doubles is also a default way
                    % of defining polar intervals 
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
                    obj.Abs = ciat.RealInterval(varargin{1},varargin{2});
                    obj.Angle = ciat.RealInterval(varargin{3},varargin{4});
            end
        end
        
        %% Defining properties
        
        % Abs
        function value = abs(obj)
        % Absolute value of polar intervals
        %
        % This function creates the real interval representing the 
        % absolute value of a set of polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = abs(ciat.PolarInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Abs;
        end
        
        % Angle
        function value = angle(obj)
         % Angle of polar intervals
        %
        % This function creates the real interval representing the 
        % angle of a set of polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = abs(ciat.PolarInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Angle;
        end
        
        %% Dependent properties
        
        % Real
        function value = get.Real(obj)
            value = obj.Abs .* cos( obj.Angle );
        end
        function value = real(obj)
        % Real value of polar intervals
        %
        % This function creates the real interval representing the 
        % real values of a set of polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = real(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = real(ciat.PolarInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Real;
        end
        
        % Imag
        function value = get.Imag(obj)
            value = obj.Abs * sin( obj.Angle );
        end
        function value = imag(obj)
        % Imaginary value of polar intervals
        %
        % This function creates the real interval representing the 
        % imaginary values of a set of polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = imag(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = imag(ciat.PolarInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Imag;
        end
        
        % Area
        function value = get.Area(obj)
            value = 0.5 * width(obj.Angle) .* ( obj.Abs.Supremum.^2 - ...
                                                obj.Abs.Infimum.^2 );
        end 
        function value = area(obj)
        % Area of polar intervals
        %
        % This function returns the area valuea of a set of 
        % polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = abs(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = abs(ciat.PolarInterval(0,1,2,3));
        % _________________________________________________________________________
            value = obj.Area;
        end
       
        %% Other methods

        function r = eq(obj1,obj2)
            % Equality of polar intervals
            %
            % This function checks if two polar intervals are equal
            % _________________________________________________________________________
            % USAGE
            %   r = eq(obj1,obj2)
            % _________________________________________________________________________
            % NECESSARY ARGUMENTS
            %   obj1      : array of objects from the ciat.PolarInterval class
            %   obj2      : array of objects from the ciat.PolarInterval class
            % _________________________________________________________________________
            % OPTIONS
            % _________________________________________________________________________
            % EXAMPLES
            %   r = eq(ciat.PolarInterval(0,1,2,3),ciat.PolarInterval(0,1,2,3));
            % _________________________________________________________________________
                    [M1,N1] = size(obj1);
                    [M2,N2] = size(obj2);
                    assert(M1 == M2 && N1 == N2)
                    r = all(obj1.Abs == obj2.Abs, "all") && ...
                        all(obj1.Angle == obj2.Angle, "all");
        end

        function r = ne(obj1,obj2)
            % Inequality of polar intervals
            %
            % This function checks if two polar intervals are not equal
            % _________________________________________________________________________
            % USAGE
            %   r = ne(obj1,obj2)
            % _________________________________________________________________________
            % NECESSARY ARGUMENTS
            %   obj1      : array of objects from the ciat.PolarInterval class
            %   obj2      : array of objects from the ciat.PolarInterval class
            % _________________________________________________________________________
            % OPTIONS
            % _________________________________________________________________________
            % EXAMPLES
            %   r = ne(ciat.PolarInterval(0,1,2,3),ciat.PolarInterval(0,1,2,3));
            % _________________________________________________________________________
            [M1,N1] = size(obj1);
            [M2,N2] = size(obj2);
            assert(M1 == M2 && N1 == N2)
            r = any(obj1.Abs ~= obj2.Abs, "all") || ...
                any(obj1.Angle ~= obj2.Angle, "all");
        end

        function r = transpose(obj)
            r = obj;
            r.Abs = r.Abs.';
            r.Angle = r.Angle.';
        end

        function r = ctranspose(obj)
            r = obj;
            r.Abs = r.Abs.';
            r.Angle = -r.Angle.';
        end

        % In interval
        function r = ininterval(obj, points)
        % Check if points are in polar intervals
        %
        % This function checks if a set of points are in a set of
        % polar intervals
        % _________________________________________________________________________
        % USAGE
        %   r = ininterval(obj, points)
        % _________________________________________________________________________
        % NECESSARY ARGUMENTS
        %   obj       : array of objects from the ciat.PolarInterval class
        %   points    : array of complex numbers
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   r = ininterval(ciat.PolarInterval(0,1,2,3),[1,2,3,4,5]);
        % _________________________________________________________________________
            % r = obj.Abs.ininterval(abs(points)) & ...
            %     obj.Angle.ininterval(angle(points));
            % wrap angle interval to [-pi,pi]
            r = obj.Abs.ininterval(abs(points)) & ...
                sin(obj.Angle).ininterval(sin(angle(points))) & ...
                cos(obj.Angle).ininterval(cos(angle(points)));
        end

        % Union
        function r = union(obj)
        % Union of polar intervals
        %
        % This function creates the polar interval representing the 
        % union of a set of polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = union(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = union([ciat.PolarInterval(0,1,2,3), ...
        %                    ciat.PolarInterval(2,3,4,5)]);
        % _________________________________________________________________________
            N = length(obj(:));
            assert(N>1)
            r = ciat.PolarInterval(union([obj.Abs]) , union([obj.Angle]));
        end
        
        % Intersection
        function r = intersection(obj)
        % Intersection of polar intervals
        %
        % This function creates the polar interval representing the 
        % intersection of a set of polar intervals
        % _________________________________________________________________________
        % USAGE        
        %   r = intersection(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = intersection([ciat.PolarInterval(0,1,2,3), ...
        %                    ciat.PolarInterval(2,3,4,5)]);
        % _________________________________________________________________________
            N = length(obj(:));
            assert(N>1)
            r = ciat.PolarInterval(intersection([obj.Abs]) , ...
                                   intersection([obj.Angle]));
        end
        
        % Negative (uminus)
        function r = uminus(obj)
        % Negative of polar intervals (- operator)
        %
        % This function creates the polar interval representing the 
        % negative of a set of polar intervals (see MATLAB uminus
        % function)
        % _________________________________________________________________________
        % USAGE        
        %   r = -obj
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   polarInt = -ciat.PolarInterval(0,1,2,3);
        % _________________________________________________________________________
            r = obj;
            for n = 1:length(r(:))
                r(n).Angle = r(n).Angle + pi;
            end
        end    

        % IsNaN
        function r = isnan(obj)
            r = isnan(obj.Area);
        end 

        % Inside
        function r = isin(obj,x)
            % Check if the point is in between the circles
            inCircle = obj.Abs.isin( abs(x) );

            % Check if the point is in the sector
            % inSector = abs(wrapToPi( angle(x) - obj.Angle.Infimum ) ) ...
                                    % <= obj.Angle.Width;
            xAng = angle(x);
            if any(obj.Angle.isin([-pi,pi]))
                inSector = (xAng>=0 & xAng >= wrapToPi(obj.Angle.Infimum))| ...
                           (xAng <0 & xAng <= wrapToPi(obj.Angle.Supremum));
            else
                inSector = xAng >= obj.Angle.Infimum & ...
                           xAng <= obj.Angle.Supremum;
            end
            

            % Check if the point is in the rectangle around the circle
            arcBox = ciat.RectangularInterval(obj);
            inBox = arcBox.isin(x);

            % Combine conditions
            r = inCircle & inSector & inBox;
        end
        
        % Plot
        function h = plot(obj, varargin)
        % Plot polar intervals 
        %
        % This function plots a set of polar intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.PolarInterval class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   h = plot(ciat.PolarInterval(0,1,2,3));
        % _________________________________________________________________________
            tf = ishold; 
            if tf == false 
                clf
            end
            hold on
            [M,N] = size(obj); 
            h = [];
            for n = 1:length(obj(:))
                angles = linspace(obj(n).Angle.Infimum, ...
                                  obj(n).Angle.Supremum, 200);
                curve_outer = obj(n).Abs.Supremum * exp(1j*angles);
                curve_inner = obj(n).Abs.Infimum * exp(1j*angles);
                startPoint = obj(n).Abs.Supremum * exp(1j * ... 
                                                obj(n).Angle.Infimum);
                vec = [ curve_outer, flip(curve_inner), startPoint];
                h = [h;plot(real(vec), imag(vec), varargin{:})];
            end
            h = reshape(h,M,N);
            if tf == false 
                hold off
            end
        end  
        
        %% Function headers
        points = sample(obj, n_points)
        r = times(obj1,obj2)
        r = mtimes(obj1,obj2)
        r = exp(obj)
    end
    
    methods (Static)
       outObj = cast(inObj)
    end

    %% Vectorization
    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Abs = obj.Abs.(indexOp(1));
            obj.Angle = obj.Angle.(indexOp(1));
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
                % such as doing obj(5,2) = ciat.PolarInterval
                % Might not be the place or the way to do it

                % Instanciate object with zero values of correct size.
                obj = ciat.PolarInterval;
                obj.Abs = ciat.RealInterval(nan([indexOp.Indices{:}]), ...
                                            nan([indexOp.Indices{:}]));
                obj.Angle = ciat.RealInterval(nan([indexOp.Indices{:}]), ...
                                              nan([indexOp.Indices{:}]));

                % obj = varargin{1};
                varargin{1} = obj.(indexOp);
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.Abs.(indexOp) = rhs.Abs;
                obj.Angle.(indexOp) = rhs.Angle;
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
            obj.Abs.(indexOp) = [];
            obj.Angle.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.PolarInterval')
                    newArgs{ix} = varargin{ix}.Abs;
                    newArgs2{ix} = varargin{ix}.Angle;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ciat.PolarInterval(cat(dim,newArgs{:}), ...
                                     cat(dim,newArgs2{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Abs,varargin{:});
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
            obj.Abs = reshape(obj.Abs,varargin{:});
            obj.Angle = reshape(obj.Angle,varargin{:});
        end
    end
end

