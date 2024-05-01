classdef Arc < matlab.mixin.indexing.RedefinesParen

	properties
        Center          % Arc center as a complex number
        Radius          % Arc radius as a real number (negative means concave arc, zero means vertex)
	end	

	properties (Dependent)
    	ArcAngles			% Arc angle as a real interval (always counter-clockwise)
		StartPoint		% First point of the arc as complex numbers
        MidPoint        % MidPoint of the arc as complex numbers
        StopPoint		% Last point of the arc as complex numbers
		Length			% Length of the arc curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 centered arc (zero is zero centered)
		CurveParameter	% Parameter interval of the normalized curve as real interval
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
    end

    properties (Access = private)
       ArcAngleStore        % Storage property for the arc angle
    end

	methods
		%% Constructor
        function obj = Arc(varargin)
            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 3
                    center = varargin{1};
                    radius = varargin{2};
                    angles = varargin{3};
                    mustBeA(center,'double')
                    mustBeA(radius,'double')
                    mustBeA(angles,{'double','ciat.RealInterval'})
                    assert(all([size(center) == size(radius),...
                                size(center) == size(angles)]))
                    obj.Center = center;
                    obj.Radius = radius;
                    if strcmp(class(angles),'double')
                        obj.ArcAngles = ciat.RealInterval(angles);
                    else
                        obj.ArcAngles = angles;
                    end
                otherwise
                    error('incorrect number of input')
            end
        end

        %% Defining properties
                   
        % Set angles (store in the hidden property ArcAngleStore after wrapping)
        function obj = set.ArcAngles(obj,angleArray)

            [M,N] = size(obj);
            [Ma,Na] = size(angleArray);
            assert(M==Ma && N==Na)

            % Initialize ArcAngles property if necessary
            if isempty(obj.ArcAngles)
                obj.ArcAngleStore = repmat(ciat.RealInterval(0),M,N);
            end

            for m = 1:M
            for n = 1:N

                % Select angle interval
                angles = angleArray(m,n);

                % Check the width of the interval
                widthPi = floor(angles.Width / pi);
    
                % Check if the (-pi,pi) wrapped order is the same
                wrapInf = wrapToPi(angles.Infimum);
                wrapSup = wrapToPi(angles.Supremum);
                wrapOrder = ( wrapInf < wrapSup );
    
                switch widthPi
                    case 0  % Less than a Pi width
                        if wrapOrder    % Normal order
                            % Wrap both angles to (-pi,pi)
                            angInf = wrapInf;
                            angSup = wrapSup;
                        else            % wrapInf is pos and wrapSup is neg
                            % Wrap the infimum to (-2pi,-pi),
                            % wrap the supremum to (-pi,0)
                            angInf = wrapInf - 2*pi;
                            angSup = wrapSup;
                        end
    
                    case 1 % More than a Pi width but less than 2Pi
                        if wrapOrder    % Normal order
                            % Wrap both angles to (-pi,pi)
                            angInf = wrapInf;
                            angSup = wrapSup;
                        else            % wrapped angles have the same sign 
                                        % and wrapSup is lower than wrapInf
                            % Wrap the infimum to (-2pi,0),
                            % wrap the supremum to (0,2*pi)
                            angInf = wrapTo2Pi(angles.Infimum) - 2*pi;
                            angSup = wrapToPi(angles.Supremum);
                        end
                    otherwise           % full circle
                        angInf = -pi;
                        angSup = pi;
                end
                % Assign angle values
                obj.ArcAngleStore(m,n) = ciat.RealInterval(angInf,angSup);
            end
            end
        end

        % Get points (retrieve ArcAngles from hidden property ArcAngleStore)
        function value = get.ArcAngles(obj)
            value = obj.ArcAngleStore; 
        end

		%% Dependent properties
        

        % StartPoint
        function value = get.StartPoint(obj)
            angMask = (obj.Radius < 0);
            angInf = obj.ArcAngles.Infimum;
            angSup = obj.ArcAngles.Supremum;
            angStart = angInf.*angMask + angSup.*~angMask;
            value = obj.Center + obj.Radius .* exp(1i*angStart);
        end

        % MidPoint
        function value = get.MidPoint(obj)
            value = obj.Center + obj.Radius * exp(1i*obj.ArcAngles.MidPoint);
        end

        % StopPoint
        function value = get.StopPoint(obj)
            angMask = (obj.Radius < 0);
            angInf = obj.ArcAngles.Infimum;
            angSup = obj.ArcAngles.Supremum;
            angStop = angInf.*~angMask + angSup.*angMask;
            value = obj.Center + obj.Radius .* exp(1i*angStop);
        end

        % Length
        function value = get.Length(obj)
            value = 2*pi*obj.Radius * obj.ArcAngles.Width/(2*pi);
        end

        % Gauss map angle interval
        function value = get.GaussMap(obj)
            value = obj.ArcAngles;
            if ~isempty(value) && value.isin(pi) && all(pi~=value.Bounds)
                value = [ciat.RealInterval(wrapToPi(value.Infimum),pi);...
                         ciat.RealInterval(-pi,wrapToPi(value.Supremum))];
            end
        end

        % Log-Gauss map angle interval
        function value = get.LogGaussMap(obj)
            value = obj.GaussMap - ciat.RealInterval(angle(obj.StartPoint),...
                                                     angle(obj.StopPoint));
        end

        % Normalization factor
        function value = get.NormFactor(obj)
            if obj.Center == 0
                value = 0;
            else
                value = 1/obj.Center;
            end
        end

        % Curve parameter
        function value = get.CurveParameter(obj)
            value = obj.ArcAngles + angle(obj.NormFactor);
        end

        % Real
        function value = get.Real(obj)
            
            % Create boolean masks
            isConcave = obj.Radius < 0;
            crossInf = obj.ArcAngles.isin((isConcave-3)*pi) | ...
                       obj.ArcAngles.isin((isConcave-1)*pi) | ...
                       obj.ArcAngles.isin((isConcave+1)*pi) ;
            crossSup = obj.ArcAngles.isin((isConcave-2)*pi) | ...
                       obj.ArcAngles.isin((isConcave+0)*pi);
            
            % Calculate real bounds of the endpoints
            pntInf = min(real(obj.StartPoint),real(obj.StopPoint));
            pntSup = max(real(obj.StartPoint),real(obj.StopPoint));

            % Calculate real bounds of the envelope
            envInf = real(obj.Center) - abs(obj.Radius);
            envSup = real(obj.Center) + abs(obj.Radius);

            % Pick the correct value according to the masks
            realInf = (pntInf .* ~crossInf) + (envInf .* crossInf);
            realSup = (pntSup .* ~crossSup) + (envSup .* crossSup);

            value = ciat.RealInterval(realInf,realSup);

        end
        function value = real(obj)
            value = obj.Real;
        end
        
        % Imag
        function value = get.Imag(obj)
            
            % Create boolean masks
            isConcave = obj.Radius < 0;
            crossInf = obj.ArcAngles.isin((isConcave-2.5)*pi) | ...
                       obj.ArcAngles.isin((isConcave-0.5)*pi) | ...
                       obj.ArcAngles.isin((isConcave+1.5)*pi) ;
            crossSup = obj.ArcAngles.isin((isConcave-1.5)*pi) | ...
                       obj.ArcAngles.isin((isConcave+0.5)*pi);
            
            % Calculate real bounds of the endpoints
            pntInf = min(imag(obj.StartPoint),imag(obj.StopPoint));
            pntSup = max(imag(obj.StartPoint),imag(obj.StopPoint));

            % Calculate imag bounds of the envelope
            envInf = imag(obj.Center) - abs(obj.Radius);
            envSup = imag(obj.Center) + abs(obj.Radius);

            % Pick the correct value according to the masks
            imagInf = (pntInf .* ~crossInf) + (envInf .* crossInf);
            imagSup = (pntSup .* ~crossSup) + (envSup .* crossSup);

            value = ciat.RealInterval(imagInf,imagSup);

        end
        function value = imag(obj)
            value = obj.Imag;
        end
        
        % Abs
        function value = get.Abs(obj)

            % Create boolean masks
            isConcave = obj.Radius < 0;
            angOffset = angle(obj.Center);
            crossInf = obj.ArcAngles.isin(angOffset+(isConcave-3)*pi) | ...
                       obj.ArcAngles.isin(angOffset+(isConcave-1)*pi) | ...
                       obj.ArcAngles.isin(angOffset+(isConcave+1)*pi);
            crossSup = obj.ArcAngles.isin(angOffset+(isConcave-2)*pi) | ...
                       obj.ArcAngles.isin(angOffset+(isConcave+0)*pi);
            
            % Calculate real bounds of the endpoints
            pntInf = min(abs(obj.StartPoint),abs(obj.StopPoint));
            pntSup = max(abs(obj.StartPoint),abs(obj.StopPoint));

            % Calculate abs bounds of the envelope
            envInf = abs(obj.Center) - abs(obj.Radius);
            envSup = abs(obj.Center) + abs(obj.Radius);

            % Pick the correct value according to the masks
            absInf = (pntInf .* ~crossInf) + (envInf .* crossInf);
            absSup = (pntSup .* ~crossSup) + (envSup .* crossSup);

            value = ciat.RealInterval(absInf,absSup);
        end
        function value = abs(obj)
            value = obj.Abs;
        end
        
        % Angle
        function value = get.Angle(obj)

            % Calculate real bounds of the endpoints
            pntInf = min(angle(obj.StartPoint),angle(obj.StopPoint));
            pntSup = max(angle(obj.StartPoint),angle(obj.StopPoint));

            % Calculate angle bounds of the envelope
            envInf = angle(obj.Center) - asin(abs(obj.Radius./obj.Center));
            envSup = angle(obj.Center) + asin(abs(obj.Radius./obj.Center));

            % Create boolean masks
            isConcave = obj.Radius < 0;
            crossInf = obj.ArcAngles.isin(envInf+(isConcave-2.5)*pi) | ...
                       obj.ArcAngles.isin(envInf+(isConcave-0.5)*pi) | ...
                       obj.ArcAngles.isin(envInf+(isConcave+1.5)*pi);
            crossSup = obj.ArcAngles.isin(envSup+(isConcave-3.5)*pi) | ...
                       obj.ArcAngles.isin(envSup+(isConcave-1.5)*pi) | ...
                       obj.ArcAngles.isin(envSup+(isConcave+0.5)*pi) | ...
                       obj.ArcAngles.isin(envSup+(isConcave+2.5)*pi);

            % Pick the correct value according to the masks
            angleInf = (pntInf .* ~crossInf) + (envInf .* crossInf);
            angleSup = (pntSup .* ~crossSup) + (envSup .* crossSup);

            value = ciat.RealInterval(angleInf,angleSup);
        end
        function value = angle(obj)
            value = obj.Angle;
        end

		%% Other methods

        % Plot
        function h = plot(obj, varargin)
        % Plot arcs
        %
        % This function plots a set of polygonal intervals 
        % (see MATLAB plot function)
        % _________________________________________________________________________
        % USAGE        
        %   r = plot(obj)
        % _________________________________________________________________________
        % NECESSARY ARGUMENT
        %   obj       : array of objects from the ciat.Arc class
        % _________________________________________________________________________
        % OPTIONS
        % _________________________________________________________________________
        % EXAMPLES
        %   
        % _________________________________________________________________________
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                % Extract parameters
                cr = real(obj(n).Center);
                ci = imag(obj(n).Center);
                r = obj(n).Radius;
                % Plot arc by angle quadrants between in (-2Pi,2Pi)
                % q = [-pi,-pi/2 ; -pi/2,0 ; 0,pi/2 ; pi/2,pi];
                q = ((-2:0.5:1.5)'+[0,0.5])*pi;
                if r~=0
                    for idx = 1:length(q)
                        qi = intersection(obj(n).ArcAngles, ...
                                        ciat.RealInterval(q(idx,1),q(idx,2)));
                        if ~isempty(qi)
                            xBound = sort(cr + r*cos(qi.Bounds) );
                            yBound = sort(ci + r*sin(qi.Bounds) );
                            h = [h; ...
                                fimplicit(@(x,y) (x-cr).^2 + (y-ci).^2 - r^2 , ...
                                            [xBound yBound],varargin{:})];    
                        end
                    end
                else
                    h = [h; plot(cr,ci,varargin{:},'Marker','.')];
                end
            end
            if tf == false 
                hold off
            end
        end

        % Plot Gauss maps
        function h = plotMap(obj, logMap, arrowSize, varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];

            % Plot normal vectors
            for n = 1:length(obj(:))
                % Extract variables
                arc = obj(n);
                if logMap == 0
                    map = arc.GaussMap;
                else
                    map = arc.LogGaussMap;
                end

                % Set arrow positions
                p = zeros(1,3);
                p(1) = arc.StartPoint;
                p(2) = arc.MidPoint;
                p(3) = arc.StopPoint;

                % Set vector lengths
                if length(map) == 1
                    a = zeros(1,3);
                else
                    a = zeros(1,6);
                    p = reshape(repmat(p,1,2),1,[]);
                end

                % Set arrow angle and length
                for m = 1:length(map)
                    a(1+3*(m-1)) = arrowSize * exp(1i*map(m).Infimum);
                    a(2+3*(m-1)) = arrowSize * exp(1i*map(m).MidPoint);
                    a(3+3*(m-1)) = arrowSize * exp(1i*map(m).Supremum);
                end

                % Flip vectors if arc is concave 
                if obj(n).Radius < 0
                    a = flip(a);
                end

                % Plot arrows
                h=[h; ...
                    quiver(real(p),imag(p),real(a),imag(a),...
                           varargin{:},'AutoScale','off')];
            end

            if tf == false 
                hold off
            end
        end
        function h = plotGaussMap(obj, arrowSize, varargin)
            h = obj.plotMap(0,arrowSize,varargin{:});
        end
        function h = plotLogGaussMap(obj, arrowSize, varargin)
            h = obj.plotMap(1,arrowSize,varargin{:});
        end
        
        %% Function headers
        r = plus(obj1,obj2)

	end


%% Vectorizing the object

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Center = obj.Center.(indexOp(1));
            obj.Radius = obj.Radius.(indexOp(1));
            obj.ArcAngleStore = obj.ArcAngleStore.(indexOp(1));
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
                    obj.Center = zeros(sz);
                    obj.Radius = zeros(sz);
                    obj.ArcAngleStore = repmat(ciat.RealInterval,sz);
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
                
                obj = ciat.Arc;
                obj.Center = zeros(sz);
                obj.Radius = zeros(sz);
                obj.ArcAngleStore = repmat(ciat.RealInterval,sz);
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
                    obj.Center.(indexOp(1)) = rhs.Center;
                    obj.Radius.(indexOp(1)) = rhs.Radius;
                    obj.ArcAngleStore.(indexOp(1)) = rhs.ArcAngleStore;
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
                    newArgs{ix} = varargin{ix}.Center;
                    newArgs2{ix} = varargin{ix}.Radius;
                    newArgs3{ix} = varargin{ix}.ArcAngleStore;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ciat.Arc(cat(dim,newArgs{:}), ...
                           cat(dim,newArgs2{:}), ...
                           cat(dim,newArgs3{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Center,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            %disp('empty')
            obj = ciat.Arc;
        end
    end

    methods
        function obj = reshape(obj,varargin)
            obj.Center = reshape(obj.Center,varargin{:});
            obj.Radius = reshape(obj.Radius,varargin{:});
            obj.ArcAngleStore = reshape(obj.ArcAngleStore,varargin{:});
        end
    end
end