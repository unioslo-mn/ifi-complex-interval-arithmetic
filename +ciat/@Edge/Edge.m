classdef Edge < matlab.mixin.indexing.RedefinesParen

	properties
        Start           % First point of the edge
        Stop            % Last point of the edge
    end	

	properties (Dependent)
        Midpoint        % Midpoint of the arc as complex numbers
        Vector          % Endpoint difference as a complex vector
        Slope           % Slope of the line the edge fits on
        ZeroCrossing    % Real axis crossing point of the line the edge fits on
		Length			% Length of the edge curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 crossing vertical edge (zero is zero crossing)
		CurveParameter	% Parameter interval of the normalized curve as real interval
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
	end

	methods
		%% Constructor
        function obj = Edge(varargin)
            for varIdx = 1:length(varargin)
                mustBeA(varargin{varIdx},'double')
            end

            switch length(varargin)
                case 0
                    % This is for initializing an array of objects
                case 1
                    obj.Start = varargin{1};
                    obj.Stop = obj.Start;
                case 2
                    obj.Start = varargin{1};
                    obj.Stop = varargin{2};
            end
		end

		%% Dependent properties

        % Midpoint
        function value = get.Midpoint(obj)
            value = (obj.Start+obj.Stop)/2;
        end

        % Vector
        function value = get.Vector(obj)
            value = obj.Stop - obj.Start;
        end

        % Slope
		function value = get.Slope(obj)
            vd = obj.Vector;
            value = imag(vd)./real(vd);
        end

        % Zero crossing
		function value = get.ZeroCrossing(obj)
            value = real(obj.Start) - imag(obj.Start) ./ obj.Slope;
        end

        % Length
		function value = get.Length(obj)
            value = abs(obj.Vector);
        end

        % GaussMap
		function value = get.GaussMap(obj)
            % value = ciat.RealInterval(ciat.wrapToPi(angle(obj.Vector)-pi/2));
            value = ciat.wrapToPi(angle(obj.Vector)-pi/2);
        end

        % Log-GaussMap
		function value = get.LogGaussMap(obj)
            value = obj.GaussMap - ciat.RealInterval(angle(obj.Start), ...
                                                     angle(obj.Stop));
        end

        % Normalization factor
		function value = get.NormFactor(obj)
            if obj.ZeroCrossing == 0
                value = 0;
            else
                % Extract line parameters
                v1 = obj.Start;
                va = obj.Slope;        
                
                % Find scale-rotate factor
                vr = exp(-1i*(atan(va)+pi/2));   % rotation factor
                vs = 1 ./ real(vr .* v1);            % scale factor
                value = vr .* vs;                 % complex factor

            end
        end

        function value = get.CurveParameter(obj)
            value = ciat.RealInterval(imag(obj.Start .* obj.NormFactor),...
                                        imag(obj.Stop .* obj.NormFactor));
        end

        % Real
        function value = get.Real(obj)
            value = ciat.RealInterval(min(real(obj.Start),real(obj.Stop)),...
                                      max(real(obj.Start),real(obj.Stop)));
        end
        function value = real(obj)
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
        end
        
        % Imag
        function value = get.Imag(obj)
            value = ciat.RealInterval(min(imag(obj.Start),imag(obj.Stop)),...
                                      max(imag(obj.Start),imag(obj.Stop)));
        end
        function value = imag(obj)
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            if obj.CurveParameter.isin(0)
                minAbs = 1/abs(obj.NormFactor);
            else
                minAbs = min(abs(obj.Start),abs(obj.Stop));
            end
            maxAbs = max(abs(obj.Start),abs(obj.Stop));
            value = ciat.RealInterval(minAbs,maxAbs);
        end
        function value = abs(obj)
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            value = ciat.RealInterval(min(angle(obj.Start),min(angle(obj.Stop))),...
                                      max(angle(obj.Start),max(angle(obj.Stop))));
        end
        function value = angle(obj)
            [M,N] = size(obj);
            value = reshape([obj.Angle],M,N);
        end

        %% Other methods

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
                points = [obj(n).Start , obj(n).Stop];
                if obj(n).Length ~= 0
                    h = [h;plot(real(points), imag(points), varargin{:})];    
                else
                    h = [h;plot(real(points), imag(points), varargin{:},...
                        'Marker','.')];    
                end
            end
            if tf == false 
                hold off
            end
        end

        % Plot
        function h = plotMap(obj,logMap,arrowSize,varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];

            % Plot normal vectors
            for n = 1:length(obj(:))
                % Extact variables
                edge = obj(n);
                if logMap == 0
                    map = edge.GaussMap;
                else
                    map = edge.LogGaussMap;
                end

                % Set arrow positions
                p(1) = edge.Start;
                p(2) = edge.Midpoint;
                p(3) = edge.Stop;

                % Set vector lengths
                if length(map) == 1
                    a = zeros(1,3);
                else
                    a = zeros(1,6);
                    p = repmat(p,1,2);
                end

                % Set arrow angle and length
                for m = 1:length(map)
                    a(1+3*(m-1)) = arrowSize * exp(1i*map(m).Infimum);
                    a(2+3*(m-1)) = arrowSize * exp(1i*map(m).Midpoint);
                    a(3+3*(m-1)) = arrowSize * exp(1i*map(m).Supremum);
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

    end

    %% Vectorizing the object

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Start = obj.Start.(indexOp(1));
            obj.Stop = obj.Stop.(indexOp(1));
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
            % obj(1,1).Stop = 1;
            % Should use
            % obj.Stop(1,1) = 1;

            % Ensure object instance is the first argument of call.
            if isempty(obj)
                % Object must be of the correct size
                
                % If rhs is of size 1, then use indices to set size
                if isscalar(varargin{1})
                    sz = [indexOp.Indices{:}];
                    obj = ciat.Edge;
                    obj.Start = zeros(sz);
                    obj.Stop = zeros(sz);
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
                
                obj = ciat.Edge;
                obj.Start = zeros(sz);
                obj.Stop = zeros(sz);
                return;
            end
            if numel(indexOp) == 1
                if isscalar(indexOp(1))
                    assert(nargin==3);
                    rhs = varargin{1};
                    % If rhs is not an interval, then convert it to one.
                    if ~isa(rhs, 'ciat.Edge')
                        rhs = ciat.Edge(rhs);
                    end
                    obj.Start.(indexOp(1)) = rhs.Start;
                    obj.Stop.(indexOp(1)) = rhs.Stop;
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
            obj.Start.(indexOp) = [];
            obj.Stop.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.Edge')
                    newArgs{ix} = varargin{ix}.Start;
                    newArgs2{ix} = varargin{ix}.Stop;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            out = ciat.Edge(cat(dim,newArgs{:}), cat(dim,newArgs2{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Start,varargin{:});
        end
    end

    methods (Static, Access=public)
        function obj = empty()
            %disp('empty')
            obj = ciat.Edge;
        end
    end

    methods
        function obj = reshape(obj,varargin)
            obj.Start = reshape(obj.Start,varargin{:});
            obj.Stop = reshape(obj.Stop,varargin{:});
        end
    end

end