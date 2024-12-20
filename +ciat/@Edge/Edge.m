classdef Edge < matlab.mixin.indexing.RedefinesParen

	properties
        Startpoint           % First point of the edge
        Endpoint            % Last point of the edge
    end	

	properties (Dependent)
        Midpoint        % Midpoint of the arc as complex numbers
        Vector          % Endpoint difference as a complex vector
        Slope           % Slope of the line the edge fits on
        Offset          % Offset of the line the edge fits on
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
                    if length(varargin{1})==1
                        obj.Startpoint = varargin{1};
                        obj.Endpoint = obj.Startpoint;
                    elseif length(varargin{1})==2
                        obj.Startpoint = varargin{1}(1);
                        obj.Endpoint = varargin{1}(2);
                    else
                        error('Input array is too long')
                    end
                case 2
                    obj.Startpoint = varargin{1};
                    obj.Endpoint = varargin{2};
            end
		end

		%% Dependent properties

        % Midpoint
        function value = get.Midpoint(obj)
            value = (obj.Startpoint+obj.Endpoint)/2;
        end

        % Vector
        function value = get.Vector(obj)
            value = obj.Endpoint - obj.Startpoint;
        end

        % Slope
		function value = get.Slope(obj)
            vd = obj.Vector;
            value = imag(vd)./real(vd);
        end

        % Offset (the imaginary value at zero real value)
        function value = get.Offset(obj)
            value = imag(obj.Startpoint) - real(obj.Startpoint) .* obj.Slope;
        end

        % Zero crossing
		function value = get.ZeroCrossing(obj)
            value = real(obj.Startpoint) - imag(obj.Startpoint) ./ obj.Slope;
        end

        % Length
		function value = get.Length(obj)
            value = abs(obj.Vector);
        end

        % GaussMap
		function value = get.GaussMap(obj)
            GaussAngle = ciat.wrapToPi(angle(obj.Vector)-pi/2);
            value = ciat.RealInterval(GaussAngle,GaussAngle);
            % value = ciat.wrapToPi(angle(obj.Vector)-pi/2);
        end

        % Log-GaussMap
		function value = get.LogGaussMap(obj)
            % Calculate LGM at start and endpoint
            startLGM = ciat.wrapToPi(obj.GaussMap.mid - angle(obj.Startpoint));
            endLGM = ciat.wrapToPi(obj.GaussMap.mid - angle(obj.Endpoint));

            % Extract curve parameter direction
            curveParDir = (obj.Endpoint - obj.Startpoint) .* obj.NormFactor;

            % Create LGM interval
            if imag(curveParDir) >= 0 
                value = ciat.RealInterval(endLGM - (endLGM>startLGM)*2*pi, ...
                                          startLGM);
            else
                value = ciat.RealInterval(startLGM - (startLGM>endLGM)*2*pi, ...
                                          endLGM);
            end
        end

        % Normalization factor
		function value = get.NormFactor(obj)
            [M,N] = size(obj);
            if obj.ZeroCrossing == 0
                value = ones(M,N);
            else
                % Extract line parameters
                v1 = obj.Startpoint;
                va = obj.Slope;        
                
                % Find scale-rotate factor
                vr = exp(-1i*(atan(va)+pi/2));   % rotation factor
                vs = 1 ./ real(vr .* v1);            % scale factor
                value = vr .* vs;                 % complex factor

            end
        end

        function value = get.CurveParameter(obj)
            value = ciat.RealInterval(imag(obj.Startpoint .* obj.NormFactor),...
                                        imag(obj.Endpoint .* obj.NormFactor));
        end

        % Real
        function value = get.Real(obj)
            value = ciat.RealInterval(min(real(obj.Startpoint), ...
                                          real(obj.Endpoint)),...
                                      max(real(obj.Startpoint), ...
                                          real(obj.Endpoint)));
        end
        function value = real(obj)
            [M,N] = size(obj);
            value = reshape([obj.Real],M,N);
        end
        
        % Imag
        function value = get.Imag(obj)
            value = ciat.RealInterval(min(imag(obj.Startpoint), ...
                                          imag(obj.Endpoint)),...
                                      max(imag(obj.Startpoint), ...
                                          imag(obj.Endpoint)));
        end
        function value = imag(obj)
            [M,N] = size(obj);
            value = reshape([obj.Imag],M,N);
        end
        
        % Abs
        function value = get.Abs(obj)
            curveParamZero = obj.CurveParameter.isin(0);
            minAbs =  curveParamZero ./ abs(obj.NormFactor) + ...
                     ~curveParamZero .* min(abs(obj.Startpoint), ...
                                            abs(obj.Endpoint));
            maxAbs = max(abs(obj.Startpoint),abs(obj.Endpoint));
            value = ciat.RealInterval(minAbs,maxAbs);
        end
        function value = abs(obj)
            [M,N] = size(obj);
            value = reshape([obj.Abs],M,N);
        end
        
        % Angle
        function value = get.Angle(obj)
            value = ciat.RealInterval( min(angle(cat(3,obj.Startpoint, ...
                                                    obj.Endpoint)),[],3),...
                                         max(angle(cat(3,obj.Startpoint,...
                                                    obj.Endpoint)),[],3) );
        end
        function value = angle(obj)
            value = obj.Angle;
        end

        % Negative
        function r = uminus(obj)
            r = ciat.Edge(-obj.Startpoint, -obj.Endpoint);
        end

        % Minus
        function r = minus(obj1,obj2)
            r = obj1 + (-obj2);
        end

        % Reciprocal
        function r = recip(obj) % Needs vectorization
            P1 = obj.Startpoint;
            P2 = obj.Endpoint;
            pEdge = mean([P1,P2]);
            if abs(obj.ZeroCrossing) < 100*eps
                if obj.ison(0)
                    warning('Edge contains zero, reciprocal returns NaN')
                    r = ciat.Edge;
                else
                    r = ciat.Edge(1./P1,1./P2);
                end
            else
                norm = obj.NormFactor;
                
                center = 0.5 * norm;
                radius = 0.5 * abs(norm) * ((abs(P2)-abs(P1)<0)-0.5)*2;
                ang = ciat.RealInterval(angle(1./P1-center),...
                                        angle(1./P2-center));
                ang = ang + pi * (~ang.isin(ciat.wrapToPi(...
                                    angle(1./pEdge-center) + pi*(radius < 0))));
    
                r = ciat.Arc(center, radius, ang);
            end
        end

        %% Other methods

        % Intersection
        function r = intersection(obj1,obj2)

            [M1,N1] = size(obj1);
            [M2,N2] = size(obj2);
            assert(M1 == M2 && N1 == N2)
            
            if isa(obj2,'ciat.Edge')
                % Find the intersection point of two lines
                a1 = obj1.Slope;
                b1 = obj1.Offset;
                a2 = obj2.Slope;
                b2 = obj2.Offset;
                if isinf(a1)
                    xCoord = real(obj1.Startpoint);
                    yCoord = a2 * xCoord + b2;
                elseif isinf(a2)
                    xCoord = real(obj2.Startpoint);
                    yCoord = a1 * xCoord + b1;
                else
                    xCoord = (b2 - b1) ./ (a1 - a2);
                    yCoord = a1 * xCoord + b1;
                end

                % Assign values
                r = nan(M1,N1);
                mask = obj1.Real.isin(xCoord) & obj1.Imag.isin(yCoord) & ...
                       obj2.Real.isin(xCoord) & obj2.Imag.isin(yCoord);
                if any(mask,'all')
                    r(mask) = xCoord(mask) + 1i * yCoord(mask);
                end
            elseif isa(obj2,'ciat.Arc')
                r = obj2.intersection(obj1);
            else
                error('Invalid input type')
            end
        end
        function r = cap(obj1,obj2)
            r = intersection(obj1,obj2);
        end

        % IsNaN
        function r = isnan(obj)
            r = isnan(obj.Length);
        end 

        % Is point on edge
        function r = ison(obj,x)
            % Check if the point is on the line of the edge
            onLine = abs(real(x) * obj.Slope + obj.Offset - imag(x)) < 100*eps;

            % Check if the point is on the edge real interval
            inReal = obj.Real.isin(real(x));

            % Combine conditions
            r = onLine & inReal;
        end

        % Split
        function edgeOut = split(edgeIn,splitPoint)

            % It only works for a single arc with multiple split points
            assert( numel(edgeIn)==1 )

            % Exclude invalid points
            splitPoint = splitPoint(~isnan(splitPoint));
            splitPoint = splitPoint(edgeIn.ison(splitPoint));
            splitPoint = splitPoint(splitPoint ~= edgeIn.Startpoint & ...
                                    splitPoint ~= edgeIn.Endpoint );

            % Split edge
            if ~isempty(splitPoint)
                % Extract edge parameters
                p1 = edgeIn.Startpoint;
                p2 = edgeIn.Endpoint;

                % Sort split points
                [~,sortIdx] = sort(abs(splitPoint-p1));
                splitPoint = [p1 ; splitPoint(sortIdx); p2];
                [~,uniqueIdx,~] = uniquetol(abs(splitPoint-p1),10*eps);
                splitPoint = splitPoint(uniqueIdx);

                % Calculate edge segments
                L = length(splitPoint)-1;
                edgeOut(L,1) = ciat.Edge;
                for l = 1:L
                    edgeOut(l) = ciat.Edge(splitPoint(l),splitPoint(l+1));
                end
            else
                edgeOut = edgeIn;
            end
        end

        % Find point with a given log-Gauss map
        function point = findLGM(obj,LGM)
            if obj.LogGaussMap.isin(LGM)
                % gFunc = @(s) -pi*(obj.LogGaussMap.mid<0) - angle(1+1i*s);
                gFunc = @(s) - angle(1+1i*s);
                sSolv=zeros(1,3);
                for idx = 1:3
                    sSolv(idx) = fsolve(@(s) gFunc(s)-LGM+(idx-2)*pi,...
                        obj.CurveParameter.mid,optimset('Display','off'));
                end
                sSolv = uniquetol(sSolv,1e-6);
                sSolv = sSolv(obj.CurveParameter.isin(sSolv,'tolerance',1e-6));
                point = (1+1i*sSolv) / obj.NormFactor;
            else
                point = nan();
            end
        end

        % Sample
        function points = sample(obj, nPoints)
            [M,N] = size(obj);
            points = cell(M,N);
            for m = 1:M
            for n = 1:N
                points{m,n} = linspace(obj(m,n).Startpoint,...
                                       obj(m,n).Endpoint,nPoints).';
            end
            end

            if M*N == 1
                points = points{:};
            end
        end

        % Transpose
        function r = transpose(obj)
            [M,N] = size(obj);
            r = reshape(obj,N,M);
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
                points = [obj(n).Startpoint , obj(n).Endpoint];
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

        function h = plotLine(obj,varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                fimplicit(@(x,y) obj.Slope .* x + obj.Offset - y,varargin{:})
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
        r = times(obj1,obj2)
        r = mtimes(obj1,obj2)
        r = sum(obj,varargin)
        h = plotMap(obj,logMap,arrowSize,varargin)

    end

   
    %% Vectorizing the object

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Startpoint = obj.Startpoint.(indexOp(1));
            obj.Endpoint = obj.Endpoint.(indexOp(1));
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
            % obj(1,1).Endpoint = 1;
            % Should use
            % obj.Endpoint(1,1) = 1;

            % Ensure object instance is the first argument of call.
            if isempty(obj)
                % Object must be of the correct size
                
                % If rhs is of size 1, then use indices to set size
                if isscalar(varargin{1})
                    sz = [indexOp.Indices{:}];
                    obj = ciat.Edge;
                    obj.Startpoint = nan(sz);
                    obj.Endpoint = nan(sz);
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
                obj.Startpoint = nan(sz);
                obj.Endpoint = nan(sz);
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
                    obj.Startpoint.(indexOp(1)) = rhs.Startpoint;
                    obj.Endpoint.(indexOp(1)) = rhs.Endpoint;
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
            obj.Startpoint.(indexOp) = [];
            obj.Endpoint.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.Edge')
                    newArgs{ix} = varargin{ix}.Startpoint;
                    newArgs2{ix} = varargin{ix}.Endpoint;
                else
                    newArgs{ix} = varargin{ix};
                end
            end
            % out = ciat.Edge(vertcat(dim,newArgs{:}), ...
            %                   vertcat(dim,newArgs2{:}));
            out = ciat.Edge(cat(dim,newArgs{:}), ...
                            cat(dim,newArgs2{:}));
        end

        function varargout = size(obj,varargin)
            % disp('size')
            [varargout{1:nargout}] = size(obj.Startpoint,varargin{:});
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
            obj.Startpoint = reshape(obj.Startpoint,varargin{:});
            obj.Endpoint = reshape(obj.Endpoint,varargin{:});
        end
    end

end