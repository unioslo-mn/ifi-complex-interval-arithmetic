classdef Arc < matlab.mixin.indexing.RedefinesParen

	properties
        Center          % Arc center as a complex number
        Radius          % Arc radius as a real number (negative means concave arc, zero means vertex)
        ArcAngle		% Arc angle as a real interval (always counter-clockwise)
    end	

    % properties (Access = private)
       % ArcAngleStore        % Storage property for the arc angle
    % end

	properties (Dependent)
    	% ArcAngle		% Arc angle as a real interval (always counter-clockwise)
		Startpoint		% First point of the arc as complex numbers
        Midpoint        % Midpoint of the arc as complex numbers
        Endpoint		% Last point of the arc as complex numbers
		Length			% Length of the arc curve
		GaussMap        % Gauss map of the curve as real interval
		LogGaussMap     % Gauss map of the log curve as real interval
		NormFactor		% Normalization factor to a real 1 centered arc (zero is zero centered)
		CurveParameter	% Parameter interval of the normalized curve as real interval
        Real;           % Projection of the polygonal interval to the real axis
        Imag;           % Projection of the polygonal interval to the imaginary axis
        Abs;            % Projection of the polygonal interval to the absolute value axis
        Angle;          % Projection of the polygonal interval to the angle axis
        Area;           % Area of the arc with its endpoints connected by an edge
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
                    if isa(angles,'double')
                        angles = ciat.RealInterval(angles);
                    end
                    obj.ArcAngle = ciat.Arc.wrapArcAngle(angles);
                otherwise
                    error('incorrect number of input')
            end
        end

        % %% Defining properties
        % 
        % % Set angles (store in the hidden property ArcAngleStore)
        % % wrap angles and split at pi, 
        % % return vector of arcs no matter the input format
        % function obj = set.ArcAngle(obj,angleArray)
        % 
        %     [M,N] = size(obj);
        %     [Ma,Na] = size(angleArray);
        %     assert(M==Ma && N==Na)
        % 
        %     % Initialize ArcAngle property if necessary
        %     if isempty(obj.ArcAngle)
        %         obj.ArcAngleStore = repmat(ciat.RealInterval(0),M,N);
        %     end
        % 
        %     for m = 1:M
        %     for n = 1:N
        %         obj.ArcAngleStore(m,n) = ciat.Arc.wrapArcAngle(angleArray(m,n));
        %     end
        %     end
        % end
        % 
        % % Get points (retrieve ArcAngle from hidden property ArcAngleStore)
        % function value = get.ArcAngle(obj)
        %     value = obj.ArcAngleStore; 
        % end

		%% Dependent properties
        

        % Startpoint
        function value = get.Startpoint(obj)
            angMask = (obj.Radius < 0);
            angInf = obj.ArcAngle.Infimum;
            angSup = obj.ArcAngle.Supremum;
            angStart = angInf.*~angMask + angSup.*angMask;
            value = obj.Center + obj.Radius .* exp(1i*angStart);
        end

        % Midpoint
        function value = get.Midpoint(obj)
            value = obj.Center + obj.Radius .* exp(1i*obj.ArcAngle.Midpoint);
        end

        % Endpoint
        function value = get.Endpoint(obj)
            angMask = (obj.Radius < 0);
            angInf = obj.ArcAngle.Infimum;
            angSup = obj.ArcAngle.Supremum;
            angStop = angInf.*angMask + angSup.*~angMask;
            value = obj.Center + obj.Radius .* exp(1i*angStop);
        end

        % Length
        function value = get.Length(obj)
            value = abs(2*pi*obj.Radius .* obj.ArcAngle.Width/(2*pi));
        end

        % Gauss map angle interval
        function value = get.GaussMap(obj)
            % value = obj.ArcAngleStore;
            value = obj.ArcAngle;
        end

        % Log-Gauss map angle interval
        function value = get.LogGaussMap(obj)
            % value = obj.GaussMap - ciat.RealInterval(angle(obj.Startpoint),...
                                                     % angle(obj.Endpoint));
            isConvex = obj.Radius >= 0;
            startGauss = obj.ArcAngle.Bounds(2-isConvex);
            endGauss = obj.ArcAngle.Bounds(1+isConvex);
            value = ciat.RealInterval( startGauss - angle(obj.Startpoint),...
                                       endGauss - angle(obj.Endpoint) );

            % value.Infimum = ciat.wrapToPi(value.Infimum);
            % value.Supremum= ciat.wrapToPi(value.Supremum);
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
            value = obj.ArcAngle + angle(obj.NormFactor);
        end

        % Real
        function value = get.Real(obj)
            
            % Create boolean masks
            isConcave = obj.Radius < 0;
            crossInf = obj.ArcAngle.isin((isConcave-3)*pi) | ...
                       obj.ArcAngle.isin((isConcave-1)*pi) | ...
                       obj.ArcAngle.isin((isConcave+1)*pi) ;
            crossSup = obj.ArcAngle.isin((isConcave-2)*pi) | ...
                       obj.ArcAngle.isin((isConcave+0)*pi);
            
            % Calculate real bounds of the endpoints
            pntInf = min(real(obj.Startpoint),real(obj.Endpoint));
            pntSup = max(real(obj.Startpoint),real(obj.Endpoint));

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
            crossInf = obj.ArcAngle.isin((isConcave-2.5)*pi) | ...
                       obj.ArcAngle.isin((isConcave-0.5)*pi) | ...
                       obj.ArcAngle.isin((isConcave+1.5)*pi) ;
            crossSup = obj.ArcAngle.isin((isConcave-1.5)*pi) | ...
                       obj.ArcAngle.isin((isConcave+0.5)*pi);
            
            % Calculate real bounds of the endpoints
            pntInf = min(imag(obj.Startpoint),imag(obj.Endpoint));
            pntSup = max(imag(obj.Startpoint),imag(obj.Endpoint));

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
            crossInf = obj.ArcAngle.isin(angOffset+(isConcave-3)*pi) | ...
                       obj.ArcAngle.isin(angOffset+(isConcave-1)*pi) | ...
                       obj.ArcAngle.isin(angOffset+(isConcave+1)*pi);
            crossSup = obj.ArcAngle.isin(angOffset+(isConcave-2)*pi) | ...
                       obj.ArcAngle.isin(angOffset+(isConcave+0)*pi);
            originIn = abs(obj.Radius) > abs(obj.Center);
            
            % Calculate real bounds of the endpoints
            pntInf = min(abs(obj.Startpoint),abs(obj.Endpoint));
            pntSup = max(abs(obj.Startpoint),abs(obj.Endpoint));

            % Calculate abs bounds of the envelope
            envInf = abs( abs(obj.Center) - abs(obj.Radius) );
            envSup = abs( abs(obj.Center) + abs(obj.Radius) );

            % Pick the correct value according to the masks
            % absInf = (pntInf .* (~crossInf | ~originIn) ) + ...
            %          (envInf .* (crossInf & originIn) );
            % absSup = (pntSup .* (~crossSup | ~originIn) ) + ...
            %          (envSup .* (crossSup & originIn) );
            absInf = (pntInf .* ~crossInf + envInf .* crossInf );
            absSup = (pntSup .* ~crossSup + envSup .* crossSup );

            value = ciat.RealInterval(absInf,absSup);
        end
        function value = abs(obj)
            value = obj.Abs;
        end
        
        % Angle
        function value = get.Angle(obj)

            % Calculate real bounds of the endpoints
            isConcave = obj.Radius < 0;
            originIn = abs(obj.Radius) > abs(obj.Center);
            % pntInf = min(angle(obj.Startpoint),angle(obj.Endpoint));
            % pntSup = max(angle(obj.Startpoint),angle(obj.Endpoint));
            angStart = angle(obj.Startpoint);
            angEnd = angle(obj.Endpoint);
            pntInf = (angStart .* ~isConcave + angEnd .* isConcave) .* originIn + ...
                     min(angle(obj.Startpoint),angle(obj.Endpoint)) .* ~originIn;
            pntSup = (angStart .* isConcave + angEnd .* ~isConcave) .* originIn + ...
                     max(angle(obj.Startpoint),angle(obj.Endpoint)) .* ~originIn;
            pntSup = pntSup + originIn .* (pntSup < pntInf)*2*pi;

            % Calculate angle bounds of the envelope
            if obj.Radius <= abs(obj.Center)
                envAngle = asin(abs(obj.Radius ./ obj.Center));
                envInf = angle(obj.Center) - envAngle;
                envSup = angle(obj.Center) + envAngle;
            else
                envInf = -pi;
                envSup = pi;
            end
            
            % Create boolean masks
            crossInf = obj.ArcAngle.isin(envInf+(isConcave-2.5)*pi) | ...
                       obj.ArcAngle.isin(envInf+(isConcave-0.5)*pi) | ...
                       obj.ArcAngle.isin(envInf+(isConcave+1.5)*pi);
            crossSup = obj.ArcAngle.isin(envSup+(isConcave-3.5)*pi) | ...
                       obj.ArcAngle.isin(envSup+(isConcave-1.5)*pi) | ...
                       obj.ArcAngle.isin(envSup+(isConcave+0.5)*pi) | ...
                       obj.ArcAngle.isin(envSup+(isConcave+2.5)*pi);
            
            % angOrder = wrapToPi(obj.ArcAngle.Infimum + isConcave*pi) > ...
            %            wrapToPi(obj.ArcAngle.Supremum + isConcave*pi);

            % Pick the correct value according to the masks
            % angleInf = (pntInf .* ~crossInf) + (envInf .* crossInf);
            % angleSup = (pntSup .* ~crossSup) + (envSup .* crossSup);
            angleInf = (pntInf .* (~crossInf | originIn) ) + ...
                       (envInf .* (crossInf & ~originIn) );
            angleSup = (pntSup .* (~crossSup | originIn) ) + ...
                       (envSup .* (crossSup & ~originIn) );

            % % Correct angles for cases when the origin is inside
            % angleInf = angleInf + originIn*angOrder*2*pi;

            value = ciat.RealInterval(angleInf,angleSup);
        end
        function value = angle(obj)
            value = obj.Angle;
        end

        % Area
        function value = get.Area(obj)
            % Area of the sector
            sectorArea = obj.Radius.^2 .* obj.ArcAngle.Width / 2;
            
            % Area of triangle (center and arc endpoints)
            triBase = abs(obj.Startpoint - obj.Endpoint);
            triSide = abs(obj.Startpoint - obj.Center);
            triHeight = sqrt(triSide.^2 - (triBase/2).^2);
            triArea = triBase .* triHeight / 2;

            % Combine the triangle area and the sector area
            combSign = sign(obj.ArcAngle.Width - pi);
            value = (sectorArea + combSign .* triArea) .* sign(obj.Radius);
        end
        

		%% Other methods

        % Unary negative operator
        function r = uminus(obj)
            r = ciat.Arc(-obj.Center,obj.Radius,obj.ArcAngle+pi);
        end

        % Intersection
        function r = intersection(obj1,obj2)

            [M1,N1] = size(obj1);
            [M2,N2] = size(obj2);
            assert(M1 == M2 && N1 == N2)
            
            if isa(obj2,'ciat.Arc')
                % Find the intersection point of two arcs
                % Source: https://mathworld.wolfram.com/Circle-CircleIntersection.html
                x1 = real(obj1.Center);
                y1 = imag(obj1.Center);
                r1 = obj1.Radius;
                x2 = real(obj2.Center);
                y2 = imag(obj2.Center);
                r2 = obj2.Radius;
                R = sqrt( (x2 - x1).^2 + (y2 - y1).^2 );
                coeff1 = (r1.^2 - r2.^2) ./ (2*R.^2);
                coeff2 = 0.5 * sqrt( 2*(r1.^2+r2.^2)/R.^2 - ...
                                     (r1^2-r2.^2).^2/R.^4 - 1);

                xCoord1 = 0.5*(x1+x2) + coeff1*(x2-x1) - coeff2*(y2-y1);
                xCoord2 = 0.5*(x1+x2) + coeff1*(x2-x1) + coeff2*(y2-y1);
                yCoord1 = 0.5*(y1+y2) + coeff1*(y2-y1) - coeff2*(x1-x2);
                yCoord2 = 0.5*(y1+y2) + coeff1*(y2-y1) + coeff2*(x1-x2);

                % Check if the intersections exist
                mask0 = (abs(r1)+abs(r2)) >= R; 

            elseif isa(obj2,'ciat.Edge')
                % Find the intersection point of an arc and an edge
                % Source: https://mathworld.wolfram.com/Circle-LineIntersection.html
                xc = real(obj1.Center);
                yc = imag(obj1.Center);
                r = obj1.Radius;
                x1 = real(obj2.Startpoint) - xc;
                y1 = imag(obj2.Startpoint) - yc;
                x2 = real(obj2.Endpoint) - xc;
                y2 = imag(obj2.Endpoint) - yc;
                
                dx = x2 - x1;               % Edge vector real
                dy = y2 - y1;               % Edge vector imag
                dr = sqrt(dx.^2 + dy.^2);   % Edge length
                D = x1.*y2 - x2.*y1;
                discr = r.^2 .* dr.^2 - D.^2;

                xCoord1 = (D.*dy - sign(dy).*dx * sqrt(discr)) ./ dr.^2 + xc;
                xCoord2 = (D.*dy + sign(dy).*dx * sqrt(discr)) ./ dr.^2 + xc;
                yCoord1 = (-D.*dx - abs(dy) * sqrt(discr)) ./ dr.^2 + yc;
                yCoord2 = (-D.*dx + abs(dy) * sqrt(discr)) ./ dr.^2 + yc;

                % Check if the intersections exist
                mask0 = discr >= 0; 

            else
                error('Invalid input type')
            end

            % Assign values
            r = nan(M1,N1);
            mask1 = obj1.Real.isin(xCoord1) & obj1.Imag.isin(yCoord1) & ...
                    obj2.Real.isin(xCoord1) & obj2.Imag.isin(yCoord1) & mask0;
            mask2 = obj1.Real.isin(xCoord2) & obj1.Imag.isin(yCoord2) & ...
                    obj2.Real.isin(xCoord2) & obj2.Imag.isin(yCoord2) & mask0;
            if any(mask1 & mask2,'all')
                if M1*N1 ==1
                    r = [xCoord1 + 1i*yCoord1 ; xCoord2 + 1i*yCoord2];
                else
                    warning('Not all results could be returned')
                end
            else
                if any(mask1,'all')
                    r(mask1) = xCoord1(mask1) + 1i * yCoord1(mask1);
                end
                if any(mask2,'all')
                    r(mask2) = xCoord2(mask2) + 1i * yCoord2(mask2);
                end
            end
        end
        function r = cap(obj1,obj2)
            r = intersection(obj1,obj2);
        end

        % Inside
        function r = isin(obj,x)
            % Check if the point is in the circle of the arc
            inCircle = abs(x - obj.Center) <= abs(obj.Radius);

            % Check if the point is in the sector
            inSector = abs(wrapToPi( angle(x - obj.Center) - ...
                                     obj.ArcAngle.Infimum  + ...
                                    (obj.Radius<0)*pi ) ) ...
                            <= obj.ArcAngle.Width;

            % Check if the point is outside the chord connecting the end
            % points
            rotateCoeff = exp(-1j*(obj.ArcAngle.Midpoint + ...
                                   (obj.Radius < 0) * pi) );
            outChord = real(x.* rotateCoeff) >= ...
                       real(obj.Startpoint .* rotateCoeff);

            % Combine conditions
            r = inCircle & inSector & outChord;
        end

        % Is point on the arc
        function r = ison(obj,x)
            % Check if the point is in the circle of the arc
            onCircle = abs(abs(x - obj.Center) - abs(obj.Radius)) < 10*eps;

            % Check if the point is in the sector
            inSector = abs(wrapToPi( angle(x - obj.Center) - ...
                                     obj.ArcAngle.Infimum  + ...
                                    (obj.Radius<0)*pi ) ) ...
                            <= obj.ArcAngle.Width;

            % Combine conditions
            r = onCircle & inSector;
        end

        % Split
        function arcOut = split(arcIn,splitPoint)

            % It only works for a single arc with multiple split points
            assert( numel(arcIn)==1 )

            % Exclude invalid split points 
            splitPoint = splitPoint(~isnan(splitPoint));
            splitPoint = splitPoint(arcIn.ison(splitPoint));
            splitPoint = splitPoint(splitPoint ~= arcIn.Startpoint & ...
                                    splitPoint ~= arcIn.Endpoint );

            % Split arc
            if ~isempty(splitPoint)
                % Get arc parameters
                arcCenter = arcIn.Center;
                arcRadius = arcIn.Radius;
                arcAngInf = arcIn.ArcAngle.Infimum;
                arcAngSup = arcIn.ArcAngle.Supremum;

                % Calculate split angles
                splitAngle = angle(splitPoint-arcCenter) + pi*(arcRadius<0);
                splitAngle = wrapTo2Pi(splitAngle)-2*pi;
                splitAngle = splitAngle( abs(splitAngle-arcAngInf)>10*eps & ...
                                         abs(splitAngle-arcAngSup)>10*eps);
                splitAngle = splitAngle + 2*pi*(splitAngle<arcAngInf);
                splitAngle = sort(splitAngle);
                splitAngle = [arcAngInf ; splitAngle ; arcAngSup ];
                splitAngle = uniquetol(splitAngle,10*eps);

                % Add arc segments
                L = length(splitAngle)-1;
                arcOut(L,1) = ciat.Arc;
                for l = 1:L
                    arcOut(l) = ciat.Arc(arcCenter,arcRadius,...
                                         ciat.RealInterval(splitAngle(l), ...
                                                           splitAngle(l+1)));
                end
            else
                arcOut = arcIn;
            end
        end

        % IsNaN
        function r = isnan(obj)
            r = isnan(obj.Length);
        end 

        % Sample
        function points = sample(obj, nPoints)
            [M,N] = size(obj);
            points = cell(M,N);
            for m = 1:M
            for n = 1:N
                center = obj(m,n).Center;
                radius = obj(m,n).Radius;
                angInf = obj(m,n).ArcAngle.Infimum;
                angSup = obj(m,n).ArcAngle.Supremum;
                points{m,n} = center + radius * ...
                                   exp(1j*linspace(angInf,angSup,nPoints));
                if radius < 0
                    points{m,n} = flip(points{m,n});
                end
            end
            end
        end

        % Transpose
        function r = transpose(obj)
            [M,N] = size(obj);
            r = reshape(obj,N,M);
        end

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
                if ~isnan(obj(n))
                    % Extract parameters
                    cr = real(obj(n).Center);
                    ci = imag(obj(n).Center);
                    r = obj(n).Radius;
                    % Plot arc by angle quadrants between in (-2Pi,2Pi)
                    q = ((-2:0.5:1.5)'+[0,0.5])*pi;
                    if r~=0
                        for idx = 1:length(q)
                            qi = intersection(obj(n).ArcAngle, ...
                                            ciat.RealInterval(q(idx,1),q(idx,2)));
                            if ~isnan(qi)
                                xBound = sort(cr + r*cos(qi.Bounds) );
                                yBound = sort(ci + r*sin(qi.Bounds) );
                                h = [h; fimplicit(@(x,y) ...
                                     (x-cr).^2 + (y-ci).^2 - r^2 , ...
                                            [xBound yBound],varargin{:})];    
                            end
                        end
                    else
                        h = [h; plot(cr,ci,varargin{:},'Marker','.')];
                    end
                end
            end
            if tf == false 
                hold off
            end
        end

        function h = plotCircle(obj,varargin)
            tf = ishold;
            if tf == false 
                clf
            end
            hold on
            h = [];
            for n = 1:length(obj(:))
                fimplicit(@(x,y) (x-real(obj.Center)).^2 + ...
                                 (y-imag(obj.Center)).^2 - ...
                                 obj.Radius.^2 , varargin{:})
            end
            if tf == false 
                hold off
            end
        end

        % Plot Gauss maps
        
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
        h = plotMap(obj, logMap, arrowSize, varargin)

    end

    methods(Static)
        angles = wrapArcAngle(angles)
    end

%% Vectorizing the object

    methods (Access=protected)
        function varargout = parenReference(obj, indexOp)
            % disp('parenReference')
            obj.Center = obj.Center.(indexOp(1));
            obj.Radius = obj.Radius.(indexOp(1));
            obj.ArcAngle = obj.ArcAngle.(indexOp(1));
            % obj.ArcAngleStore = obj.ArcAngleStore.(indexOp(1));
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
                    obj.Center = nan(sz);
                    obj.Radius = nan(sz);
                    % ang(sz) = ciat.RealInterval;
                    % obj.ArcAngleStore = ang;
                    % obj.ArcAngleStore = repmat(ciat.RealInterval,sz);
                    % obj.ArcAngle = ang;
                    obj.ArcAngle = repmat(ciat.RealInterval,sz);
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
                obj.Center = nan(sz);
                obj.Radius = nan(sz);
                ang(sz(1),sz(2)) = ciat.RealInterval;
                % obj.ArcAngleStore = ang;
                obj.ArcAngle = ang;
                return;
            end
            if numel(indexOp) == 1
                if isscalar(indexOp(1))
                    assert(nargin==3);
                    rhs = varargin{1};
                    % If rhs is not an interval, then convert it to one.
                    if ~isa(rhs, 'ciat.Arc')
                        rhs = ciat.RealInterval(rhs);
                    end
                    obj.Center.(indexOp(1)) = rhs.Center;
                    obj.Radius.(indexOp(1)) = rhs.Radius;
                    % obj.ArcAngleStore.(indexOp(1)) = rhs.ArcAngleStore;
                    obj.ArcAngle.(indexOp(1)) = rhs.ArcAngle;
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
            obj.Center.(indexOp) = [];
            obj.Radius.(indexOp) = [];
        end
    end

    methods (Access=public)
        function out = cat(dim,varargin)
            numCatArrays = nargin-1;
            newArgs = cell(numCatArrays,1);
            newArgs2 = cell(numCatArrays,1);
            newArgs3 = cell(numCatArrays,1);
            for ix = 1:numCatArrays
                if isa(varargin{ix},'ciat.Arc')
                    newArgs{ix} = varargin{ix}.Center;
                    newArgs2{ix} = varargin{ix}.Radius;
                    % newArgs3{ix} = varargin{ix}.ArcAngleStore;
                    newArgs3{ix} = varargin{ix}.ArcAngle;
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
            % obj.ArcAngleStore = reshape(obj.ArcAngleStore,varargin{:});
            obj.ArcAngle = reshape(obj.ArcAngle,varargin{:});
        end
    end
end