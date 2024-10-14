function outObj = cast(inObj,inObj2)

    arguments
        inObj
        inObj2               (:,:)   = []
    end

% Cast complex intervals of other types to polyarx interval type
   
    switch class(inObj)
        case 'double'
            arx = [real(inObj),imag(inObj),0,pi];
        case 'ciat.RealInterval'
            arx = zeros(2,4);
            arx(1,:) = [inObj.Infimum,0,0,-pi/2];
            arx(2,:) = [inObj.Supremum,0,0,pi/2];
        case 'ciat.RectangularInterval'
            inReal = [inObj.Real];
            inImag = [inObj.Imag];
            
            arx = zeros(4,4);
            arx(1,:) = [inReal.Infimum,inImag.Infimum,0,-pi/2];
            arx(2,:) = [inReal.Supremum,inImag.Infimum,0,0];
            arx(3,:) = [inReal.Supremum,inImag.Supremum,0,pi/2];
            arx(4,:) = [inReal.Infimum,inImag.Supremum,0,pi];
            
        case 'ciat.CircularInterval'
            if isempty(inObj2)
                arx = [real(inObj.Center),imag(inObj.Center),inObj.Radius,pi];
            else
                if isa(inObj2,'ciat.PolarInterval')
                    arx = castPolarTimesCircular(inObj2,inObj);
                else
                    error('Invalid input type at position 2')
                end
            end
            
        case 'ciat.PolarInterval'
            if isempty(inObj2)
                % Extract parameters
                absInf = inObj.Abs.Infimum;
                absSup = inObj.Abs.Supremum;
                angInf = ciat.wrapToPi(inObj.Angle.Infimum);
                angSup = ciat.wrapToPi(inObj.Angle.Supremum);
                % angSup = angSup + (angSup<angInf)*2*pi; %???
                
                % Calculate vertex locations
                v1 = absSup * exp(1j*angSup);
                v2 = absInf * exp(1j*angSup);
                v3 = absInf * exp(1j*angInf);
                v4 = absSup * exp(1j*angInf);
    
                % Calculate normal angles vertices
                a1 = ciat.wrapToPi(angSup+pi/2);
                a2 = ciat.wrapToPi(angle(v2+v3)+pi);
                a3 = ciat.wrapToPi(angInf-pi/2);
                a4 = angInf;
    
                % Generate arx
                arx = zeros(5,4);
                arx(1,:) = [ 0 , 0 , absSup , angSup ];
                arx(2,:) = [ real(v1) , imag(v1) , 0 , a1 ];
                arx(3,:) = [ real(v2) , imag(v2) , 0 , a2 ];
                arx(4,:) = [ real(v3) , imag(v3) , 0 , a3 ];
                arx(5,:) = [ real(v4) , imag(v4) , 0 , a4 ];
            else
                if isa(inObj2,'ciat.CircularInterval')
                    arx = castPolarTimesCircular(inObj,inObj2);
                else
                    error('Invalid input type at position 2')
                end
            end
        case 'ciat.PolyarcularInterval'
            cxObj = inObj.convexify;
            arcs = [cxObj.Arcs{:} ; cxObj.Vertices{:}];

            % Set polyarx
            arx = [real(arcs.Center) , imag(arcs.Center), ...
                   arcs.Radius , arcs.ArcAngle.Supremum];
            [~,idx] = sort(arx(:,4));
            arx = arx(idx,:);
        otherwise
            error('Invalid input type at position 1')
    end  
    outObj = ciat.PolyarxInterval(arx);       
end

%% Utility function for casting the product of a polar and circular interval

function outArx = castPolarTimesCircular(pInt, cInt)

    % SETUP AND CHECKS
    arguments
        pInt   ciat.PolarInterval
        cInt   ciat.CircularInterval
    end
         
    pAng = pInt.Angle;
    pAbs = pInt.Abs;
    cCen = cInt.Center;
    cAng = angle(cCen);
    cRad = cInt.Radius;
    
    % Handle exceptions
    
    % If the polar interval is just a point
    if pInt.Area == 0
        outArx = [real(cCen),imag(cCen),cRad,pi];
        return
    end
    
    % If the circle is just a point
    if cInt.Area == 0
        % If the circle center is at zero creat a zero polygon
        if cCen == 0
            outArx = [0,0,0,pi];
            return 
        end

        xInt = ciat.PolyarxInterval( pInt * cCen );
        outArx = xInt.Arx;
        
        return
    end

    if cInt.isin(0) 
        error('This algorithm does not work correctly for circles including the origin.')
    end

    % If the product is expected to be concave give a warning
    if (width(angle(pInt) + angle(cInt)) > pi) || ...
       (abs(cCen) <= cRad)
        warning('Circle times polar may not be convex.');
        % This must be fixed for apodization windows where weights can be
        % 0, e.g., hann window.
        % Also a case if: if ( abs(cCen) <= cRad)
    end

    % Generate product shape from sampled arcs

    % Initialize output
    outArx = zeros(5,4);

    % Calculate parameters for corner elements (c: center, r:radius, a: angle)
    rotAng = asin( cRad / abs(cCen) ); % shift angle from slope over two circles

    % Segment 1: Outer curve (top center)
    outArx(1,:) = [0,0,(abs(cCen)+cRad)* pAbs.sup,...
                        cAng + pAng.sup];

    % Segment 2: max amplitude, max phase (top left) corner
    cArx = cCen * pAbs.sup * exp(1j*pAng.sup);
    outArx(2,:) = [ real(cArx) , imag(cArx) , cRad * pAbs.sup , ...
                                cAng + pAng.sup + rotAng + pi/2 ];
    
    % Segment 3: min amplitude, max phase (bottom left) corner
    cArx = cCen * pAbs.inf * exp(1j*pAng.sup);
    outArx(3,:) = [ real(cArx) , imag(cArx) , cRad * pAbs.inf , ...
                                cAng + pAng.sup - pAng.width/2 + pi ];

    % Segment 4: min amplitude, min phase (bottom right) corner
    cArx = cCen * pAbs.inf *exp(1j*pAng.inf);
    outArx(4,:) = [ real(cArx) , imag(cArx) , cRad * pAbs.inf , ...
                                cAng + pAng.inf - rotAng - pi/2 ];

    % Segment 5: max amplitude, min phase (top right) corner
    cArx = cCen * pAbs.sup * exp(1j*pAng.inf);
    outArx(5,:) = [ real(cArx) , imag(cArx) , cRad * pAbs.sup , ...
                                cAng + pAng.inf ];

    % Wrap supremum angles to Pi
    outArx(:,4) = ciat.wrapToPi(outArx(:,4));

end