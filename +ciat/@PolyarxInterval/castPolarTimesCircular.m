function outArx = castPolarTimesCircular(pInt, cInt)

    %% SETUP AND CHECKS
    arguments
        pInt   ciat.PolarInterval
        cInt   ciat.CircularInterval
    end
         
    pAngle = pInt.Angle;
    pAbs   = pInt.Abs;
    cAngle   = angle(cInt.Center);
    cRadius     = cInt.Radius;
    
    %% Handle exceptions
    
    % If the polar interval is just a point
    if pAngle.inf == pAngle.sup && pAbs.inf == pAbs.sup
        outArx = [real(cInt.Center),imag(cInt.Center),cRadius,pi];
        return
    end
    
    % If the circle is just a point
    if (cRadius == 0)
        % If the circle center is at zero creat a zero polygon
        if cInt.Center == 0
            outArx = [0,0,0,pi];
            return 
        end
        
        % Otherwise multiply the polar by the center and cast to polygon
        absInf = inObj.Abs.Infimum;
        absSup = inObj.Abs.Supremum;
        angInf = ciat.wrapToPi(inObj.Angle.Infimum);
        angSup = ciat.wrapToPi(inObj.Angle.Supremum);
        angSup = angSup + (angSup<angInf)*2*pi;
        
        % Calculate vertex locations
        v1 = absSup * exp(1j*angSup);
        v2 = absInf * exp(1j*angSup);
        v3 = absInf * exp(1j*angInf);
        v4 = absSup * exp(1j*angInf);

        % Calculate normal angles vertices
        a1 = ciat.wrapToPi(angSup+pi/2);
        a2 = wrapToPi(angle(v2+v3)+pi);
        a3 = ciat.wrapToPi(angInf-pi/2);
        a4 = angInf;

        % Generate outArx
        outArx = zeros(5,4);
        outArx(1,:) = [ 0 , 0 , absSup , angSup ];
        outArx(2,:) = [ real(v1) , imag(v1) , 0 , a1 ];
        outArx(3,:) = [ real(v2) , imag(v2) , 0 , a2 ];
        outArx(4,:) = [ real(v3) , imag(v3) , 0 , a3 ];
        outArx(5,:) = [ real(v4) , imag(v4) , 0 , a4 ];
        return
    end

    if cInt.isin(0) 
        error('This algorithm does not work correctly for circles including the origin.')
    end

    % If the product is expected to be concave give a warning
    if (width(angle(pInt) + angle(cInt)) > pi) || ...
       (abs(cInt.Center) <= cRadius)
        warning('Circle times polar may not be convex.');
        % This must be fixed for apodization windows where weights can be
        % 0, e.g., hann window.
        % Also a case if: if ( abs(cInt.Center) <= cRadius)
    end

    %% Generate product shape from sampled arcs

    % Initialize output
    outArx = zeros(5,4);

    % Outer curve
    % outCurveArc = ciat.Arc(0 , (abs(cInt.Center)+cRadius)* pAbs.sup, ...
                                % cAngle + pAngle );
    outArx(2,:) = [0,0,(abs(cInt.Center)+cRadius)* pAbs.sup,...
                        cAngle + pAngle.sup];
       
    
    % Calculate parameters for corner elements (c: center, r:radius, a: angle)
    aShift = asin( cInt.Radius / abs(cInt.Center) ); % shift angle from slope over two circles
    
    % Corner 1: min phase
    cArx = cInt.Center * pAbs.sup * exp(1j*pAngle.inf);
    outArx(1,:) = [ real(cArx) , imag(cArx) , cRadius * pAbs.sup , ...
                                cAngle + pAngle.inf ];
    
    % Corner 2: max phase
    cArx = cInt.Center * pAbs.sup * exp(1j*pAngle.sup);
    outArx(3,:) = [ real(cArx) , imag(cArx) , cRadius * pAbs.sup , ...
                                cAngle + pAngle.sup + aShift + pi/2 ];
    
    % Corner 3: max phase
    cArx = cInt.Center * pAbs.inf * exp(1j*pAngle.sup);
    outArx(4,:) = [ real(cArx) , imag(cArx) , cRadius * pAbs.inf , ...
                                cAngle + pAngle.sup - pAngle.width/2 + pi ];

    % Corner 4: min phase
    cArx = cInt.Center * pAbs.inf *exp(1j*pAngle.inf);
    outArx(5,:) = [ real(cArx) , imag(cArx) , cRadius * pAbs.inf , ...
                                cAngle + pAngle.inf - aShift - pi/2 ];

    outArx(:,4) = ciat.wrapToPi(outArx(:,4));


end


