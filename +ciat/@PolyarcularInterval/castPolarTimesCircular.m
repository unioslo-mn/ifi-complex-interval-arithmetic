function outArcs = castPolarTimesCircular(pInt, cInt)

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
        outArcs = ciat.Arc(cInt.Center,cRadius,ciat.RealInterval(-pi,pi));
        return
    end
    
    % If the circle is just a point
    if (cRadius == 0)
        % If the circle center is at zero creat a zero polygon
        if cInt.Center == 0
            outArcs = ciat.Arc(0,0,0);
            return 
        end
        
        % Otherwise multiply the polar by the center and cast to polygon
        outArcs(2,1) = ciat.Arc;
        outArcs(1) = ciat.Arc(0,-pAbs.Infimum,pAngle+pi);
        outArcs(2) = ciat.Arc(0,pAbs.Supremum,pAngle);
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

    % Outer curve
    outCurveArc = ciat.Arc(0 , (abs(cInt.Center)+cRadius)* pAbs.sup, ...
                                cAngle + pAngle );
    
    % Inner curve
    inCurveArc = ciat.Arc(0 , -(abs(cInt.Center)-cRadius)* pAbs.inf , ...
                                cAngle + pAngle + pi );

    
    
    % Calculate parameters for corner elements (c: center, r:radius, a: angle)
    cMax   = cInt.Center * pAbs.sup; % unrotated center of outer corners
    cMin   = cInt.Center * pAbs.inf; % unrotated center of inner corners
    rMax = cRadius * pAbs.sup;      % radius of outer corners
    rMin = cRadius * pAbs.inf;      % radius of inner corners
    aMax = cAngle + pAngle.sup;
    aMin = cAngle + pAngle.inf;
    aShift = asin((rMax-rMin)/(abs(cMax)-abs(cMin))); % shift angle from slope over two circles
    

    % Corner 1: min phase
    start_ang = -pi/2 - aShift;
    stop_ang = 0;
    cornerArc1 = ciat.Arc(cMax*exp(1j*pAngle.inf) , rMax , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMin );
    
    % Corner 2: max phase
    start_ang = 0;
    stop_ang = pi/2 + aShift;
    cornerArc2 = ciat.Arc(cMax*exp(1j*pAngle.sup) , rMax , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMax );
    
    % Corner 3: max phase
    start_ang = pi/2 + aShift;
    stop_ang = pi;
    cornerArc3 = ciat.Arc(cMin*exp(1j*pAngle.sup) , rMin , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMax );

    % Corner 4: min phase
    start_ang = pi;
    stop_ang = 3*pi/2 - aShift;
    cornerArc4 = ciat.Arc(cMin*exp(1j*pAngle.inf) , rMin , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMin);
            
    %% Return as collumn vector and make a polygon
    outArcs = [cornerArc1, outCurveArc cornerArc2, cornerArc3, inCurveArc cornerArc4].';

end


