function points = castPolarTimesCircular(pInt, cInt, dR)

    %% SETUP AND CHECKS
    arguments
        pInt   ciat.PolarInterval
        cInt   ciat.CircularInterval
        dR      double = 1e-6 
    end
         
    pAngleSup = pInt.Angle.Supremum;
    pAngleInf = pInt.Angle.Infimum;
    pAbsSup   = pInt.Abs.Supremum;
    pAbsInf   = pInt.Abs.Infimum;
    cAngle   = angle(cInt.Center);
    cRadius     = cInt.Radius;
    
    %% Handle exceptions
    
    % If the polar interval is just a point
    if pAngleInf == pAngleSup && pAbsInf == pAbsSup
        % Multiply the circle by it and then turn the circle into a polygon
        temp_polar = pAbsSup * exp(1i*pAngleInf);
        temp_circle = ciat.CircularInterval(cInt.Center * temp_polar, ...
                                            cInt.radius * abs(temp_polar));
        polygon = ciat.PolygonalInterval(temp_circle,'tolerance',dR);
        points = polygon.Points{:};
        return
    end
    
    % If the circle is just a point
    if (cRadius == 0)
        % If the circle center is at zero creat a zero polygon
        if cInt.Center == 0
            polygon = ciat.PolygonalInterval(0);
            points = polygon.Points{:};
            return 
        end
        
        % Otherwise multiply the polar by the center and cast to polygon
        temp_polar = pInt;
        temp_polar.Abs = temp_polar.Abs * abs(cInt.Center);
        temp_polar.Angle = temp_polar.Angle + angle(cInt.Center);
        polygon = ciat.PolygonalInterval(temp_polar,'tolerance', dR);
        points = polygon.Points{:};
        return
    end

    if cInt.isin(0) 
        error('This algorithm does not work correctly for circles including the origin.')
    end

    % If the product is expected to be concave give a warning
    if (width(angle(pInt) + angle(cInt)) > pi) || ...
       (abs(cInt.Center) <= cRadius)
        warning('Circle times polar may not be convex, although assumed to be (angle interval larger than pi)! Consider taking convex hull.');
        % This must be fixed for apodization windows where weights can be
        % 0, e.g., hann window.
        % Also a case if: if ( abs(cInt.Center) <= cRadius)
    end
    %% Generate product shape from sampled arcs
    
    % Corners (C = centers, R = radius)
    C_u   = cInt.Center * pAbsSup; % unrotated center of outer corners
    C_l   = cInt.Center * pAbsInf; % unrotated center of inner corners
    R_max = cRadius * pAbsSup;      % radius of outer corners
    R_min = cRadius * pAbsInf;      % radius of inner corners
    shift_ang = asin((R_max-R_min)/(abs(C_u)-abs(C_l))); % shift angle from slope over two circles
    
    % Outer curve
    max_ang = cAngle + pAngleSup;
    min_ang = cAngle + pAngleInf;
    R_outer = (cRadius + abs(cInt.Center))*pAbsSup; % radius, centered on (0+i0)
    
    % Outer curve (between two outer corner circles)
    angularResolution = 2*acos(R_outer/(R_outer+dR));
    n_points = ceil((max_ang - min_ang)/angularResolution) + 1; % at least two points
    
    angs = linspace(min_ang, max_ang, n_points);
    angs = angs(2:end-1); % don't need overlapping points w/ corners
    outer_curve_points = (R_outer + dR) * (cos(angs) + 1j*sin(angs));
    
    %% Corner 1 & 2: max R
    
    % Corner 1: min phase
    start_ang = -pi/2 - shift_ang;
    stop_ang = 0;
    
    angularResolution = 2*acos(R_max/(R_max+dR));
    n_points = ceil((stop_ang - start_ang)/angularResolution) + 1; % at least two points
    
    angs = linspace(start_ang, stop_ang, n_points); 
    Corner1_points = (R_max + dR) * (cos(angs) + 1j*sin(angs)) * exp(1j*min_ang) + C_u*exp(1j*pAngleInf);
    
    % Corner 2: max phase
    start_ang = 0;
    stop_ang = pi/2 + shift_ang;
    
    angs = linspace(start_ang, stop_ang, n_points); 
    Corner2_points = (R_max + dR) * (cos(angs) + 1j*sin(angs)) * exp(1j*(max_ang)) + C_u*exp(1j*pAngleSup);
        
    %% Corner 3 & 4: min R
    G = (pAngleSup - pAngleInf)/2; % Mid angle
    
    % Corner 3: max phase
    start_ang = pi/2 + shift_ang;
    stop_ang = pi - G;

    angularResolution = 2*acos(R_min/(R_min+dR));
    n_points = ceil((stop_ang - start_ang)/angularResolution) + 1; % at least two points
    
    angs = linspace(start_ang, stop_ang, n_points);     
    Corner3_points = (R_min + dR) * (cos(angs) + 1j*sin(angs)) * exp(1j*(max_ang)) + C_l*exp(1j*pAngleSup);
    
    % Corner 4: min phase
    start_ang = pi + G ;
    stop_ang = 3*pi/2 - shift_ang;
    
    angs = linspace(start_ang, stop_ang, n_points);     
    Corner4_points = (R_min + dR) * (cos(angs) + 1j*sin(angs)) * exp(1j*(min_ang)) + C_l*exp(1j*pAngleInf);
        
    %% Return as collumn vector and make a polygon
    points = [Corner1_points, outer_curve_points Corner2_points, Corner3_points, Corner4_points].';
end


