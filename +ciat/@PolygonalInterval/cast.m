function outObj = cast(inObj,inObj2,options)

% Cast complex intervals of other types to polygonal interval type
%
% This function takes one or more complex intervals of another type
% and creates the smallest inclusive convex interval(s) of the polygonal
% interval type with the given tolerance. It also allows casting
% the combination of polar and a circular intervals into their product 
% as polygonal intervals.
% _________________________________________________________________________
% USAGE        
%   outObj = ciat.PolygonalInterval.cast(inObj,inObj2,'tolerance',tol)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   inObj       : object of one of the following types:
%                   - ciat.RectangularInterval
%                   - ciat.CircularInterval               
%                   - ciat.PolarInterval
% _________________________________________________________________________
% OPTIONS
%   inObj2       : object of one of the ciat.CircularInterval type
%                  (works only if inObj is of ciat.PolarInterval type)
%   tolerance    : the maximum allowed deviation from the boundary of the
%                  input interval (except where the input is concave)
% _________________________________________________________________________
% EXAMPLES
%   polyInt = ciat.PolygonalInterval(ciat.RectangularInterval(1,3,2,4));
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    arguments
       inObj
       inObj2               (:,:)   = []
       options.tolerance    (1,1)   {mustBeNumeric}     = 1e-6
    end
    
    [M,N] = size(inObj);
    dR = options.tolerance;
    
    switch class(inObj)
        case 'double'
            outPoints = inObj;
        case 'ciat.RealInterval'
            outPoints(1) = inObj.Infimum;
            outPoints(2) = inObj.Supremum;
        case 'ciat.RectangularInterval'
            inReal = [inObj.Real];
            inImag = [inObj.Imag];
            
            outPoints = zeros(4,1);
            outPoints(1) = complex([inReal.Infimum] , [inImag.Infimum]);
            outPoints(2) = complex([inReal.Infimum] , [inImag.Supremum]);
            outPoints(3) = complex([inReal.Supremum] , [inImag.Supremum]);
            outPoints(4) = complex([inReal.Supremum] , [inImag.Infimum]);
            
        case 'ciat.CircularInterval'
            inCenter = [inObj.Center];
            inRadius = [inObj.Radius];
            angRes = 2*acos(inRadius ./ (inRadius + dR));
            cntPoints = ceil(2*pi/angRes); 
            angles = linspace( 0, 2*pi, cntPoints+1 )';
            
            outPoints = (inRadius + dR) .* exp(1j*angles(2:end)) + inCenter;
            
        case 'ciat.PolarInterval'
            if isempty(inObj2)
                inAbs = [inObj.Abs];
                inAngle = [inObj.Angle];
                maxAbs = [inAbs.Supremum];
                minAbs = [inAbs.Infimum];
                maxAngle = [inAngle.Supremum];
                minAngle = [inAngle.Infimum];
                angRes = 2 * acos( maxAbs ./ (maxAbs + dR) );
                cntPoints = ceil((maxAngle - minAngle) ./ angRes) + 1;
                halfCircle = [inAngle.Width] >= pi;
                fullCircle = [inAngle.Width] >= 2*pi;

                % Put together outer arc and inner points: COUNTER-CLOCKWISE ORDER
                angs = linspace(minAngle, maxAngle, cntPoints);
                pL = minAbs * exp( 1j * angs(1) ); % inner corner 1
                rH = (maxAbs + dR) * exp(1j * angs); % outer arc
                pH = minAbs * exp( 1j * angs(end));  % inner corner 2

                % Checks and corrections in case a half or full circle is made 
                maxAngle(fullCircle) = 2*pi;
                minAngle(fullCircle) = 0;
                pL(halfCircle) = [];
                pH(halfCircle) = [];

                % Compile points
                outPoints = [pL,rH, pH].';
            else
                if isa(inObj2,'ciat.CircularInterval')
                    outPoints = timesPolarCircular(inObj,inObj2,dR);
                else
                    error('Invalid input type at position 2')
                end
            end
        otherwise
            error('Invalid input type at position 1')
    end  
    outObj = ciat.PolygonalInterval(outPoints);       
    outObj = reshape(outObj,M,N);  
end

%% Utility function

function points = timesPolarCircular(pInt, cInt, dR)

    %% SETUP AND CHECKS
    arguments
        pInt   ciat.PolarInterval
        cInt   ciat.CircularInterval
        dR  double = 1e-6 
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
        return
    end
    
    % If the circle is just a point
    if (cRadius == 0)
        % If the circle center is at zero creat a zero polygon
        if cInt.Center == 0
            polygon = ciat.PolygonalInterval(0);
            return 
        end
        
        % Otherwise multiply the polar by the center and cast to polygon
        temp_polar = pInt;
        temp_polar.Abs = temp_polar.Abs * abs(cInt.Center);
        temp_polar.Angle = temp_polar.Angle + angle(cInt.Center);
        polygon = ciat.PolygonalInterval(temp_polar,'tolerance', dR);
        return
    end

    % If the product is expected to be concave give a warning
    if (width(angle(pInt) + angle(cInt)) > pi) || (abs(cInt.Center) <= cRadius)
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


