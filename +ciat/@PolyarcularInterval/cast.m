function outObj = cast(inObj,inObj2)

% Cast complex intervals of other types to polyarcular interval type
%
% This function takes one or more complex intervals of another type
% and creates the smallest inclusive convex interval(s) of the polyarcular
% interval type with the given tolerance. 
% _________________________________________________________________________
% USAGE        
%   outObj = ciat.Polyarcular.Interval.cast(inObj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   inObj       : object of one of the following types:
%                   - ciat.RectangularInterval
%                   - ciat.CircularInterval               
%                   - ciat.PolarInterval
%                   - ciat.PolygonalInterval
% _________________________________________________________________________
% OPTIONS
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

    % [M,N] = size(inObj);
    
    arguments
        inObj
        inObj2               (:,:)   = []
    end

    switch class(inObj)
        case 'double'
            outArcs = ciat.Arc(inObj,0,ciat.RealInterval(-pi,pi));
        case 'ciat.RealInterval'
            outArcs(1) = ciat.Arc(inObj.Infimum,0,...
                                    ciat.RealInterval(-3*pi/2,-pi/2));
            outArcs(2) = ciat.Arc(inObj.Supremum,0,...
                                    ciat.RealInterval(-pi/2,pi/2));
        case 'ciat.RectangularInterval'
            inReal = [inObj.Real];
            inImag = [inObj.Imag];
            
            outArcs(4,1) = ciat.Arc;
            outArcs(1) = ciat.Arc(complex([inReal.Infimum],[inImag.Infimum]),...
                                    0,ciat.RealInterval(-pi,-pi/2));
            outArcs(2) = ciat.Arc(complex([inReal.Supremum],[inImag.Infimum]),...
                                    0,ciat.RealInterval(-pi/2,0));
            outArcs(3) = ciat.Arc(complex([inReal.Supremum],[inImag.Supremum]),...
                                    0,ciat.RealInterval(0,pi/2));
            outArcs(4) = ciat.Arc(complex([inReal.Infimum],[inImag.Supremum]),...
                                    0,ciat.RealInterval(pi/2,pi));
            
        case 'ciat.CircularInterval'
            if isempty(inObj2)
                inCenter = [inObj.Center];
                inRadius = [inObj.Radius];
                outArcs = ciat.Arc(inCenter,inRadius,ciat.RealInterval(-pi,pi));
            else
                if isa(inObj2,'ciat.PolarInterval')
                    outArcs = castPolarTimesCircular(inObj2,inObj);
                else
                    error('Invalid input type at position 2')
                end
            end
            
        case 'ciat.PolarInterval'
            if isempty(inObj2)
                inAbs = [inObj.Abs];
                inAngle = [inObj.Angle];
                maxAbs = [inAbs.Supremum];
                minAbs = [inAbs.Infimum];
                outArcs(2,1) = ciat.Arc;
                outArcs(1) = ciat.Arc(0,-minAbs,inAngle+pi);
                outArcs(2) = ciat.Arc(0,maxAbs,inAngle);
            else
                if isa(inObj2,'ciat.CircularInterval')
                    outArcs = castPolarTimesCircular(inObj,inObj2);
                else
                    error('Invalid input type at position 2')
                end
            end
        case 'ciat.PolygonalInterval'
            N = inObj.PointCount;
            outArcs(N,1) = ciat.Arc;
            for n = 1:N
                outArcs(n) = ciat.Arc(inObj.Points(n),0,0);
            end
        otherwise
            error('Invalid input type at position 1')
    end  
    outObj = ciat.PolyarcularInterval(outArcs);       
    % outObj = reshape(outObj,M,N);  
end

%% Utility function to cast the product of a polar and circular interval

function outArcs = castPolarTimesCircular(pInt, cInt)

    % SETUP AND CHECKS
    arguments
        pInt   ciat.PolarInterval
        cInt   ciat.CircularInterval
    end
         
    pAngle = pInt.Angle;
    pAbs   = pInt.Abs;
    cAngle   = angle(cInt.Center);
    cRadius     = cInt.Radius;
    
    % Handle exceptions
    
    % If the polar interval is just a point
    if pInt.Area == 0
        outArcs = ciat.Arc(cInt.Center,cRadius,ciat.RealInterval(-pi,pi));
        return
    end
    
    % If the circle is just a point
    if cInt.Area == 0
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

    % Generate product shape from sampled arcs

    % Initialize output
    outArcs(6,1) = ciat.Arc;

    % Calculate parameters for corner elements (c: center, r:radius, a: angle)
    cMax   = cInt.Center * pAbs.sup; % unrotated center of outer corners
    cMin   = cInt.Center * pAbs.inf; % unrotated center of inner corners
    rMax = cRadius * pAbs.sup;      % radius of outer corners
    rMin = cRadius * pAbs.inf;      % radius of inner corners
    aMax = cAngle + pAngle.sup;
    aMin = cAngle + pAngle.inf;
    aShift = asin((rMax-rMin)/(abs(cMax)-abs(cMin))); % shift angle from slope over two circles
    
    % Segment 1: Outer curve (top center)
    outArcs(1) = ciat.Arc(0 , (abs(cInt.Center)+cRadius)* pAbs.sup, ...
                                cAngle + pAngle );
    
    % Segment 2: max amplitude, max phase (top left) corner
    start_ang = 0;
    stop_ang = pi/2 + aShift;
    outArcs(2) = ciat.Arc(cMax*exp(1j*pAngle.sup) , rMax , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMax );
    
    % Segment 3: min amplitude, max phase (bottom left) corner
    start_ang = pi/2 + aShift;
    stop_ang = pi;
    outArcs(3) = ciat.Arc(cMin*exp(1j*pAngle.sup) , rMin , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMax );

    % Segment 4: inner curve (bottom center)
    outArcs(4) = ciat.Arc(0 , -(abs(cInt.Center)-cRadius)* pAbs.inf , ...
                                cAngle + pAngle + pi );
    
    % Segment 5: min amplitude, min phase (bottom right) corner
    start_ang = pi;
    stop_ang = 3*pi/2 - aShift;
    outArcs(5) = ciat.Arc(cMin*exp(1j*pAngle.inf) , rMin , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMin);

    % Segment 6: max amplitude, min phase (top right) corner
    start_ang = -pi/2 - aShift;
    stop_ang = 0;
    outArcs(6) = ciat.Arc(cMax*exp(1j*pAngle.inf) , rMax , ...
                          ciat.RealInterval(start_ang,stop_ang) + aMin );
            
end




