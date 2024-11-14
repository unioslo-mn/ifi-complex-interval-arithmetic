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
            if isempty(inObj2)
                inCenter = [inObj.Center];
                inRadius = [inObj.Radius];
                angRes = 2*acos(inRadius ./ (inRadius + dR));
                cntPoints = ceil(2*pi/angRes); 
                angles = linspace( 0, 2*pi, cntPoints+1 )';
                
                %outPoints = (inRadius + dR) .* exp(1j*angles(2:end)) + inCenter;
                outPoints = (inRadius/cos(2*pi/(2*cntPoints))) .* ...
                                exp(1j*angles(2:end)) + inCenter;
            else
                if isa(inObj2,'ciat.PolarInterval')
                    outPoints = castPolarTimesCircular(inObj2,inObj,dR);
                else
                    error('Invalid input type at position 2')
                end
            end
            
        case 'ciat.PolarInterval'
            if isempty(inObj2)
                outPoints = castPolar(inObj,dR);
                
            else
                if isa(inObj2,'ciat.CircularInterval')
                    outPoints = castPolarTimesCircular(inObj,inObj2,dR);
                else
                    error('Invalid input type at position 2')
                end
            end
        otherwise
            error('Invalid input type at position 1')
    end  
    outObj = ciat.PolygonalInterval(outPoints);       
    % outObj = reshape(outObj,M,N);  
end

%%

function points = castPolar(inObj,dR)

    % Extract parameters
    inAbs = inObj.Abs;
    inAngle = inObj.Angle;
    maxAbs = inAbs.Supremum;
    minAbs = inAbs.Infimum;
    maxAng = inAngle.Supremum;
    minAng = inAngle.Infimum;

    if minAng == maxAng % Degenerate interval with zero angle width

        % Put together outer arc and inner points: COUNTER-CLOCKWISE ORDER
        points = [minAbs * exp(1j*minAng) , maxAbs * exp(1j*minAng)];
    else
        % Put together outer arc and inner points: COUNTER-CLOCKWISE ORDER
        pL = minAbs * exp( 1j * minAng ); % inner corner 1
        rH = ciat.Arc(0,maxAbs,minAng,maxAng).polyWrap(dR);
        pH = minAbs * exp( 1j * maxAng);  % inner corner 2
    
        % Compile points
        if inAngle.Width >= pi
            points = rH;
        else
            points = [pL ; rH ; pH]; 
        end

    end
end


%% Utility function for casting the product of a polar and circular interval

function points = castPolarTimesCircular(pInt, cInt, dR)

    % SETUP AND CHECKS
    arguments
        pInt   ciat.PolarInterval
        cInt   ciat.CircularInterval
        dR      double = 1e-6 
    end
         
    pAng = pInt.Angle;
    pAbs = pInt.Abs;
    cCen = cInt.Center;
    cAng = angle(cCen);
    cRad = cInt.Radius;
    
    % Handle exceptions
    
    % If the polar interval is just a point
    if pInt.Area == 0
        % Multiply the circle by it and then turn the circle into a polygon
        polygon = ciat.PolygonalInterval(cInt,'tolerance',dR);
        points = polygon.Points{:} * pAbs.inf * exp(1j*pAng.inf);
        return
    end
    
    % If the circle is just a point
    if cInt.Area == 0
        % If the circle center is at zero create a zero polygon
        if cCen == 0
            points = 0;
            return 
        end
        
        % Otherwise multiply the polar by the center and cast to polygon
        polygon = ciat.PolygonalInterval(pInt,'tolerance', dR);
        points = polygon.Points{:} * cCen;
        return
    end

    % Check if the circle includes the origin (results unexpected behavior)
    if cInt.isin(0) 
        error('This algorithm does not work correctly for circles including the origin.')
    end

    % If the product is expected to be concave give a warning
    if (width(angle(pInt) + angle(cInt)) > pi) || ...
       (abs(cCen) <= cRad)
        warning('Circle times polar may not be convex, although assumed to be (angle interval larger than pi)! Consider taking convex hull.');
        % This must be fixed for apodization windows where weights can be
        % 0, e.g., hann window.
        % Also a case if: if ( abs(cCen) <= cRad)
    end
    
    % Generate product shape from sampled arcs
    
    % Calculate outer angles of the result shape
    minAng = cAng + pAng.inf;        % unrotated center of outer corners
    maxAng = cAng + pAng.sup;        % unrotated center of outer corners
    rotAng = asin(cRad/abs(cCen)); % rotation angle from slope over two circles

    % Segment 1: Outer curve (top center)
    arcCenter = 0;
    arcAngInf = minAng;
    arcAngSup = maxAng;
    arcRad = (abs(cCen) + cRad) * pAbs.sup;
    points{1} = ciat.Arc(arcCenter,arcRad,arcAngInf,arcAngSup).polyWrap(dR);
    points{1} = points{1}(2:end-1);% don't need overlapping points w/ corners
    
    % Segment 2: max amplitude, max phase (top left) corner
    arcCenter = cCen*pAbs.sup*exp(1j*pAng.sup);
    arcAngInf = maxAng;
    arcAngSup = maxAng + pi/2 + rotAng;
    arcRad = cRad * pAbs.sup;
    points{2} = ciat.Arc(arcCenter,arcRad,arcAngInf,arcAngSup).polyWrap(dR);
        
    % Segment 3: min amplitude, max phase (bottom left) corner
    arcCenter = cCen*pAbs.inf*exp(1j*pAng.sup);
    arcAngInf = maxAng + pi/2 + rotAng;
    arcAngSup = maxAng + pi - (pAng.sup - pAng.inf)/2;
    arcRad = cRad * pAbs.inf;
    points{3} = ciat.Arc(arcCenter,arcRad,arcAngInf,arcAngSup).polyWrap(dR);
    
    % Segment 4: min amplitude, min phase (bottom right) corner
    arcCenter = cCen*pAbs.inf*exp(1j*pAng.inf);
    arcAngInf = minAng + pi + (pAng.sup - pAng.inf)/2;
    arcAngSup = minAng + 3*pi/2 - rotAng;
    arcRad = cRad * pAbs.inf;
    points{4} = ciat.Arc(arcCenter,arcRad,arcAngInf,arcAngSup).polyWrap(dR);

    % Segment 5: max amplitude, min phase (top right) corner
    arcCenter = cCen*pAbs.sup*exp(1j*pAng.inf);
    arcAngInf = minAng;
    arcAngSup = minAng - pi/2 - rotAng;
    arcRad = cRad * pAbs.sup;
    points{5} = ciat.Arc(arcCenter,arcRad,arcAngInf,arcAngSup).polyWrap(dR);

    % Calculate the four tangent edges
    bndEdge(4,1) = ciat.Edge;
        % outer curve edge
    edgeCenter = (pAbs.sup*(1+cRad)) * exp(1j*pAng.mid);
    % edgeVector = (pAbs.sup*(1+cRad))*exp(1j*(pAng.mid+pi/2));
    edgeVector = pAbs.sup*exp(1j*pAng.inf) - pAbs.sup*exp(1j*pAng.sup);
    bndEdge(1) = ciat.Edge( edgeCenter + [-1 1]*edgeVector);  
        % max phase edge
    edgeCenter = pAbs.mid*exp(1j*pAng.sup) + ...
                 pAbs.mid*cRad*exp(1j*(pAng.sup+pi/2));
    edgeVector = pAbs.width*exp(1j*(pAng.sup+rotAng));
    bndEdge(2) = ciat.Edge( edgeCenter + [-1 1]*edgeVector);  
        % inner curve edge
    edgeCenter = mean([pAbs.inf*exp(1j*pAng.inf),pAbs.inf*exp(1j*pAng.sup)]) + ...
                 pAbs.inf*cRad*exp(1j*(pAng.mid+pi));
    edgeVector = pAbs.inf*exp(1j*pAng.inf) - pAbs.inf*exp(1j*pAng.sup);
    bndEdge(3) = ciat.Edge( edgeCenter + [-1 1]*edgeVector); 
        % min phase edge
    edgeCenter = pAbs.mid*exp(1j*pAng.inf) + ...
                 pAbs.mid*cRad*exp(1j*(pAng.inf-pi/2));
    edgeVector = pAbs.width*exp(1j*(pAng.inf-rotAng));
    bndEdge(4) = ciat.Edge( edgeCenter + [-1 1]*edgeVector);  

    % Adjust the first and last arc wrapping vertices
    points{2}(1) = cap(bndEdge(1) , ciat.Edge(points{2}(1:2)));
    points{2}(end) = cap(bndEdge(2) , ciat.Edge(points{2}(end-1:end)));
    points{3}(1) = cap(bndEdge(2) , ciat.Edge(points{3}(1:2)));
    points{3}(end) = cap(bndEdge(3) , ciat.Edge(points{3}(end-1:end)));
    points{4}(1) = cap(bndEdge(3) , ciat.Edge(points{4}(1:2)));
    points{4}(end) = cap(bndEdge(4) , ciat.Edge(points{4}(end-1:end)));
    points{5}(1) = cap(bndEdge(4) , ciat.Edge(points{5}(1:2)));
    points{5}(end) = cap(bndEdge(1) , ciat.Edge(points{5}(end-1:end)));
    
            
    % Return as column vector
    points = vertcat(points{:});
end


