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
            
            %outPoints = (inRadius + dR) .* exp(1j*angles(2:end)) + inCenter;
            outPoints = (inRadius/cos(2*pi/(2*cntPoints))) .* ...
                            exp(1j*angles(2:end)) + inCenter;
            
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
                %rH = (maxAbs + dR) * exp(1j * angs); % outer arc
                rH = (maxAbs/cos((maxAngle-minAngle)/(2*(cntPoints-1)))) * ...
                            exp(1j * angs); % outer arc
                pH = minAbs * exp( 1j * angs(end));  % inner corner 2

                % Checks and corrections in case a half or full circle is made 
                maxAngle(fullCircle) = 2*pi;
                minAngle(fullCircle) = 0;
                pL(halfCircle) = [];
                pH(halfCircle) = [];

                % Compile points
                if inAbs.width == 0
                    outPoints = rH.';
                elseif inAngle.width == 0
                    outPoints = [pL,maxAbs * exp(1j*maxAngle)].';
                else
                    outPoints = [pL,rH, pH].';                    
                end


            else
                if isa(inObj2,'ciat.CircularInterval')
                    outPoints = ciat.PolygonalInterval.castPolarTimesCircular(...
                                                        inObj,inObj2,dR);
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
