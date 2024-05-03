function outObj = cast(inObj)

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
            inCenter = [inObj.Center];
            inRadius = [inObj.Radius];
            outArcs = ciat.Arc(inCenter,inRadius,ciat.RealInterval(-pi,pi));
            
        case 'ciat.PolarInterval'
            inAbs = [inObj.Abs];
            inAngle = [inObj.Angle];
            maxAbs = [inAbs.Supremum];
            minAbs = [inAbs.Infimum];
            outArcs(2,1) = ciat.Arc;
            outArcs(1) = ciat.Arc(0,-minAbs,inAngle+pi);
            outArcs(2) = ciat.Arc(0,maxAbs,inAngle);
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

