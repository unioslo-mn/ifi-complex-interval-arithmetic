function points = sample(obj, nPoints)

% Sample polar interval boundary
%
% This function creates an ordered array of complex values 
% of the requested size representing the boundary points of 
% the sampled polar interval.
% _________________________________________________________________________
% USAGE        
%   r = sample(obj, nPoints)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.PolarInterval class
% _________________________________________________________________________
% OPTIONS
%   nPoints   : total number of points on the boundary (default: 10)
% _________________________________________________________________________
% EXAMPLES
%   points = sample(ciat.PolarInterval(2,3,2,3), 8);
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
        obj 
        nPoints = 10;
    end

    % Extract parameters
    maxAbs = obj.Abs.Supremum;
    minAbs = obj.Abs.Infimum;
    maxAng = obj.Angle.Supremum;
    minAng = obj.Angle.Infimum;                     
    
    % Calculate length of each side
    lenArcMin = minAbs * (maxAng-minAng);
    lenArcMax = maxAbs * (maxAng-minAng);
    lenEdge = maxAbs - minAbs;
    lenTotal = lenArcMin + lenArcMax + 2*lenEdge;
    
    % Divide the points between the sides
    cntEdge = floor(nPoints * lenEdge / lenTotal);
    cntArcMax = floor(nPoints * lenArcMax / lenTotal);
    cntArcMin = nPoints - cntArcMax - 2*cntEdge;   
    

    % Generate points in counter-clockwise order
    pntEdgeMin = linspace(minAbs, maxAbs, cntEdge+1) * exp(1j*minAng);
    pntArcMax = (maxAbs) * exp(1j * linspace(minAng, maxAng, cntArcMax+1));
    pntEdgeMax = linspace(maxAbs, minAbs, cntEdge+1) * exp(1j*maxAng);
    pntArcMin = (minAbs) * exp(1j * linspace(maxAng, minAng, cntArcMin+1));
    
    % Assemble the point array and remove douplicate points
    points = [pntEdgeMin(1:end-1), ...
              pntArcMax(1:end-1), ...
              pntEdgeMax(1:end-1),...
              pntArcMin(1:end-1)].';

end        