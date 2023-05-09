function points = backtrack(obj,trackAngle)

% Backtrack polygonal intervals
%
% This function finds the point on the boundary of the given 
% polygonal intervals that has the selected normal angle.
% Polygonal vertices have an interval of normal angles ranging
% between the normal angles of the adjaces edges. It can happen
% that the selected angle is excatly the normal angle of an edge
% in which case both adjacent vertices will include that angle.
% In case of ambiguity the one first in counter-clockwise order
% is selected.
% _________________________________________________________________________
% USAGE        
%   r = backtrack(obj, trackAngle)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.CircularInterval class
%   trackAngle: angle selected for backtracking
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   points = backtrack([ciat.PolygonalInterval(complex(rand(3),rand(3))), ...
%                       ciat.PolygonalInterval(complex(rand(3),rand(3)))],...
%                       pi/3);
%_________________________________________________________________________
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
        trackAngle      (1,1)   {mustBeNumeric} = 0
    end

    points = zeros(1, length(obj));
    eps10 = 10*eps;

    for i = 1:length(obj) % loop over all obj
        v = obj(i).Points;
        
        % init prev angle
        prev_angle = wrapTo2Pi( angle( v(1) - v(end) )- pi/2); 
        
        % just assume it is the final point to avoid extra condition later
        points(i) = v(end); 
        
        for j = 1 : length(v)-1
            current_angle = wrapTo2Pi( angle( v(j+1) - v(j) )- pi/2);

            if (current_angle - prev_angle) < 0 % we crossed 0
                condition1 = (trackAngle >= prev_angle       -eps10) && ...
                                  (trackAngle <= 2*pi+current_angle+eps10);
                condition2 = (trackAngle >= prev_angle -2*pi -eps10) && ...
                                  (trackAngle <=      current_angle+eps10);
                condition = condition1 || condition2;
            else
                condition = (trackAngle >= prev_angle-eps10) && ...
                                 (trackAngle <= current_angle+eps10);
            end

            if condition == true
                points(i) = v(j);
                break;
            end
            prev_angle = current_angle;
        end
    end 

end

