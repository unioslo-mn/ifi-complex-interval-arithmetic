function points_backtracked = backtrack(obj,trackAngle)

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
%   points_backtracked = backtrack([ciat.PolygonalInterval(complex(rand(3),rand(3))), ...
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
    
    trackAngle = wrapTo2Pi(trackAngle);
    points_backtracked = zeros(1, length(obj));

    for i = 1:length(obj)  % loop over all polygons to find a backtracked point
        v = obj(i).Points;       % vertices of current polygon
        v = unique(v, 'stable'); % needed in cases where vertices are duplicated (may be appropriate to warn user, although it rarely matters)
        
        prev_angle = wrapTo2Pi( angle(v(1) - v(end)) - pi/2 ); % init prev angle
        points_backtracked(i) = v(end); % "assume" it is the final point to avoid an extra condition later (***)
        
        for j = 1 : length(v)-1 % loop over every vertex except the last (***)
            current_angle = wrapTo2Pi( angle(v(j+1) - v(j)) - pi/2 );

            if (current_angle - prev_angle) < 0 % we crossed 0 (x-axis)
                condition1 = (trackAngle >= prev_angle) && (trackAngle <=          2*pi);
                condition2 = (trackAngle >=          0) && (trackAngle <= current_angle);
                condition = condition1 || condition2;
            else
                condition = (trackAngle >= prev_angle) && (trackAngle <= current_angle);
            end

            if condition == true
                points_backtracked(i) = v(j);
                break;
            end
            
            prev_angle = current_angle; % we go to next vertex (shares edge, so this statement holds)
        
        end 
    end 
end

