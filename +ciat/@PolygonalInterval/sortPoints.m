function points = sortPoints(obj)

% Sort the vertex points of polygonal intervals
%
% This function finds the convex hull of the input polygonal interval
% and sorts the points in counter-clockwise order. This is a utility
% function used in the set method of the Points property. It works only
% on a single interval object, not an array.
% _________________________________________________________________________
% USAGE        
%   r = sortPoints(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.CircularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%_________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    
    % Extract points from the object
    points = obj.Boundary;
    
    % Remove non-unique points
    points = uniquetol([real(points), imag(points)] ,...
                      obj.Tolerance,'ByRows',true, 'DataScale', 1);

    % Find convex hull
    if (length(points) > 3) && ~(isColinear(points(:,1) , points(:,2)))
        points = points( convhull(points(:,1) , points(:,2)) , : );
    end
    points = points(:,1) + 1j*points(:,2);
    
    % Sort points to form a counter-clockwise boundary and assign to object
    points = sortCounterClockwise(points);
    
end

%% Utility functions

function value = isColinear(x,y)

    % See if (v2 - v1) is parallel to all of (vi - v1), for all i. If so
    % all points lie on the same line
    almost_one = (1-5*eps);
    n_points = length(x);
    
    value = true;
    % if loop above does trigger a rewrite below, either:
    %   1) there are less than three points _OR_
    %   2) ALL points are colinear
    
    if n_points >= 3
        v2  = [x(2)-x(1); y(2)-y(1)];
        nv2 = norm(v2);
        v2 = v2./nv2;
        
        for i = 3:n_points
            vi  = [x(i)-x(1); y(i)-y(1)];
            nvi = norm(vi);
            if dot(v2, vi)/nvi < almost_one
                value = false;
                break;
            end
        end  
    end 
end

function points = sortCounterClockwise(points)
    % Sort points counter-clockwise
    center = mean(points);
    angles = ciat.wrapTo2Pi( angle(points - center) );
    [~, sortIdx] = sort(angles,'ascend');
    points = points(sortIdx);
    
    % Shift starting point to closest to zero
    [minIm, minIdx] = min(imag(points));
    points = circshift(points, -minIdx+1);

    % In case of a tie
    tieIdx = ( imag(points) == minIm );
    [~, tieSortIdx] = sort(real(points(tieIdx)),'ascend');
    points(tieIdx) = points(tieSortIdx);
end