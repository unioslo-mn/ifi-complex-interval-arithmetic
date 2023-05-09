function r = union(obj)

% Union of circular intervals
%
% This function creates the circular interval representing the 
% union of a set of circular intervals
% _________________________________________________________________________
% USAGE        
%   r = union(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.CircularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   circInt = union([ciat.CircularInterval(0,1), ...
%                    ciat.CircularInterval(2,2)]);
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    N = length(obj(:));
    assert(N>1)
    r = obj(1);
    for n = 2:N
        r = uniteCircular(r,obj(n));
    end
end

%%

function r = uniteCircular(obj1,obj2)
    
    % Calculate parameters
    cDiff = obj2.Center - obj1.Center;
    cDist = abs(cDiff);
    rDiff = obj2.Radius - obj1.Radius;
    rDist = abs(rDiff);
    
    % Check if the circles are concentric or if one includes the other
    if cDiff == 0 || cDist <= rDist
        if obj1.Radius > obj2.Radius
            r = ciat.CircularInterval(obj1.Center,obj1.Radius);
        else
            r = ciat.CircularInterval(obj2.Center,obj2.Radius);
        end
    else
        rAxis = (cDist + rDiff)/2;
        r = ciat.CircularInterval(rAxis * cDiff / cDist + obj1.Center,...
                                  cDist - abs(rAxis) + obj2.Radius);
    end
end