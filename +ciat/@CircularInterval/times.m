function r = times(obj1,obj2)

% Element-wise multiplication of circular intervals (.* operator)
%
% This function creates the circular intervals representing the 
% element-wise product of two sets of circular intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.CircularInterval class
%   obj2       : array of objects from the ciat.CircularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   circInt = ciat.CircularInterval(0,1) .* ciat.CircularInterval(2,2);
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    % Check input class
    mustBeA(obj1,["ciat.CircularInterval","double"]);
    mustBeA(obj2,["ciat.CircularInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.CircularInterval(obj1, zeros(M1,N1));
    end
    if isa(obj2, 'double')
        obj2 = ciat.CircularInterval(obj2, zeros(M2,N2));
    end 

    % Loop throught the arrays
    r(M,N) = ciat.CircularInterval;
    r.Center = obj1.Center .* obj2.Center;
    r.Radius = abs(obj1.Center) .* obj2.Radius + ...
               abs(obj2.Center) .* obj1.Radius + ...
               obj1.Radius .* obj2.Radius;
end