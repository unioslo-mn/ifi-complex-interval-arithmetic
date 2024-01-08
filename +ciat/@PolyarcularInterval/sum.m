function r = plus(obj1,obj2)

% Addition of polyarcular intervals (+ operator)
%
% This function creates the polyarcular intervals representing the sum
% of two sets of polyarcular intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolyarcularInterval class
%   obj2       : array of objects from the ciat.PolyarcularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
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
    mustBeA(obj1,["ciat.PolyarcularInterval","double"]);
    mustBeA(obj2,["ciat.PolyarcularInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.PolyarcularInterval(obj1);
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolyarcularInterval(obj2);
    end 
            
    % Loop throught the arrays
    r(M,N) = ciat.PolygonalInterval;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            r(M,N) = add( obj1(m1,n1) , obj2(m2,n2) );
            r(M,N).Points = r(M,N).Boundary;
        end
    end
end

%% Function for adding two polygons

function r = add(obj1,obj2)
    % Ensure counter-clockwise order and that v1 and w1
    % being the vertices with smallest y-coordinate 
    % (and smallest x-coordinate in case of ties)
    v = obj1.Points;
    w = obj2.Points;
    
    % Handle exception when one of the inputs is a degenerate interval
    if length(v)==1 || length(w)==1 
        arcs = 
        r = ciat.PolyarcularInterval(arcs);
        return
    end

    % Use Gauss map matching and add arcs
    arcs = 
    r = ciat.PolygonalInterval(arcs);

end