function r = times(obj1,obj2)
   
% Element-wise multiplication of polyarcular intervals (.* operator)
%
% This function creates the polyarcular intervals representing the 
% element-wise product of two sets of polyarcular intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
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
        obj1 = ciat.PolyarcularInterval(obj1, 0);
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolyarcularInterval(obj2, 0);
    end 

    % Loop throught the arrays
    r(M,N) = ciat.PolyarcularInterval;
    for m = 1:M
        for n = 1:N
            % Calculate product
            r(M,N) = multiply( obj1(m1,n1) , obj2(m2,n2) );
        end
    end
end

%%
function r = multiply(obj1,obj2)
    % Ensure counter-clockwise order and that v1 and w1
    % being the vertices with smallest y-coordinate 
    % (and smallest x-coordinate in case of ties)
    v = obj1.Arcs;
    w = obj2.Arcs;
    
    % Handle exception when one of the inputs is a degenerate interval
    if length(v)==1 || length(w)==1 
        arcs = 
        r = ciat.PolyarcularInterval(arcs);
        return
    end

    % Use log Gauss map matching and multiply arcs
    % Use segmentProduct(v(index),w(index)) function to get the new arc
    r = ciat.PolygonalInterval(arcs);

end