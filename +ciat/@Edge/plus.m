function r = plus(obj1,obj2)

% Addition of an edge segment with an edge, arc or vertex (+ operator)
%
% This function creates the edge representing the sum of two edges, an 
% edge and an arc or an edge and a vertex (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.Edge class
%   obj2       : array of objects from the ciat.Edge, ciat.Arc or double class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   
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
    mustBeA(obj1,["ciat.Edge","ciat.Arc","double"]);
    mustBeA(obj2,["ciat.Edge","ciat.Arc","double"]);
        
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    
    % Check if one of the objects is an arc
    assert(isa(obj1,'ciat.Edge') || isa(obj2,'ciat.Edge'))

    % If the first object is a double, flip the objects
    if isa(obj1,'double')
        objTemp = obj1;
        obj1 = obj2;
        obj2 = objTemp;
    end
        
    % Calculate parameters
    switch class(obj2)
        case 'ciat.Edge'
            mask = (obj1.GaussMap == obj2.GaussMap);
            p1 = obj1.Endpoints(1) + obj2.Endpoints(2);
            p2 = obj1.Endpoints(1) + obj2.Endpoints(2);

        case 'double'
            p1 = obj1.Startpoint + obj2;
            p2 = obj1.Endpoint + obj2;
            mask = ones(size(p1,1),size(p1,2));
        
        case 'ciat.Arc'
            mask = isin(obj2.ArcAngles,obj1.GaussMap.Midpoint);
            offset = obj2.Radius * exp(1i*obj1.GaussMap.Midpoint);
            p1 = obj1.Startpoint + offset;
            p2 = obj1.Endpoint + offset;
    end

    % Create new object        
    r(M1,N1) = ciat.Edge;
    if any(mask)
        r(mask) = ciat.Edge(p1(mask),p2(mask));
    end
end

        