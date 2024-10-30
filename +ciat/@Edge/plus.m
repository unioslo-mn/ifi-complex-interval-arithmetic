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
    
    % Check if one of the objects is an edge
    assert(isa(obj1,'ciat.Edge') || isa(obj2,'ciat.Edge'))
    

    % If the first object is a double, flip the objects
    if isa(obj1,'double')
        objTemp = obj1;
        obj1 = obj2;
        obj2 = objTemp;
    end

    % If the inputs are two arrays of different orientation
    % form matrices
    if (M1 == M2) && (N1 == N2)
        M = size(obj1,1);
        N = size(obj1,2);
    elseif (N1 == 1) && (M2 == 1)
        obj1 = repmat(obj1,1,N2);
        obj2 = repmat(obj2,M1,1);
        M = M1;
        N = N2;
    elseif (M1 == 1) && (N2 == 1)
        obj1 = repmat(obj1,M2,1);
        obj2 = repmat(obj2,1,N1);
        M = M2;
        N = N1;
    else
        error('Incorrect input size for operation.')
    end
        
    % Calculate parameters
    switch class(obj2)
        case 'ciat.Edge'
            mask = (obj1.GaussMap == obj2.GaussMap);
            if any(mask,'all')
                p1 = obj1.Startpoint + obj2.Startpoint;
                p2 = obj1.Endpoint + obj2.Endpoint;
            end
        case 'double'
            p1 = obj1.Startpoint + obj2;
            p2 = obj1.Endpoint + obj2;
            mask = ~isnan(p1) & ~isnan(p2);
        
        case 'ciat.Arc'
            arcGauss = obj2.GaussMap;
            edgeGauss = obj1.GaussMap.Midpoint;
            mask = isin( arcGauss , edgeGauss ) | ...
                   isin( arcGauss , edgeGauss - sign(edgeGauss)*2*pi);
            offset = obj2.Center + obj2.Radius .* exp(1i*edgeGauss);
            p1 = obj1.Startpoint + offset;
            p2 = obj1.Endpoint + offset;
    end

    % Create new object        
    r(M,N) = ciat.Edge;
    if any(mask,'all')
        r.Startpoint(mask) = p1(mask);
        r.Endpoint(mask) = p2(mask);
    end
end

        