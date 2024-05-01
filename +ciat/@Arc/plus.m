function r = plus(obj1,obj2)

% Addition of an arc segment with an arc, edge or vertex (+ operator)
%
% This function creates the arc representing the sum of two arcs, an 
% arc and an edge or an arc and a vertex (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.Arc class
%   obj2       : array of objects from the ciat.Arc, ciat.Edge or double class
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
    mustBeA(obj1,["ciat.Arc","ciat.Edge","double"]);
    mustBeA(obj2,["ciat.Arc","ciat.Edge","double"]);
        
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    
    % Check if one of the objects is an arc
    assert(isa(obj1,'ciat.Arc') || isa(obj2,'ciat.Arc'))

    % If any of the objects is an edge, forward to the edge plus function
    if isa(obj1,'ciat.Edge') || isa(obj2,'ciat.Edge')
        r = ciat.Edge.plus(obj1,obj2);
    end

    % If the first object is a double, flip the objects
    if isa(obj1,'double')
        objTemp = obj1;
        obj1 = obj2;
        obj2 = objTemp;
    end
        
    % Calculate parameters
    switch class(obj2)

        case 'ciat.Arc'
            center = obj1.Center + obj2.Center;
            radius = obj1.Radius + obj2.Radius;
            
            % Intersect Gauss-maps
            G1 = obj1.GaussMap;
            G2 = obj2.GaussMap;
            if size(G1,3) > 1 && size(G2,3) > 1
                angleCap = cat(3,cap(G1(:,:,1),G2(:,:,1)), ...
                               cap(G1(:,:,1),G2(:,:,2)), ...
                               cap(G1(:,:,2),G2(:,:,1)), ...
                               cap(G1(:,:,2),G2(:,:,2)) );
                angleCap = sort(angleCap,3);
                nanLayer = all(isnan(angleCap),[1,2]);
                nanLayer = find(nanLayer,1,'first');
                angles = angleCap(:,:,1:nanLayer-1);

            elseif size(G1,3) == 1 && size(G2,3) > 1
                angleCap = cat(3,cap(G1(:,:,1),G2(:,:,1)), ...
                               cap(G1(:,:,1),G2(:,:,2)) );
                angleCap = sort(angleCap,3);
                nanLayer = all(isnan(:,:,2));
                angles = angleCap(:,:,1:nanLayer+1);
            elseif size(G1,3) > 1 && size(G2,3) == 1
                angleCap = cat(3,cap(G1(:,:,1),G2(:,:,1)), ...
                               cap(G1(:,:,2),G2(:,:,1)) );
                angleCap = sort(angleCap,3);
                nanLayer = all(isnan(:,:,2));
                angles = angleCap(:,:,1:nanLayer+1);
            else
                angles = cap(G1,G2);
            end

        case 'double'
            center = obj1.Center + obj2;
            radius = obj1.Radius;
            angles = obj1.ArcAngle;
    end

    % Create new object        
    r(M1,N1) = ciat.Arc;
    mask = ~isnan(angles);
    if any(mask)
        % r(mask) = ciat.Arc(center(mask),radius(mask),angles(mask));
        r.Center(mask) = center(mask);
        r.Radius(mask) = radius(mask);
        r.ArcAngle(mask) = angles(mask);
    end
end

        