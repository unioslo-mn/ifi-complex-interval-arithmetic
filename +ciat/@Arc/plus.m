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

    % If the second objects is an edge, forward to the edge plus function
    if isa(obj2,'ciat.Edge')
        r = obj2 + obj1;
        return
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
            angles = ciat.Arc.capGaussMap(obj1.GaussMap,obj2.GaussMap);

            % Handle cases of two disconnected intersections
            if size(angles,3) > 1
                priAng = angles(:,:,1);
                secAng = angles(:,:,2);
                priMask = ~isnan(priAng);
                secMask = ~isnan(secAng);
                center = [ center(priMask) ; center(secMask) ];
                radius = [ radius(priMask) ; radius(secMask) ];
                angles = [ priAng(priMask) ; secAng(secMask) ];
            end

        case 'double'
            center = obj1.Center + obj2;
            radius = obj1.Radius;
            angles = obj1.ArcAngle;
    end

    % Create new object        
    [M,N] = size(angles);
    r(M,N) = ciat.Arc;
    mask = ~isnan(angles);
    if any(mask)
        % r(mask) = ciat.Arc(center(mask),radius(mask),angles(mask));
        r.Center(mask) = center(mask);
        r.Radius(mask) = radius(mask);
        r.ArcAngle(mask) = angles(mask);
    end
end

