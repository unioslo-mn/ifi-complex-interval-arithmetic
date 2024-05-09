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
            angles = capGaussMap(obj1.GaussMap,obj2.GaussMap);

            if size(angles,3) > 1
                priMask = ~isnan(angles(:,:,1));
                secMask = ~isnan(angles(:,:,2));
                center = [ center(priMask) ; center(secMask) ];
                radius = [ radius(priMask) ; radius(secMask) ];
                angles = [ angles(priMask) ; angles(secMask) ];
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

%% Function for intersecting the Gauss map intervals

function output = capGaussMap(input1, input2)
    
    % Extract sizes
    [M1,N1] = size(input1);
    [M2,N2] = size(input2);

    % Split angles
    input1 = splitAngle(input1);
    input2 = splitAngle(input2);
    L1 = size(input1,3) > 1;
    L2 = size(input2,3) > 1;

    % If the inputs are two arrays of different orientation
    % form matrices
    if (M1 == M2) && (N1 == N2)
        M = size(input1,1);
        N = size(input1,2);
    elseif (N1 == 1) && (M2 == 1)
        input1 = repmat(input1,1,N2);
        input2 = repmat(input2,M1,1);
        M = M1;
        N = N2;
    elseif (M1 == 1) && (N2 == 1)
        input1 = repmat(input1,M2,1);
        input2 = repmat(input2,1,N1);
        M = M2;
        N = N1;
    else
        error('Incorrect input size for operation.')
    end
  
    % Intersect angles
    primaryCap = cap(input1(:,:,1),input2(:,:,1));
    if ~L1 && L2
        secondaryCap = cap(input1(:,:,1),input2(:,:,2));
    elseif L1 && ~L2
        secondaryCap = cap(input1(:,:,2),input2(:,:,1));
    elseif L1 && L2
        % Merge the adjacent intersections 
        infPositive = primaryCap.Infimum > 0;
        tempCap = cap(input1(:,:,2),input2(:,:,2));
        primaryCap = union( primaryCap - infPositive*2*pi, ...
                            tempCap + ~infPositive*2*pi);

        % There can be only one secondary intersection
        tempCap = cat( 3 , cap(input1(:,:,1),input2(:,:,2)) , ...
                           cap(input1(:,:,2),input2(:,:,1)) );
        secondaryCap(M,N) = ciat.RealInterval;
        secondaryCap( sum(~isnan(tempCap),3)>0 ) = tempCap(~isnan(tempCap));

    end

    % Join primary and secondary caps
    if ~L1 && ~L2
        output = primaryCap;
    else
        moveMask = isnan(primaryCap) & ~isnan(secondaryCap);
        if any(moveMask,'all')
            emptyInterval(1) = ciat.RealInterval;
            primaryCap(moveMask) = secondaryCap(moveMask);
            secondaryCap(moveMask) = emptyInterval;
        end
        
        if all(isnan(secondaryCap),'all')
            output = primaryCap;
        else
            output = cat(3 , primaryCap , secondaryCap);
            lastwarn('Surplus angle intersections found!')
        end
    end
end


%% Function for splitting the angles including the -pi or pi value

function output = splitAngle(input)

    % Check if any elements contain the -pi or pi value
    mask = input.Infimum < -pi | input.Supremum > pi;
    
    % Create a second layer of the array in the 3rd dimension
    if any(mask,'all')
        [M,N] = size(input);
        input2(M,N) = ciat.RealInterval;
        input2(mask) = ciat.RealInterval(...
                             -pi*ones(sum(mask,'all'),1) , ...
                             wrapToPi(input.Supremum(mask)));
        input(mask) = ciat.RealInterval( ...
                            wrapToPi(input.Infimum(mask)) , ...
                            pi*ones(sum(mask,'all'),1));
        output = cat(3,input,input2);
    else
        output = ciat.RealInterval(...
                            wrapToPi(input.Infimum) , ...
                            wrapToPi(input.Supremum) );
    end

end
        