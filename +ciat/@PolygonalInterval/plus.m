function r = plus(obj1,obj2)

% Addition of polygonal intervals (+ operator)
%
% This function creates the polygonal intervals representing the sum
% of two sets of polygonal intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolygonalInterval class
%   obj2       : array of objects from the ciat.PolygonalInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polycInt = ciat.PolygonalInterval([0,1,1i]) + ...
%              ciat.PolygonalInterval([0,-1,-1i]);
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
    mustBeA(obj1,["ciat.PolygonalInterval","double"]);
    mustBeA(obj2,["ciat.PolygonalInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.PolygonalInterval(obj1);
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolygonalInterval(obj2);
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
            r(m,n) = add( obj1(m1,n1) , obj2(m2,n2) );
            r(m,n).Points = r(m,n).Boundary;
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
        points = reshape(v + w.' ,[],1);
        r = ciat.PolygonalInterval(points);
        return
    end

    i = 1; j = 1; 
    eps10 = 10*eps; % needed so avoid skipping vertices due to numerical precision

    I = length(obj1.Points) + 1; 
    v = [v; v(1); v(2)];
    v_arg = ciat.wrapTo2Pi( angle( v(2:end) - v(1:end-1) ));
    v_arg(end) = v_arg(end) + 2*pi; % otherwise it wraps around

    J = length(obj2.Points) + 1; 
    w = [w; w(1); w(2)];
    w_arg = ciat.wrapTo2Pi( angle( w(2:end) - w(1:end-1) ));
    w_arg(end) = w_arg(end) + 2*pi; % otherwise it wraps around

    p = zeros( I + J, 1);
    n = 0;

    % v_arg(i) = angle(v_i+1 - v_i), which is why we also
    % repeat when i=I and j=J (as opposed to in DeBerg2008). 
    % The loop breaks in the following turn.
    while (i <= I) && (j <= J) % continue finding more points
        n = n+1;
        p(n) = v(i) + w(j); % add what we found

        if     v_arg(i) < w_arg(j) + eps10 
            i = i + 1;
        elseif v_arg(i) > w_arg(j) + eps10
            j = j + 1;
        else
            i = i + 1;
            j = j + 1;
        end            
    end

    r = ciat.PolygonalInterval(p(1:n-1));

end

