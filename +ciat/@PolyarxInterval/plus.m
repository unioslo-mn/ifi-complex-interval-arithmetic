function r = plus(obj1,obj2)

% Addition of convex polyarcular intervals (+ operator)
%
% This function creates the convex polyarcular intervals representing the 
% sum of two sets of convex polyarcular intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolyarxInterval class
%   obj2       : array of objects from the ciat.PolyarxInterval class
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
    mustBeA(obj1,["ciat.PolyarxInterval","double"]);
    mustBeA(obj2,["ciat.PolyarxInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % % Turn scalars to degenerate intervals
    % if isa(obj1, 'double')
    %     obj1 = ciat.PolyarxInterval(obj1);
    % end
    % if isa(obj2, 'double')
    %     obj2 = ciat.PolyarxInterval(obj2);
    % end 
            
    % Loop throught the arrays
    r(M,N) = ciat.PolyarxInterval;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            if isa(obj1,'double')
                r(m,n) = obj2(m2,n2);
                r(m,n).Arx(:,1) = r(m,n).Arx(:,1) + real(obj1(m1,n1));
                r(m,n).Arx(:,2) = r(m,n).Arx(:,2) + imag(obj1(m1,n1));
            elseif isa(obj2,'double')
                r(m,n) = obj1(m1,n1);
                r(m,n).Arx(:,1) = r(m,n).Arx(:,1) + real(obj2(m2,n2));
                r(m,n).Arx(:,2) = r(m,n).Arx(:,2) + imag(obj2(m2,n2));
            else
                r(m,n) = add( obj1(m1,n1) , obj2(m2,n2) );
            end
        end
    end
end

%% Function for adding two polygons

function r = add(obj1,obj2)

    % Extrac arx
    arx1 = obj1.Arx;
    arx2 = obj2.Arx;

    % Extract parameters
    cen1 = complex(arx1(:,1),arx1(:,2));
    rad1 = arx1(:,3);
    ang1 = arx1(:,4);
    cen2 = complex(arx2(:,1),arx2(:,2));
    rad2 = arx2(:,3);
    ang2 = arx2(:,4);
    
    % Taken from the polygonal plus function
    N1 = size(ang1,1);
    N2 = size(ang2,1);
    N3 = N1 + N2;
    cen3 = zeros(N3,1);
    rad3 = zeros(N3,1);
    ang3 = zeros(N3,1);
    n1 = 1;
    n2 = 1;
    n3 = 0;
    eps10 = eps*10;
    while (n1 <= N1) && (n2 <= N2) % continue finding more points
        
        n3 = n3 + 1;
    
        % Sum arcs
        cen3(n3) = cen1(n1) + cen2(n2);
        rad3(n3) = rad1(n1) + rad2(n2);
        ang3(n3) = min(ang1(n1),ang2(n2));
        
        % Increment index
        if  ang1(n1) < ang2(n2) + eps10 
            n1 = n1 + 1;
        elseif ang1(n1) > ang2(n2) + eps10
            n2 = n2 + 1;
        else
            n1 = n1 + 1;
            n2 = n2 + 1;
        end            

        if n3 > 1 && ang3(n3) <= ang3(n3-1) + eps10
            n3 = n3 - 1;
        end
    end
    cen3 = cen3(1:n3,:);
    rad3 = rad3(1:n3,:);
    ang3 = ang3(1:n3,:);
       
    % Generate polyarx
    arxSum = [real(cen3),imag(cen3),rad3,ang3];
    r = ciat.PolyarxInterval(arxSum);

end