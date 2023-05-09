function r = mtimes(obj1,obj2)

% Matrix multiplication of polygonal intervals (* operator)
%
% This function creates the polygonal intervals representing the matrix 
% product of two sets of polygonal intervals (see MATLAB mtimes function),
% which it achieves by using the times (.*) and plus (+) functions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 * obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolygonalInterval class
%   obj2       : array of objects from the ciat.PolygonalInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polycInt = ciat.PolygonalInterval([0,1,1i]) * ...
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
    mustBeA(obj1,"ciat.PolygonalInterval");
    mustBeA(obj2,"ciat.PolygonalInterval");
   
    % Get input sizes and check if they can be multiplied
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(N1 == M2 || (M1 == 1 && N1 == 1) || (M2 == 1 && N2 == 1) )
    if N1 == M2
        M = M1;
        N = N2;
    else
       M = max(M1,M2);
       N = max(N1,N2);
    end
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.PolygonalInterval(obj1, 0);
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolygonalInterval(obj2, 0);
    end 

    % Loop through the arrays
    r(M,N) = ciat.PolygonalInterval;
    if N1 == M2
        for m = 1:M
            for n = 1:N
                % Calculate product using the times function
                r(M,N) = sum( obj1(m,:) .* obj2(:,n)' );
            end
        end
    else
        % Calculate product using the times function
        r = obj1 .* obj2;
    end
end 