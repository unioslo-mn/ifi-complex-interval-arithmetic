function r = plus(obj1,obj2)

% Addition of circular intervals (+ operator)
%
% This function creates the circular intervals representing the sum
% of two sets of circular intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.CircularInterval class
%   obj2       : array of objects from the ciat.CircularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   circInt = ciat.CircularInterval(0,1) + ciat.CircularInterval(2,2);
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
    mustBeA(obj1,"ciat.CircularInterval");
    mustBeA(obj2,"ciat.CircularInterval");
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.CircularInterval(obj1, 0);
    end
    if isa(obj2, 'double')
        obj2 = ciat.CircularInterval(obj2, 0);
    end 
            
    % Loop throught the arrays
    r(M,N) = ciat.CircularInterval;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            r(M,N).Center = obj1(m1,n1).Center + obj2(m2,n2).Center;
            r(M,N).Radius = obj1(m1,n1).Radius + obj2(m2,n2).Radius;
        end
    end
end