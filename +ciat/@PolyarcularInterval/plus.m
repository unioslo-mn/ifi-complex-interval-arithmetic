function r = plus(obj1,obj2)

% Addition of polyarcular intervals (+ operator)
%
% This function creates the polyarcular intervals representing the sum
% of two sets of polyarcular intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
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
    
    % If the first object is a double make that the second
    if isa(obj1, 'double')
        objTemp = obj1;
        obj1 = obj2;
        obj2 = objTemp;
    end
            
    % Loop throught the arrays
    r(M,N) = ciat.PolyarcularInterval;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            if isa(obj2,'double')
                % r(M,N) = obj1(m1,n1);
                % r.ArcStorage{M,N}.Center = r.ArcStorage{M,N}.Center + ...
                %                             obj2(m2,n2);
                r(m,n) = ciat.PolyarcularInterval(obj1(m1,n1).DefArcs{:} + ...
                                                  obj2(m2,n2));
            else
                if obj1(m1,n1).isconvex && obj2(m2,n2).isconvex
                    r(m,n) = ciat.PolyarcularInterval.plusConvex( ...
                                            obj1(m1,n1) , obj2(m2,n2) );
                else
                    r(m,n) = ciat.PolyarcularInterval.plusConcave( ...
                                            obj1(m1,n1) , obj2(m2,n2) );
                end
            end
        end
    end
end
