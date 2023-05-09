function r = times(obj1,obj2)

% Element-wise multiplication of rectangular intervals (.* operator)
%
% This function creates the rectangular intervals representing the 
% element-wise product of two sets of rectangular intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.RectangularInterval class
%   obj2       : array of objects from the ciat.RectangularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectangularInt = ciat.RectangularInterval(0,1,2,3) .* ...
%                    ciat.RectangularInterval(2,3,4,5);
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
    mustBeA(obj1,"ciat.RectangularInterval");
    mustBeA(obj2,"ciat.RectangularInterval");
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.RectangularInterval(obj1);
    end
    if isa(obj2, 'double')
        obj2 = ciat.RectangularInterval(obj2);
    end 

    % Loop throught the arrays
    r(M,N) = ciat.RectangularInterval;
    for m = 1:M
        for n = 1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate product
            r(m,n).Real = obj1(m1,n1).Real * obj2(m2,n2).Real - ...
                          obj1(m1,n1).Imag * obj2(m2,n2).Imag;
            r(m,n).Imag = obj1(m1,n1).Real * obj2(m2,n2).Imag + ...
                          obj1(m1,n1).Imag * obj2(m2,n2).Real;
        end
    end
end


        