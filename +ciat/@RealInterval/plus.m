function r = plus(obj1,obj2)

% Addition of real intervals (+ operator)
%
% This function creates the real intervals representing the sum
% of two sets of real intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.RealInterval class
%   obj2       : array of objects from the ciat.RealInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = ciat.RealInterval(0,1) + ciat.RealInterval(2,2);
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
    mustBeA(obj1,["ciat.RealInterval","double"]);
    mustBeA(obj2,["ciat.RealInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    % M = max([M1,M2]);
    % N = max([N1,N2]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'ciat.RealInterval') == 0
                obj1 = ciat.RealInterval(obj1, obj1);
    end
    if isa(obj2, 'ciat.RealInterval') == 0
        obj2 = ciat.RealInterval(obj2, obj2);
    end 
    
    r = ciat.RealInterval(obj1.Infimum + obj2.Infimum, obj1.Supremum + obj2.Supremum);
end