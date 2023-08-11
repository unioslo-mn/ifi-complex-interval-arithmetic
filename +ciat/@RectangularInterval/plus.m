function r = plus(obj1,obj2)

% Addition of rectangular intervals (+ operator)
%
% This function creates the rectangular intervals representing the sum
% of two sets of rectangular intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.RectangularInterval class
%   obj2       : array of objects from the ciat.RectangularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = ciat.RectangularInterval(0,1,2,4) + ...
%             ciat.RectangularInterval(2,3,4,5);
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
    mustBeA(obj1,["ciat.RectangularInterval","double"]);
    mustBeA(obj2,["ciat.RectangularInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.RectangularInterval(obj1);
    end
    if isa(obj2, 'double')
        obj2 = ciat.RectangularInterval(obj2);
    end 
            
    r = ciat.RectangularInterval(obj1.Real + obj2.Real, ...
                                 obj1.Imag + obj2.Imag)

    % If the interval has a probability grid, compute the sum
    if ~isempty(obj1.ProbaGrid) && ~isempty(obj2.ProbaGrid)
        r.ProbaGrid = obj1.ProbaGrid + obj2.ProbaGrid;
        % Fit the grid to the new interval
        r.ProbaGrid = r.ProbaGrid.fitToInterval(r);
    end
end

        