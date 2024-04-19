function r = times(obj1,obj2)
   
% Element-wise multiplication of polygonal intervals (.* operator)
%
% This function creates the polygonal intervals representing the 
% element-wise product of two sets of polygonal intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolygonalInterval class
%   obj2       : array of objects from the ciat.PolygonalInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   circInt = ciat.PolygonalInterval([0,1,1i]) .* ...
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
        obj1 = ciat.PolygonalInterval(obj1, 0);
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolygonalInterval(obj2, 0);
    end 

    % Loop throught the arrays
    r(M,N) = ciat.PolygonalInterval;
    for m = 1:M
        for n = 1:N
            % Calculate product
            r(m,n).Points = reshape(obj1.Points * obj2.Points.' ,[],1);
            r(m,n).Points = r(m,n).Boundary;
        end
    end

    % If the interval has a probability grid, compute the product
    if ~isempty(obj1.ProbaGrid) && ~isempty(obj2.ProbaGrid)
        r.ProbaGrid = obj1.ProbaGrid .* obj2.ProbaGrid;
        % Fit the grid to the new interval
        r.ProbaGrid = r.ProbaGrid.fitToInterval(r);
    end
end