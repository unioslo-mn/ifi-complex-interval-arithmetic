function points = sample(obj, nPoints)

% Sample polar interval boundary
%
% This function creates an ordered array of complex values 
% of the requested size representing the boundary points of 
% the sampled polar interval.
% _________________________________________________________________________
% USAGE        
%   r = sample(obj, nPoints)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.PolarInterval class
% _________________________________________________________________________
% OPTIONS
%   nPoints   : total number of points on the boundary (default: 10)
% _________________________________________________________________________
% EXAMPLES
%   points = sample(ciat.PolarInterval(2,3,2,3), 8);
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    arguments
        obj 
        nPoints (1,1)   {mustBeInteger,mustBePositive} = 10;
    end

    [M,N] = size(obj);
    points = cell(M,N);

    for m = 1:M
        for n = 1:N
            % Assemble the point array and remove douplicate points
            points{m,n} = obj(m,n).Center + obj(m,n).Radius * ...
                          exp(1j*linspace(-pi,pi,nPoints)');
        end
    end

end        