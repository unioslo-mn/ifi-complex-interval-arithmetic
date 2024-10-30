function r = sqrt(obj)

% Square root of polar intervals
%
% This function creates the polar intervals representing the square
% root of a polar intervals (see MATLAB sqrt function),
% _________________________________________________________________________
% USAGE        
%   r = sqrt(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.PolarInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = sqrt(ciat.PolarInterval(0,1,2,4))
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    r = ciat.PolarInterval(sqrt(obj.Abs), obj.Angle./2);
end
