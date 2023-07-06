function r = cosh(obj)

% Hyperbolic cosine of real intervals
%
% This function creates the real intervals representing the hyperbolic cosine 
% of the given set of real intervals (see MATLAB cos function).
% _________________________________________________________________________
% USAGE        
%   r = cosh    (obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj        : array of objects from the ciat.RealInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = cosh(ciat.RealInterval(0,1));
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    % Get input size
    [M,N] = size(obj);
    
    % Calculate sum
    r(M,N) = ciat.RealInterval;
    % cosh function is decreasing on negative numbers and increasing on positive numbers
    r.Supremum = max(cosh(obj.Supremum),cosh(obj.Infimum));
    r.Infimum = min(cosh(obj.Supremum),cosh(obj.Infimum));

    r.Infimum(obj.Infimum <= 0 & obj.Supremum >= 0) = 1;
end 