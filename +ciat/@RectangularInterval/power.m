function r = power(obj, exponent)

% Exponentiation of rectangular intervals
%
% This function creates the rectangular intervals representing the power
% of a rectangular intervals (see MATLAB cos function),
% _________________________________________________________________________
% USAGE        
%   r = pow(obj, exponent)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.RectangularInterval class
%   exponent  : exponent of the power function
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = pow(ciat.RectangularInterval(0,1,2,4), 2)
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
    mustBeA(obj, 'ciat.RectangularInterval');
    mustBeNumeric(exponent);
    mustBeInteger(exponent);
    if exponent ~= 2
        error('Only implemented for exponent = 2');
    end
    
    r = ciat.RectangularInterval(obj.Real.^2 - obj.Imag.^2, 2*obj.Real.*obj.Imag);
end