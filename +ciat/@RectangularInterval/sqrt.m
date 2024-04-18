function r = sqrt(obj)

% Square root of rectangular intervals
%
% This function creates the rectangular intervals representing the square
% root of a rectangular intervals (see MATLAB sqrt function),
% _________________________________________________________________________
% USAGE        
%   r = sqrt(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.RectangularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = sqrt(ciat.RectangularInterval(0,1,2,4))
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________
    
    r_inf = obj.Real.Infimum;
    r_sup = obj.Real.Supremum;
    i_inf = obj.Imag.Infimum;
    i_sup = obj.Imag.Supremum;

    xmin = min(sqrt((r_inf + sqrt(r_inf^2 + i_sup^2))/2), sqrt((r_inf + sqrt(r_inf^2 + i_inf^2))/2));
    xmax = max(sqrt((r_sup + sqrt(r_sup^2 + i_sup^2))/2), sqrt((r_sup + sqrt(r_sup^2 + i_inf^2))/2));

    xmin(i_inf <= 0 & i_sup >= 0) = min(xmin(i_inf <= 0 & i_sup >= 0), sqrt((r_inf(i_inf <= 0 & i_sup >= 0) + sqrt(r_inf(i_inf <= 0 & i_sup >= 0).^2))/2));

    ymin = min(sign(i_inf)*sqrt((-r_inf + sqrt(r_inf^2 + i_inf^2))/2), sign(i_inf)*sqrt((-r_sup + sqrt(r_sup^2 + i_inf^2))/2));
    ymax = max(sign(i_sup)*sqrt((-r_inf + sqrt(r_inf^2 + i_sup^2))/2), sign(i_sup)*sqrt((-r_sup + sqrt(r_sup^2 + i_sup^2))/2));

    r = ciat.RectangularInterval(ciat.RealInterval(xmin, xmax), ciat.RealInterval(ymin, ymax));
end