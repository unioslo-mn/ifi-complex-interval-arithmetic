function r = power(obj, exponent)
% Exponentiation of polar intervals (.^ operator)
%
% This function creates the polar intervals representing the exponentiation
% of a polar interval to a integer (see MATLAB plus function),
% it now works only with integer as exponent.
% _________________________________________________________________________
% USAGE        
%   r = obj .^ exponent
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.PolarInterval class
%   exponent       : integer
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polarInt = ciat.PolarInterval(0,1).^3;
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
    mustBeA(obj, 'ciat.PolarInterval');
    mustBeNumeric(exponent);
    mustBeInteger(exponent);
    
    % Turn scalars to degenerate intervals
    if isa(obj, 'double')
        obj = ciat.PolarInterval(obj, 0);
    end

    r = ciat.PolarInterval(obj.Abs.^exponent, obj.Angle.*exponent);
end
