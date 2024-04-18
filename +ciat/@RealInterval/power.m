function r = power(obj, exponent)

% Exponentiation of real intervals (.^ operator)
%
% This function creates the real intervals representing the exponentiation
% of a real interval to a integer (see MATLAB plus function),
% it now works only with integer as exponent.
% _________________________________________________________________________
% USAGE        
%   r = obj .^ exponent
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.RealInterval class
%   exponent       : integer
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = ciat.RealInterval(0,1).^3;
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
    mustBeA(obj,["ciat.RealInterval","double"]);
    mustBeNumeric(exponent);
    mustBeInteger(exponent);
    
    % Get input sizes and check if they can be added
    [M,N] = size(obj);
    
    % Turn scalars to degenerate intervals
    if isa(obj, 'ciat.RealInterval') == 0
        obj = ciat.RealInterval(obj, obj);
    end
    
    % Loop throught the arrays
    r(M,N) = ciat.RealInterval;
    alt1 = obj.Infimum.^exponent;
    alt2 = obj.Supremum.^exponent;
    if mod(exponent,2) == 0
        r.Supremum = max(alt1,alt2);
        r.Infimum = min(alt1,alt2);
        r.Infimum(obj.Infimum <= 0 & obj.Supremum >= 0) = 0;
    else
        r.Infimum = alt1;
        r.Supremum = alt2;
    end
end