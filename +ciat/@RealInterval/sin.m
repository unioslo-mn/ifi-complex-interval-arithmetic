function r = sin(obj)

% Sine of real intervals
%
% This function creates the real intervals representing the sine 
% of the given set of real intervals (see MATLAB cos function).
% _________________________________________________________________________
% USAGE        
%   r = sin(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj        : array of objects from the ciat.RealInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = sin(ciat.RealInterval(0,1));
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
    mustBeA(obj,"ciat.RealInterval");
    
    % Get input size
    [M,N] = size(obj);
    
    % Calculate sum
    r(M,N) = ciat.RealInterval;

    % With vectorized operations
    btm1 = ceil((obj.Infimum + pi/2)/(2*pi));
    btm2 = ceil((obj.Supremum + pi/2)/(2*pi));
    btmEnv = btm2-btm1>=1;
    top1 = ceil((obj.Infimum - pi/2)/(2*pi));
    top2 = ceil((obj.Supremum - pi/2)/(2*pi));
    topEnv = top2-top1>=1;

    r.Infimum(btmEnv) = -1;
    r.Infimum(~btmEnv) = min(sin(wrapToPi(obj.Infimum(~btmEnv))), ...
                             sin(wrapToPi(obj.Supremum(~btmEnv))));
    r.Supremum(topEnv) = 1;
    r.Supremum(~topEnv) = max(wrapToPi(sin(obj.Infimum(~topEnv))), ...
                              wrapToPi(sin(obj.Supremum(~topEnv))));
    
end 