function r = cos(obj)

% Cosine of real intervals
%
% This function creates the real intervals representing the cosine 
% of the given set of real intervals (see MATLAB cos function).
% _________________________________________________________________________
% USAGE        
%   r = cos(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj        : array of objects from the ciat.RealInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = cos(ciat.RealInterval(0,1));
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
    
    % Initialize value
    r(M,N) = ciat.RealInterval;

    % Conditions of including the envelope
    btm1 = ceil((obj.Infimum + pi)/(2*pi));
    btm2 = ceil((obj.Supremum + pi)/(2*pi));
    btmEnv = btm2-btm1>=1;
    top1 = ceil((obj.Infimum)/(2*pi));
    top2 = ceil((obj.Supremum)/(2*pi));
    topEnv = top2-top1>=1;
    


    r.Infimum(btmEnv) = -1;
    r.Infimum(~btmEnv) = min(cos(wrapToPi(obj.Infimum(~btmEnv))),...
                             cos(wrapToPi(obj.Supremum(~btmEnv))));
    r.Supremum(topEnv) = 1;
    r.Supremum(~topEnv) = max(cos(wrapToPi(obj.Infimum(~topEnv))),...
                              cos(wrapToPi(obj.Supremum(~topEnv))));

end  