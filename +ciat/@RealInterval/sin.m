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
    bot1 = ceil((obj.Infimum + pi/2)/(2*pi));
    bot2 = ceil((obj.Supremum + pi/2)/(2*pi));
    top1 = ceil((obj.Infimum - pi/2)/(2*pi));
    top2 = ceil((obj.Supremum - pi/2)/(2*pi));

    r(bot2-bot1>=1).Infimum = -1;
    r(bot2-bot1<1).Infimum = min(sin(obj(bot2-bot1<1).Infimum),sin(obj(bot2-bot1<1).Supremum));
    r(top2-top1>=1).Supremum = 1;
    r(top2-top1<1).Supremum = max(sin(obj(top2-top1<1).Infimum),sin(obj(top2-top1<1).Supremum));
end 