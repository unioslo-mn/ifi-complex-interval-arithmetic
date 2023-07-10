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
    
    % Calculate sum
    r(M,N) = ciat.RealInterval;

    bot1 = ceil((obj.Infimum + pi)/(2*pi));
    bot2 = ceil((obj.Supremum + pi)/(2*pi));
    top1 = ceil((obj.Infimum)/(2*pi));
    top2 = ceil((obj.Supremum)/(2*pi));

    r.Infimum(bot2-bot1>=1) = -1;
    r.Infimum(bot2-bot1<1) = min(cos(obj.Infimum(bot2-bot1<1)),cos(obj.Supremum(bot2-bot1<1)));
    r.Supremum(top2-top1>=1) = 1;
    r.Supremum(top2-top1<1) = max(cos(obj.Infimum(top2-top1<1)),cos(obj.Supremum(top2-top1<1)));

end  