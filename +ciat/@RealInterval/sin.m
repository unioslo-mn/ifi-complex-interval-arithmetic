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
    for m = 1:M
        for n = 1:N
            % Check if the envelope is part of the interval
            bot1 = ceil((obj.Infimum + pi/2)/(2*pi));
            bot2 = ceil((obj.Supremum + pi/2)/(2*pi));
            top1 = ceil((obj.Infimum - pi/2)/(2*pi));
            top2 = ceil((obj.Supremum - pi/2)/(2*pi));
            
            % Calculate sine of infimum and supremum
            sinInf = sin(obj.Infimum);
            sinSup = sin(obj.Supremum);
            
            % Calculate the interval
            if (bot2-bot1)>=1
                r(m,n).Infimum = -1;
            else
                r(m,n).Infimum = min(sinInf,sinSup);
            end
            if (top2-top1)>=1
                r(m,n).Supremum = 1;
            else
                r(m,n).Supremum = max(sinInf,sinSup);
            end
        end
    end
end 