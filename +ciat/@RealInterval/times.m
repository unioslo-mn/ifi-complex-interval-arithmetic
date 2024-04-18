function r = times(obj1,obj2)

% Element-wise multiplication of real intervals (.* operator)
%
% This function creates the real intervals representing the 
% element-wise product of two sets of real intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.RealInterval class
%   obj2       : array of objects from the ciat.RealInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = ciat.RealInterval(0,1) .* ciat.RealInterval(2,3);
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
    mustBeA(obj1,["ciat.RealInterval","double"]);
    mustBeA(obj2,["ciat.RealInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'ciat.RealInterval') == 0
           obj1 = ciat.RealInterval(obj1, obj1);
    end
    if isa(obj2, 'ciat.RealInterval') == 0
           obj2 = ciat.RealInterval(obj2, obj2);
    end 
    
    alt1 = obj1.Infimum .* obj2.Infimum; % LL
    alt2 = obj1.Infimum .* obj2.Supremum; % LU
    alt3 = obj1.Supremum .* obj2.Infimum; % UL
    alt4 = obj1.Supremum .* obj2.Supremum; % UU

    % Take element-wise minimum and maximum of the four alternatives along every dimension
    r = ciat.RealInterval(min(min(alt1,alt2),min(alt3,alt4)),max(max(alt1,alt2),max(alt3,alt4)));
end

