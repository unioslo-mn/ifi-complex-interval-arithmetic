function r = times(obj1,obj2)
   
% Element-wise multiplication of polar intervals (.* operator)
%
% This function creates the polar intervals representing the 
% element-wise product of two sets of polar intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolarInterval class
%   obj2       : array of objects from the ciat.PolarInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polarInt = ciat.PolarInterval(0,1,2,3) .* ciat.PolarInterval(2,3,4,5);
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
    mustBeA(obj1,["ciat.PolarInterval","double"]);
    mustBeA(obj2,["ciat.PolarInterval","ciat.PolyarxInterval","double"]);

    % Reroute if second is a polyarx interval
    if isa(obj2,"ciat.PolyarxInterval")
        r = obj2 .* obj1;
        return
    end

    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.PolarInterval(abs(obj1), angle(obj1));
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolarInterval(abs(obj2), angle(obj2));
    end 

    r = ciat.PolarInterval(obj1.Abs .* obj2.Abs, obj1.Angle + obj2.Angle);
end
