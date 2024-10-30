function r = rdivide(obj1,obj2)

% Element-wise right-division of rectangular intervals (./ operator)
%
% This function creates the rectangular intervals representing the 
% element-wise right division of two sets of rectangular intervals 
% (see MATLAB rdivide function),
% _________________________________________________________________________
% USAGE        
%   r = obj1 ./ obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.RectangularInterval class
%   obj2       : array of objects from the ciat.RectangularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = ciat.RectangularInterval(0,1,2,3) ./ ciat.RectangularInterval(2,4,6,8);
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
    mustBeA(obj1,["ciat.RectangularInterval","double"]);
    mustBeA(obj2,["ciat.RectangularInterval","double"]);
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'ciat.RectangularInterval') == 0
                obj1 = ciat.RectangularInterval(obj1);
    end
    if isa(obj2, 'ciat.RectangularInterval') == 0
        obj2 = ciat.RectangularInterval(obj2);
    end 
    
    % Calculate
    % This is not optimal
    r = obj1 .* obj2.recip;
end

