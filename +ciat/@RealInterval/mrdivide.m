function r = mrdivide(obj1,obj2)

% Matrix right-division of real intervals (/ operator)
%
% This function creates the real intervals representing the matrix 
% right division of two sets of real intervals 
% (see MATLAB mrdivide function),
% _________________________________________________________________________
% USAGE        
%   r = obj1 / obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.RealInterval class
%   obj2       : array of objects from the ciat.RealInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   realInt = ciat.RealInterval(0,1) / ciat.RealInterval(2,2);
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
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'ciat.RealInterval') == 0
                obj1 = ciat.RealInterval(obj1, obj1);
    end
    if isa(obj2, 'ciat.RealInterval') == 0
        obj2 = ciat.RealInterval(obj2, obj2);
    end 
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    if (M1+N1>2) && (M2+N2>2)
        error('Interval matrix division is not implemented')
    end
    
    % Calculate
    r = obj1 * recip(obj2);
end

