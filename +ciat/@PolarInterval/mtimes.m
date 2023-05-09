function r = mtimes(obj1,obj2)

% Matrix multiplication of polar intervals (* operator)
%
% This function creates the polar intervals representing the matrix 
% product of two sets of polar intervals (see MATLAB mtimes function),
% which it achieves by using the times (.*) and plus (+) functions.
% Since polar intervals have no implemented plus function this function
% leads to the times function (.* operator) instead.
% _________________________________________________________________________
% USAGE        
%   r = obj1 * obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolarInterval class
%   obj2       : array of objects from the ciat.PolarInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polarInt = ciat.PolarInterval(0,1,2,3) * ciat.PolarInterval(2,3,4,5);
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
    mustBeA(obj1,"ciat.PolarInterval");
    mustBeA(obj2,"ciat.PolarInterval");
   
    % Calculate product using the times function
    warning(['No plus method implemented for polar intervals, ',...
             'using times function instead (.*)'])
    r = obj1 .* obj2;
end 