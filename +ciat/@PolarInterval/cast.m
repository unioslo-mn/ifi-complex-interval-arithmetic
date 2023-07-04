function outObj = cast(inObj)

% Cast complex intervals of other types to polar interval type
%
% This function takes one or more complex intervals of another type
% and creates the smallest inclusive interval(s) of the polar
% interval type.
% _________________________________________________________________________
% USAGE        
%   outObj = ciat.PolarInterval.cast(inObj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   inObj       : object of one of the following types:
%                   - ciat.RectangularInterval
%                   - ciat.CircularInterval
%                   - ciat.PolygonalInterval               
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polarInt = ciat.PolarInterval(ciat.RectangularInterval(1,3,2,4));
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    switch class(inObj)
        case 'double'
            outAbs = ciat.RealInterval(abs(inObj));
            outAngle = ciat.RealInterval(angle(inObj));
        case {'ciat.RectangularInterval',...
              'ciat.CircularInterval',...
              'ciat.PolygonalInterval'}
            outAbs = abs(inObj);
            outAngle = angle(inObj);   
        otherwise
            error('Invalid input type')
    end
	outObj = ciat.PolarInterval(outAbs,outAngle);       
end

