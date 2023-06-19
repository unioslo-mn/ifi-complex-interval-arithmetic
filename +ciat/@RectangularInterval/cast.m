function outObj = cast(inObj)

% Cast complex intervals of other types to rectangular interval type
%
% This function takes one or more complex intervals of another type
% and creates the smallest inclusive interval(s) of the rectangular
% interval type.
% _________________________________________________________________________
% USAGE        
%   outObj = ciat.RectangularInterval.cast(inObj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   inObj       : object of one of the following types:
%                   - ciat.CircularInterval
%                   - ciat.PolarInterval
%                   - ciat.PolygonalInterval               
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = ciat.RectangularInterval(ciat.PolarInterval(1,3,2,4));
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________
    [M,N] = size(inObj);

    switch class(inObj)
        case 'double'
            outReal = ciat.RealInterval(real(inObj));
            outImag = ciat.RealInterval(imag(inObj));
        case 'ciat.RealInterval'
            outReal = inObj;
            outImag = ciat.RealInterval(zeros(size(inObj)));
        case {'ciat.CircularInterval',...
              'ciat.PolarInterval',...
              'ciat.PolygonalInterval'}
            outReal = real(inObj);
            outImag = imag(inObj);
        otherwise
            error('Invalid input type')
    end
	outObj = ciat.RectangularInterval(outReal,outImag);       
end

