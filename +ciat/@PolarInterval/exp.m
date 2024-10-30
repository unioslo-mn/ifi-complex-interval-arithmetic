function r = exp(obj)
   
% Exponential value of polar intervals.
% _________________________________________________________________________
% USAGE        
%   r = exp(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.PolarInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polarInt = exp(ciat.PolarInterval(2,3,4,5));
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    r = ciat.PolarInterval(exp(obj.Real), obj.Imag);
end
