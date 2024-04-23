function r = tan(obj)

% Tangent of rectangular intervals
%
% This function creates the rectangular intervals representing the tangent
% of a rectangular intervals (see MATLAB tan function),
% _________________________________________________________________________
% USAGE        
%   r = tan(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.RectangularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   rectInt = tan(ciat.RectangularInterval(0,1,2,4))
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
    mustBeA(obj,["ciat.RectangularInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M,N] = size(obj);
    
    % Turn scalars to degenerate intervals
    if isa(obj, 'double')
        obj = ciat.RectangularInterval(obj1);
    end
            
    % Loop throught the arrays
    r(M,N) = ciat.RectangularInterval;
    for m = 1:M
        for n = 1:N
            % Real part
            % Candidate points are the corners of the rectangle, points on the edge with zero imaginary part
            x1 = obj(m,n).Real.Infimum;
            x2 = obj(m,n).Real.Supremum;
            y1 = obj(m,n).Imag.Infimum;
            y2 = obj(m,n).Imag.Supremum;

            
            if obj(m,n).Real.Width > pi
                % If real interval is larger than pi, then the extreme values are both attained in the horizontal edges (because real(tan(x+iy)) is pi-periodic along x)
                values = [csch(2*y1), csch(2*y2), -csch(2*y1), -csch(2*y2)];
                r(m,n).Real = ciat.RealInterval(min(values),max(values));
            else
                % First assign to the value of the function at the corners
                values = [sin(2*x1)/(cos(2*x1) + cosh(2*y1)), sin(2*x1)/(cos(2*x1) + cosh(2*y2)), sin(2*x2)/(cos(2*x2) + cosh(2*y1)), sin(2*x2)/(cos(2*x2) + cosh(2*y2))];
                r(m,n).Real = ciat.RealInterval(min(values),max(values));
                
                % Else, we need to check if the real interval contains a point of the form k*pi + atan(1/tanh(y)) for some integer k
                % If it does, then the maximum values are attained in the horizontal edges
                % If it does not, then the maximum values are attained in the corners

                xstar = abs(atan(1/tanh(y1)));
                k = ceil((x1 - xstar)/pi);
                if k*pi + xstar > x1 && k*pi + xstar < x2
                    r(m,n).Real.Supremum = max(abs(csch(2*y1)), r(m,n).Real.Supremum);
                end
                % Same check for -xstar
                k = ceil((x1 + xstar)/pi);
                if k*pi - xstar > x1 && k*pi - xstar < x2
                    r(m,n).Real.Infimum = min(-abs(csch(2*y1)), r(m,n).Real.Infimum);
                end

                % Same for y2
                xstar = abs(atan(1/tanh(y2)));
                k = ceil((x1 - xstar)/pi);
                if k*pi + xstar > x1 && k*pi + xstar < x2
                    r(m,n).Real.Supremum = max(abs(csch(2*y2)), r(m,n).Real.Supremum);
                end
                % Same check for -xstar
                k = ceil((x1 + xstar)/pi);
                if k*pi - xstar > x1 && k*pi - xstar < x2
                    r(m,n).Real.Infimum = min(-abs(csch(2*y2)), r(m,n).Real.Infimum); 
                end
                

            end

            % Imaginary part

            if obj(m,n).Real.Width > pi
                values = [sinh(2*y1)/(cosh(2*y1)+1), sinh(2*y1)/(cosh(2*y1)-1), sinh(2*y2)/(cosh(2*y2)+1), sinh(2*y2)/(cosh(2*y2)-1)];
                r(m,n).Imag = ciat.RealInterval(min(values),max(values));
            else
                % First assign to the value of the function at the corners
                values = [sinh(2*y1)/(cos(2*x1)+cosh(2*y1)), sinh(2*y1)/(cos(2*x2)+cosh(2*y1)), sinh(2*y2)/(cos(2*x1)+cosh(2*y2)), sinh(2*y2)/(cos(2*x2)+cosh(2*y2))];
                r(m,n).Imag = ciat.RealInterval(min(values),max(values));

                % Else, we need to check if the real interval contains a point of the form k*pi or k*pi + pi/2 for some integer k

                % Check for k*pi
                k = ceil(x1/pi);
                if k*pi > x1 && k*pi < x2
                    if y1 > 0
                        % This is a minimum
                        r(m,n).Imag.Infimum = min(sinh(2*y1)/(cosh(2*y1)+1), r(m,n).Imag.Infimum);
                    else
                        % This is a maximum
                        r(m,n).Imag.Supremum = max(sinh(2*y1)/(cosh(2*y1)+1), r(m,n).Imag.Supremum);
                    end
                    if y2 > 0
                        % This is a minimum
                        r(m,n).Imag.Infimum = min(sinh(2*y2)/(cosh(2*y2)+1), r(m,n).Imag.Infimum);
                    else
                        % This is a maximum
                        r(m,n).Imag.Supremum = max(sinh(2*y2)/(cosh(2*y2)+1), r(m,n).Imag.Supremum);
                    end
                end

                % Check for k*pi + pi/2
                k = ceil((x1 - pi/2)/pi);
                if k*pi + pi/2 > x1 && k*pi + pi/2 < x2
                    if y1 > 0
                        % This is a maximum
                        r(m,n).Imag.Supremum = max(sinh(2*y1)/(cosh(2*y1)-1), r(m,n).Imag.Supremum);
                    else
                        % This is a minimum
                        r(m,n).Imag.Infimum = min(sinh(2*y1)/(cosh(2*y1)-1), r(m,n).Imag.Infimum);
                    end
                    if y2 > 0
                        % This is a maximum
                        r(m,n).Imag.Supremum = max(sinh(2*y2)/(cosh(2*y2)-1), r(m,n).Imag.Supremum);
                    else
                        % This is a minimum
                        r(m,n).Imag.Infimum = min(sinh(2*y2)/(cosh(2*y2)-1), r(m,n).Imag.Infimum);
                    end
                end

                % Looking at vertical edges, there are possible extreme points along it
                if cos(2*x1) < 0
                    ystar1 = acosh(-1/cos(2*x1))/2;
                    if ystar1 > y1 && ystar1 < y2
                        r(m,n).Imag.Supremum = max(sinh(2*ystar1)/(cos(2*x1)+cosh(2*ystar1)), r(m,n).Imag.Supremum);
                    end
                    ystar2 = -acosh(-1/cos(2*x1))/2;
                    if ystar2 > y1 && ystar2 < y2
                        r(m,n).Imag.Infimum = min(sinh(2*ystar2)/(cos(2*x1)+cosh(2*ystar2)), r(m,n).Imag.Infimum);
                    end
                end
                if cos(2*x2) < 0
                    ystar1 = acosh(-1/cos(2*x2))/2;
                    if ystar1 > y1 && ystar1 < y2
                        r(m,n).Imag.Supremum = max(sinh(2*ystar1)/(cos(2*x2)+cosh(2*ystar1)), r(m,n).Imag.Supremum);
                    end
                    ystar2 = -acosh(-1/cos(2*x2))/2;
                    if ystar2 > y1 && ystar2 < y2
                        r(m,n).Imag.Infimum = min(sinh(2*ystar2)/(cos(2*x2)+cosh(2*ystar2)), r(m,n).Imag.Infimum);
                    end
                end
            end



            % If 0 is in the imaginary interval
            if y1 <= 0 && y2 >= 0
                % If pi/2 + k*pi is in the real interval for some k, then it reach singular points
                % If not, then the maximum values are attained in the corners (because tan is monotone)
                k = ceil((x1 - pi/2)/pi);
                if k*pi + pi/2 > x1 && k*pi + pi/2 < x2
                    % Singular point reached
                    r(m,n).Real = ciat.RealInterval(-Inf,Inf);
                    warning('Singular point reached in tan')
                else
                    r(m,n).Real.Infimum = min(tan(x1), r(m,n).Real.Infimum);
                    r(m,n).Real.Supremum = max(tan(x2), r(m,n).Real.Supremum);
                end

                % For the imaginary part, we need to check if the real interval contains a point of the form k*pi + pi/2 for some integer k

                k = ceil((x1-pi/2)/pi);
                if k*pi + pi/2 > x1 && k*pi + pi/2 < x2
                    r(m,n).Imag = ciat.RealInterval(-Inf,Inf);
                end
            end
            % Basically is r(m,n) = union([r(m,n), tan(r(m,n).Real]))
            % Could be replaced by that

        end
    end
end

        