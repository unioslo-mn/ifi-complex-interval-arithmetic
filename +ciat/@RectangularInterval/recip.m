% Reciprocal
function r = recip(obj)
    % Reciprocal of rectangular intervals
    %
    % This function creates the rectangular interval representing the 
    % reciprocal of a set of rectangular intervals
    % _________________________________________________________________________
    % USAGE        
    %   r = recip(obj)
    % _________________________________________________________________________
    % NECESSARY ARGUMENT
    %   obj       : array of objects from the ciat.RectangularInterval class
    % _________________________________________________________________________
    % OPTIONS
    % _________________________________________________________________________
    % EXAMPLES
    %   rectInt = recip(ciat.RectangularInterval(0,1));
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

            % Start with the corners
            x1 = obj(m,n).Real.Infimum;
            x2 = obj(m,n).Real.Supremum;
            y1 = obj(m,n).Imag.Infimum;
            y2 = obj(m,n).Imag.Supremum;
            
            values =  [x1/(x1^2+y1^2), x1/(x1^2+y2^2), x2/(x2^2+y1^2), x2/(x2^2+y2^2)];

            r(m,n).Real = ciat.RealInterval(min(values),max(values));

            if x1 < abs(y1) && x2 > abs(y1)
                r(m,n).Real.Supremum = max(r(m,n).Real.Supremum, abs(y1)/(abs(y1)^2+y1^2));
            end
            if x1 < -abs(y1) && x2 > -abs(y1)
                r(m,n).Real.Infimum = min(r(m,n).Real.Infimum, -abs(y1)/(abs(y1)^2+y1^2));
            end

            if x1 < abs(y2) && x2 > abs(y2)
                r(m,n).Real.Supremum = max(r(m,n).Real.Supremum, abs(y2)/(abs(y2)^2+y2^2));
            end
            if x1 < -abs(y2) && x2 > -abs(y2)
                r(m,n).Real.Infimum = min(r(m,n).Real.Infimum, -abs(y2)/(abs(y2)^2+y2^2));
            end

            if y1 < 0 && y2 > 0
                if x1 > 0
                    r(m,n).Real.Supremum = max(r(m,n).Real.Supremum, 1/x1);
                else
                    r(m,n).Real.Infimum = min(r(m,n).Real.Infimum, 1/x1);
                end
                if x2 > 0
                    r(m,n).Real.Supremum = max(r(m,n).Real.Supremum, 1/x2);
                else
                    r(m,n).Real.Infimum = min(r(m,n).Real.Infimum, 1/x2);
                end
            end
        end

        % Imaginary part
        % It is symmetric in x and y, to a sign change (so this switches infimum and supremum)
        % Basically Im(z) = -Re(iz)

        values = -[y1/(x1^2 + y1^2), y1/(x2^2 + y1^2), y2/(x1^2 + y2^2), y2/(x2^2 + y2^2)];
        
        r(m,n).Imag = ciat.RealInterval(min(values),max(values));

        if y1 < abs(x1) && y2 > abs(x1)
            r(m,n).Imag.Infimum = min(r(m,n).Imag.Infimum, -abs(x1)/(x1^2+x1^2));
        end
        if y1 < -abs(x1) && y2 > -abs(x1)
            r(m,n).Imag.Supremum = max(r(m,n).Imag.Supremum, abs(x1)/(x1^2+x1^2));
        end

        if y1 < abs(x2) && y2 > abs(x2)
            r(m,n).Imag.Infimum = min(r(m,n).Imag.Infimum, -abs(x2)/(x2^2+x2^2));
        end
        if y1 < -abs(x2) && y2 > -abs(x2)
            r(m,n).Imag.Supremum = max(r(m,n).Imag.Supremum, abs(x2)/(x2^2+x2^2));
        end

        if x1 < 0 && x2 > 0
            if y1 > 0
                r(m,n).Imag.Infimum = min(r(m,n).Imag.Infimum, -1/y1);
            else
                r(m,n).Imag.Supremum = max(r(m,n).Imag.Supremum, -1/y1);
            end
            if y2 > 0
                r(m,n).Imag.Infimum = min(r(m,n).Imag.Infimum, -1/y2);
            else
                r(m,n).Imag.Supremum = max(r(m,n).Imag.Supremum, -1/y2);
            end
        end

        % If the interval contains 0, the values are infinite
        if x1 < 0 && x2 > 0 && y1 < 0 && y2 > 0
            r(m,n).Real = ciat.RealInterval(-Inf,Inf);
            r(m,n).Imag = ciat.RealInterval(-Inf,Inf);
            warning('The interval contains 0, the values are infinite')
        end
    end