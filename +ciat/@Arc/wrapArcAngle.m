function output = wrapArcAngle(input)
    % This function wraps the given angle interval, so the
    % infimum ends up in the (-2pi,pi) range and the
    % supremum ends up in the (-pi,2pi) range, 
    % while a full circle is identified with the (-pi,pi) interval

    % Check the width of the interval
    widthPi = floor(input.Width / pi);

    % Wrap bounds into the -pi pi range and check their order
    wrapInf = wrapToPi(input.Infimum);
    wrapSup = wrapToPi(input.Supremum);
    wrapOrder = ( wrapInf <= wrapSup );

    % Wrap interval into the -2pi +2pi range
    switch widthPi
        case 0  % Less than a Pi width
            if wrapOrder    % Normal order
                % Wrap both input to (-pi,pi)
                angInf = wrapInf;
                angSup = wrapSup;
            else            % wrapInf is pos and wrapSup is neg
                % Wrap the infimum to (-2pi,-pi),
                % wrap the supremum to (-pi,0)
                angInf = wrapInf - 2*pi;
                angSup = wrapSup;
            end

        case 1 % More than a Pi width but less than 2Pi
            if wrapOrder    % Normal order
                % Wrap both input to (-pi,pi)
                angInf = wrapInf;
                angSup = wrapSup;
            else            % wrapped input have the same sign 
                            % and wrapSup is lower than wrapInf
                % Wrap the infimum to (-2pi,0),
                % wrap the supremum to (-pi,2*pi)
                angInf = wrapTo2Pi(input.Infimum) - 2*pi;
                angSup = angInf + input.Width;
            end
        otherwise           % full circle
            angInf = -pi;
            angSup = pi;
    end

    output = ciat.RealInterval(angInf,angSup);
        
end