% function output = wrapArcAngle(input)
function angles = wrapArcAngle(angles)
    % This function wraps the given angle interval, so the
    % infimum ends up in the (-2pi,pi) range and the
    % supremum ends up in the (-pi,2pi) range, 
    % while a full circle is identified with the (-pi,pi) interval

    % Set valid ranges for properties
    infRange = ciat.RealInterval(-2*pi,pi);
    supRange = ciat.RealInterval(-pi,2*pi);
    widRange = ciat.RealInterval(0,2*pi);

    % Extract parameters
    angInf = angles.Infimum;
    angSup = angles.Supremum;
    angWid= angles.Width;

    % Check which the angle intervals needs wrapping
    maskValid = infRange.isin(angInf) & ...
                supRange.isin(angSup) & ...
                widRange.isin(angWid);
    maskWrap = ~maskValid;

    if any(maskWrap,'all')
        % Initialize output arrays
        angInf = angInf(maskWrap);
        angSup = angSup(maskWrap);
        angWid = angWid(maskWrap);

        % Extract properties
        angInfPi = wrapToPi( angInf );
        angSupPi = wrapToPi( angSup );
        angPiOrder = ( angInfPi <= angSupPi );

        % Create case masks
        wrapSimple = angPiOrder;
        wrapShift0 = angWid < pi & ~angPiOrder;
        wrapShift1 = angWid >=pi & angWid < 2*pi & ~angPiOrder;
        wrapFull   = angWid >= 2*pi;

        % Wrap interval according to the case:
        %   - width less than 2pi, wrapped order is correct
        if any(wrapSimple,'all')
            angInf(wrapSimple) = angInfPi(wrapSimple);
            angSup(wrapSimple) = angSupPi(wrapSimple);
        end
        %  - width less than pi, wrapped order is incorrect
        if any(wrapShift0,'all')
            angInf(wrapShift0) = angInfPi(wrapShift0) - 2*pi;
            angSup(wrapShift0) = angSupPi(wrapShift0);
        end
        %  - width between pi and 2pi, wrapped order is incorrect
        if any(wrapShift1,'all')
            angInf(wrapShift1) = wrapTo2Pi(angInf(wrapShift1)) - 2*pi;
            angSup(wrapShift1) = angInf(wrapShift1) + angWid(wrapShift1);
        end
        %  - width is larger or equal to 2pi
        if any(wrapFull,'all')
            angInf(wrapFull) = -pi;
            angSup(wrapFull) = pi;
        end

        % Assign wrapped values to output
        angles(maskWrap).Infimum = angInf;
        angles(maskWrap).Supremum= angSup;
    end

    


    %% Old solution
    % % Check the width of the interval
    % widthPi = floor(input.Width / pi);
    % 
    % % Wrap bounds into the -pi pi range and check their order
    % wrapInf = wrapToPi(input.Infimum);
    % wrapSup = wrapToPi(input.Supremum);
    % wrapOrder = ( wrapInf <= wrapSup );
    % 
    % % Wrap interval into the -2pi +2pi range
    % switch widthPi
    %     case 0  % Less than a Pi width
    %         if wrapOrder    % Normal order
    %             % Wrap both input to (-pi,pi)
    %             angInf = wrapInf;
    %             angSup = wrapSup;
    %         else            % wrapInf is pos and wrapSup is neg
    %             % Wrap the infimum to (-2pi,-pi),
    %             % wrap the supremum to (-pi,0)
    %             angInf = wrapInf - 2*pi;
    %             angSup = wrapSup;
    %         end
    % 
    %     case 1 % More than a Pi width but less than 2Pi
    %         if wrapOrder    % Normal order
    %             % Wrap both input to (-pi,pi)
    %             angInf = wrapInf;
    %             angSup = wrapSup;
    %         else            % wrapped input have the same sign 
    %                         % and wrapSup is lower than wrapInf
    %             % Wrap the infimum to (-2pi,0),
    %             % wrap the supremum to (-pi,2*pi)
    %             angInf = wrapTo2Pi(input.Infimum) - 2*pi;
    %             angSup = angInf + input.Width;
    %         end
    %     otherwise           % full circle
    %         angInf = -pi;
    %         angSup = pi;
    % end
    % 
    % output = ciat.RealInterval(angInf,angSup);
        
end