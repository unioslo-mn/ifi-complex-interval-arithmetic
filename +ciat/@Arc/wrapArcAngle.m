function output = wrapArcAngle(input)

    % Wrap bounds into the -pi pi range and check their order
    wrapInf = wrapToPi(input.Infimum);
    wrapSup = wrapToPi(input.Supremum);

    % Wrap and split angle interval if necessary,unless it's a full circle
    if input.Width < 2*pi
        % Check wrap order
        if wrapInf < wrapSup
            % Wrap the interval
            output = ciat.RealInterval(wrapInf,wrapSup);
        else            % wrapInf is pos and wrapSup is neg
            % Wrap and split the interval
            output = cat(3,ciat.RealInterval(wrapInf,pi), ...
                           ciat.RealInterval(-pi,wrapSup) );
        end
    else
        % Full circle
        output = ciat.RealInterval(-pi,pi);
    end
end