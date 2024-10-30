function angleOut = wrapTo2Pi(angleIn)
    % Wraps angle to the (0,2*pi) interval
    angleOut = rem(angleIn, 2*pi) + (angleIn<0)*2*pi;
    % angleOut = rem(2*pi+angleIn, 2*pi);
end