function angleOut = wrapToPi(angleIn)
    % Wraps angle to the (-pi,pi) interval
    angleOut = mod(angleIn + pi, 2*pi) - pi;
end