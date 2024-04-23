function angleOut = wrapTo2Pi(angleIn)
    angleOut = rem(angleIn, 2*pi) + (angleIn<0)*2*pi;
    % angleOut = rem(2*pi+angleIn, 2*pi);
end