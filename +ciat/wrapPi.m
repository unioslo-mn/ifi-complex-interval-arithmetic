function angleOut = wrapPi(angleIn)
    angleOut = mod(angleIn + pi, 2*pi) - pi;
end