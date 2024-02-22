function angleOut = wrapToPi(angleIn)
    angleOut = mod(angleIn + pi, 2*pi) - pi;
end