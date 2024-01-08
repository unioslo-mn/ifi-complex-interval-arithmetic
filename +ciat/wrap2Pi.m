function angleOut = wrap2Pi(angleIn)
    angleOut = rem(2*pi+angleIn, 2*pi);
end