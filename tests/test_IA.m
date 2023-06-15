% Unit test for the Matlab toolbox for complex interval arithmetic
clear all;
close all;

res = runtests(["RealIntervalTest", "RectangularIntervalTest", "PolarIntervalTest"]);
table(res)

%table(runtests("RealIntervalTest"))
%table(runtests("RectangularIntervalTest"))
%table(runtests("PolarIntervalTest"))