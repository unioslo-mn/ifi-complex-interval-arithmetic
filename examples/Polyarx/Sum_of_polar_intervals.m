clear 
% close all

%% Generate random polar intervals

% Define random polar intervals
N = 3;
absMin = rand(N,1);
absMax = absMin + rand(N,1)/10;
angMin = 2*pi*rand(N,1);
angMax = angMin + 2*pi*rand(N,1)/20;
pI = ciat.PolarInterval(absMin,absMax,angMin,angMax);

% Convert them to various other formats
pcI = ciat.RectangularInterval(pI);
pgI = ciat.PolygonalInterval(pI);
paI = ciat.PolyarcularInterval(pI);
pxI = ciat.PolyarxInterval(pI);
pIsmp = pI.sample(10);

%% Sum polar intervals

% Sum samples
pIsmpSum = 0;
for n = 1:N
    pIsmpSum = pIsmpSum(:) + pIsmp{n}.';
end
pIsmpSum = pIsmpSum(:);

% Sum polygonal intervals and measure time
tic
pcIsum = sum(pcI);
recTime = toc;

% Sum polygonal intervals and measure time
tic
pgIsum = sum(pgI); 
gonTime = toc;

% Sum polyarcular intervals and measure time
tic
paIsum = sum(paI);
arcTime = toc;

% Sum convex polygonal intervals and measure time
tic
pxIsum = sum(pxI);
arxTime = toc;

% Measure areas
recArea = pcIsum.Area;
gonArea = pgIsum.Area;
arcArea = paIsum.Area;
arxArea = pxIsum.Area;

% Check if points are inside the sum
smpInside = paIsum.isin(pIsmpSum);

% Print report
sprintf(['Rectangular interval: Area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polygonal interval area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polyarx area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n',...
         'Polyarcular area: %0.4f (tightness: %0.1f%%), Time: %0.1fms'], ...
         recArea, arcArea / recArea * 100, recTime*1e3,...
         gonArea, arcArea / gonArea * 100, gonTime*1e3,...
         arxArea, arcArea / arxArea * 100, arxTime*1e3,...
         arcArea, arcArea / arcArea * 100, arcTime*1e3)

%% Plot

% figure;clf
cla;hold on;axis equal
paI.plot('b');
paI.plotGaussMap(0.05,'c');
paIsum.plot('k','linewidth',2);
paIsum.plotGaussMap(0.05,'k');
pcIsum.plot('c','linewidth',2);
pgIsum.plot('r','linewidth',2);
pxIsum.plot('y','linewidth',2);
scatter(real([pIsmp{:}]),imag([pIsmp{:}]),10,'co');
scatter(real(pIsmpSum),imag(pIsmpSum),1,'k.');

