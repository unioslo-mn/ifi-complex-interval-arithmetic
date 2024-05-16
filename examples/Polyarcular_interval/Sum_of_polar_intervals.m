clear 
% close all

%%

% Define random polar intervals
N = 3;
absMin = rand(N,1);
absMax = absMin + rand(N,1)/10;
angMin = 2*pi*rand(N,1);
angMax = angMin + 2*pi*rand(N,1)/20;
pI = ciat.PolarInterval(absMin,absMax,angMin,angMax);
paI = ciat.PolyarcularInterval(pI);
paIc = paI.convexify;
pIsmp = pI.sample(10);

% Sum samples
pIsmpSum = 0;
for n = 1:N
    pIsmpSum = pIsmpSum(:) + pIsmp{n}.';
end
pIsmpSum = pIsmpSum(:);

% Sum intervals
paIsum = sum(paI);
paIcSum = sum(paIc);

% Check if points are inside the sum
smpInside = paIsum.isin(pIsmpSum);

%% 
% figure;clf
cla;hold on;axis equal
paI.plot('b')
paI.plotGaussMap(0.05,'c');
paIsum.plot('k','linewidth',2)
paIsum.plotGaussMap(0.05,'k');
paIcSum.plot('y','linewidth',2)
paIcSum.plotGaussMap(0.05,'k');
scatter(real([pIsmp{:}]),imag([pIsmp{:}]),10,'co')
scatter(real(pIsmpSum),imag(pIsmpSum),1,'k.')


