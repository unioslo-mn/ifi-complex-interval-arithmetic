clear 
% close all

%%
% Define array
array = biat.SensorArray(   'ElCount',5,...
                            'ElDiameterRatio',0,...
                            'Curvature',0.2,...
                            'TaperType','chebwin',...
                            'TaperParam',20,...
                            'GainError',5/100,...
                            'PhaseError',deg2rad(4),...
                            'SteeringAngle',deg2rad(5));

bp_rec = biat.BeamPattern(array,'rectangular','BeamResolutionDeg',0.35/2);
bp_gon = biat.BeamPattern(array,'polygonal','BeamResolutionDeg',0.35/2, ...
                                            'PolygonTolerance',1e-6);
bp_arc = biat.BeamPattern(array,'polyarcular','BeamResolutionDeg',0.35/2);

%% Calculate sum and measure time

recElementInt = bp_rec.ElementIntervals;
tic
recArrayInt = sum(recElementInt);
recTime = toc;

gonElementInt = bp_gon.ElementIntervals;
tic
gonArrayInt = sum(gonElementInt);
gonTime = toc;

arcElementInt = bp_arc.ElementIntervals;
tic
arcArrayInt = sum(arcElementInt);
arcTime = toc;


%% Plot

% figure;clf
cla;hold on; axis equal

recElementInt.plot('g-')
recArrayInt.plot('g-','linewidth',2)

gonElementInt.plot('r-')
gonArrayInt.plot('r-','linewidth',2)

arcElementInt.plot('b-')
arcArrayInt.plot('b-','linewidth',2)

%%
recArea = bp_rec.ArrayInterval.Area;
gonArea = bp_gon.ArrayInterval.Area;
arcArea = bp_arc.ArrayInterval.Area;

sprintf(['Rectangular area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polygonal area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polyarcular area: %0.4f (tightness: %0.1f%%), Time: %0.1fms'], ...
         recArea, arcArea / recArea * 100, recTime*1e3,...
         gonArea, arcArea / gonArea * 100, gonTime*1e3,...
         arcArea, arcArea / arcArea * 100, arcTime*1e3)


