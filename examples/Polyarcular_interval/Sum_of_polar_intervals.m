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


%% Plot

% figure;clf
cla;hold on; axis equal

bp_rec.ElementIntervals.plot('g-')
bp_rec.ArrayInterval.plot('g-','linewidth',2)

bp_gon.ElementIntervals.plot('r-')
bp_gon.ArrayInterval.plot('r-','linewidth',2)

bp_arc.ElementIntervals.plot('b-')
bp_arc.ArrayInterval.plot('b-','linewidth',2)

%%
recArea = bp_rec.ArrayInterval.Area;
gonArea = bp_gon.ArrayInterval.Area;
arcArea = bp_arc.ArrayInterval.Area;

sprintf(['Rectangular area: %0.4f (tightness: %0.1f%%)\n'...
         'Polygonal area: %0.4f (tightness: %0.1f%%)\n'...
         'Polyarcular area: %0.4f (tightness: %0.1f%%)'], ...
         recArea, arcArea / recArea * 100,...
         gonArea, arcArea / gonArea * 100, ...
         arcArea, arcArea / arcArea * 100)


