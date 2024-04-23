clear
close all

%% Define two arcs

arc1 = ciat.Arc(1+1i,1,ciat.RealInterval(1,2));
arc2 = ciat.Arc(2+2i,2,ciat.RealInterval(3,4)); 

% Plot
figure;clf;hold on;axis equal
arc1.plot;
arc2.plot;

%%  Arc array behavioour

% Create arc array
arcs = ciat.Arc([1+1i;2+2i],[1;2],ciat.RealInterval([1;2],[3;4]));

% Plot
arcs.plot('LineStyle','--','LineWidth',2);

%% Arc addition

% Add arcs
arcSum = arc1 + arc2;

% Add a vertex to the arc
arcTranslate = arc1 + 1;

arcSum.plot
arcTranslate.plot