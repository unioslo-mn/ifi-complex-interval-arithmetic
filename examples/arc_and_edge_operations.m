clear
close all

%%  Arc array

% Generate two arcs
arc1 = ciat.Arc(1+1i,1,ciat.RealInterval(1,3));
arc2 = ciat.Arc(2+2i,2,ciat.RealInterval(2,4)); 
arcs = ciat.Arc([1+1i;2+2i],[1;2],ciat.RealInterval([1,3;2,4]));

% Plot
figure;clf;hold on
arc1.plot;
arc2.plot;
arcs.plot
axis equal