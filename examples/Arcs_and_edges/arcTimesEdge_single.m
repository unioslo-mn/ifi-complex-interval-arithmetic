clear
% close all

%%  Generate two random arc with overlapping log-Gauss map

match = false;
while ~match
    % Create a random arc
    center = complex(rand,rand);
    radius = 2*rand-1;
    angMid = 2*pi*rand - pi;
    angWid = 2*pi*rand/5;
    angInf = angMid - angWid/2;
    angSup = angMid + angWid/2;
    arc = ciat.Arc(center, radius, ciat.RealInterval(angInf,angSup));

    % Create a random edge
    startpoint = complex(rand,rand);
    endpoint = complex(rand,rand);
    edge = ciat.Edge( startpoint, endpoint );
    % edge = edge - edge.ZeroCrossing;

    % Check log-Gauss map intersection
    match = ~isnan(cap(arc.LogGaussMap,edge.LogGaussMap));
end  

% Multiply arc
arcTimesEdge = arc .* edge;

%% Plot

% Normalize arc
arcNorm = arc .* arc.NormFactor;
edgeNorm = edge .* edge.NormFactor;
arcTimesEdgeNorm = arcTimesEdge * arc.NormFactor * edge.NormFactor;

% Sample arc
sCnt = 50;
arcSmp = arc.sample(sCnt);
edgeSmp = edge.sample(sCnt);
arcSmpTimesEdgeSmp = arcSmp * edgeSmp.';
arcNormSmp = arcNorm.sample(sCnt);
edgeNormSmp = edgeNorm.sample(sCnt);
arcNormSmpTimesEdgeNormSmp = arcNormSmp .* edgeNormSmp.';
arcTimesEdgeSmp = arcTimesEdge.sample(sCnt);
arcTimesEdgeNormSmp = arcTimesEdgeNorm.sample(sCnt);

% figure;clf

% Plot arc and edge
subplot(2,2,1);cla;hold on;axis equal;title('Original')
arc.plot('k','lineWidth',2);
edge.plot('k','lineWidth',2);
arcTimesEdge.plot('b','lineWidth',3);
scatter(real(arcSmpTimesEdgeSmp),imag(arcSmpTimesEdgeSmp),1,'g.')

% Plot zero LGM values when one operand is zero-centered
if center==0
    pZeroLGM = edge.findLGM(0);
    plot(real(pZeroLGM),imag(pZeroLGM),'ko')
end
if edge.ZeroCrossing<100*eps
    pZeroLGM = arc.findLGM(edge.LogGaussMap.mid);
    plot(real(pZeroLGM),imag(pZeroLGM),'ko')
end


% Plot normalized arc
subplot(2,2,2);cla;hold on;axis equal;title('Normalized')
arcNorm.plot('k','lineWidth',2);
edgeNorm.plot('k','lineWidth',2);
arcNorm.plotLogGaussMap(.1,'k');
edgeNorm.plotLogGaussMap(.1,'k');
arcTimesEdgeNorm.plot('b','lineWidth',3);
arcTimesEdgeNorm.plotLogGaussMap(.1,'b');
scatter(real(arcNormSmpTimesEdgeNormSmp),imag(arcNormSmpTimesEdgeNormSmp),1,'g.');
fimplicit(@(x,y) (x-arcNorm(1).Center).^2 + y.^2 - arcNorm(1).Radius.^2,'k:');
fimplicit(@(x,y) (x-arcNorm(2).Center).^2 + y.^2 - arcNorm(2).Radius.^2,'k:');
fimplicit(@(x,y) (x-arcTimesEdgeNorm.Center).^2 + y.^2 - arcTimesEdgeNorm.Radius.^2,'b:');
if arc.Center~=0 && abs(edge.ZeroCrossing) > 100*eps
    R = arcNorm.Radius;
    fimplicit(@(x,y) R^4 - R^2*x.^2 - R^2*y.^2 + 2*R^2*x - 2*R^2 + ...
                     x.^2-2*x + 1,'k:');
end

% Plot arc in the log-domain
subplot(2,2,3);cla;hold on;axis equal;title('Log')
plot(imag(log(arcSmp)),real(log(arcSmp)),'k.');
plot(imag(log(edgeSmp)),real(log(edgeSmp)),'k.');
plot(imag(log(arcTimesEdgeSmp)),real(log(arcTimesEdgeSmp)),'b.');
xlabel('Angle')
ylabel('Log-Range')
% Plot zero LGM values when one operand is zero-centered
if center==0
    pZeroLGM = edge.findLGM(0);
    plot(imag(log(pZeroLGM)),real(log(pZeroLGM)),'ko');
end
if edge.ZeroCrossing<100*eps
    pZeroLGM = arc.findLGM(edge.LogGaussMap.mid);
    plot(imag(log(pZeroLGM)),real(log(pZeroLGM)),'ko');
end