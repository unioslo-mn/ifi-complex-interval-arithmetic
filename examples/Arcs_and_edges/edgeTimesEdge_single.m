clear
% close all

%%  Generate two random arc with overlapping log-Gauss map

match = false;
while ~match
    % Create a random edges
    startpoint = complex(rand(2,1),rand(2,1));
    endpoint = complex(rand(2,1),rand(2,1));
    edges = ciat.Edge( startpoint, endpoint );
    % edges = edges - edges.ZeroCrossing;

    % Check log-Gauss map intersection
    match = ~isnan(cap(edges.LogGaussMap));
end  

% Multiply arc
edgeTimesEdge = edges(1) .* edges(2);

%% Plot

% Normalize arc
edgesNorm = edges .* edges.NormFactor;
edgeTimesEdgeNorm = edgeTimesEdge * prod(edges.NormFactor);

% Sample arc
sCnt = 50;
edgesSmp = cell2mat(edges.sample(sCnt).');
edgeSmpTimesEdgeSmp = edgesSmp(:,1) * edgesSmp(:,2).';
edgesNormSmp = cell2mat(edgesNorm.sample(sCnt).');
edgeNormSmpTimesEdgeNormSmp = edgesNormSmp(:,1) * edgesNormSmp(:,2).';
edgeTimesEdgeSmp = edgeTimesEdge.sample(sCnt);
edgeTimesEdgeNormSmp = edgeTimesEdgeNorm.sample(sCnt);

% figure;clf

% Plot arc and edges
subplot(2,2,1);cla;hold on;axis equal;title('Original')
edges.plot('k','lineWidth',2);
edgeTimesEdge.plot('b','lineWidth',3);
scatter(real(edgeSmpTimesEdgeSmp),imag(edgeSmpTimesEdgeSmp),1,'g.')

% Plot zero LGM values when one operand is zero-centered
if edges(1).ZeroCrossing<100*eps
    pZeroLGM = edges(2).findLGM(0);
    plot(real(pZeroLGM),imag(pZeroLGM),'ko')
end
if edges(2).ZeroCrossing<100*eps
    pZeroLGM = edges(1).findLGM(edges.LogGaussMap.mid);
    plot(real(pZeroLGM),imag(pZeroLGM),'ko')
end


% Plot normalized arc
subplot(2,2,2);cla;hold on;axis equal;title('Normalized')
edgesNorm.plot('k','lineWidth',2);
edgesNorm.plotLogGaussMap(.1,'k');
edgeTimesEdgeNorm.plot('b','lineWidth',3);
edgeTimesEdgeNorm.plotLogGaussMap(.1,'b');
scatter(real(edgeNormSmpTimesEdgeNormSmp),imag(edgeNormSmpTimesEdgeNormSmp),1,'g.');
if abs(edges(1).ZeroCrossing) > 100*eps && abs(edges(2).ZeroCrossing) > 100*eps
    fimplicit(@(x,y) y.^2 + 4*x - 4,'k:');
end

% Plot arc in the log-domain
subplot(2,2,3);cla;hold on;axis equal;title('Log')
plot(imag(log(edgesSmp)),real(log(edgesSmp)),'k.');
plot(imag(log(edgeTimesEdgeSmp)),real(log(edgeTimesEdgeSmp)),'b.');
xlabel('Angle')
ylabel('Log-Range')
% Plot zero LGM values when one operand is zero-centered
if edges(1).ZeroCrossing<100*eps
    pZeroLGM = edges(2).findLGM(0);
    plot(imag(log(pZeroLGM)),real(log(pZeroLGM)),'ko');
end
if edges(2).ZeroCrossing<100*eps
    pZeroLGM = edges(1).findLGM(edges.LogGaussMap.mid);
    plot(imag(log(pZeroLGM)),real(log(pZeroLGM)),'ko');
end