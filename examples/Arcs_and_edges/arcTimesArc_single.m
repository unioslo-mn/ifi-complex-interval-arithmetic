clear
% close all

%%  Generate two random arcs with overlapping log-Gauss map

% Create arc array
match = false;
while ~match
    centers = [1;1] .* complex(rand(2,1),rand(2,1));
    radii = 2*rand(2,1)-1;
    angMid = 2*pi*rand(2,1) - pi;
    angWid = 2*pi*rand(2,1)/5;
    angInf = angMid - angWid/2;
    angSup = angMid + angWid/2;
    arcs = ciat.Arc(centers, radii, ciat.RealInterval(angInf,angSup));
    match = ~isnan(cap(arcs.LogGaussMap));
end  

% Multiply arcs
arcTimesArc = arcs(1) .* arcs(2);

%% Plot

% Normalize arcs
arcsNorm = arcs .* arcs.NormFactor;
arcTimesArcNorm = arcTimesArc * prod(arcs.NormFactor);

% Sample arcs
arcSmp = cell2mat(arcs.sample(50).');
arcSmpTimesArcSmp = arcSmp(:,1) * arcSmp(:,2).';
arcNormSmp = cell2mat(arcsNorm.sample(50).');
arcNormSmpTimesArcNormSmp = arcNormSmp(:,1) .* arcNormSmp(:,2).';
arcTimesArcSmp = arcTimesArc.sample(50);
arcTimesArcNormSmp = arcTimesArcNorm.sample(50);

% figure;clf

% Plot arcs
subplot(2,2,1);cla;hold on;axis equal;title('Original')
arcs.plot('k','lineWidth',2);
arcTimesArc.plot('b','lineWidth',3);
scatter(real(arcSmpTimesArcSmp),imag(arcSmpTimesArcSmp),1,'g.')

% Plot zero LGM values when one operand is zero-centered
if centers(1)==0
    pZeroLGM = arcs(2).findLGM(0);
    plot(real(pZeroLGM),imag(pZeroLGM),'ko')
end
if centers(2)==0
    pZeroLGM = arcs(1).findLGM(0);
    plot(real(pZeroLGM),imag(pZeroLGM),'ko')
end


% Plot normalized arcs
subplot(2,2,2);cla;hold on;axis equal;title('Normalized')
arcsNorm.plot('k','lineWidth',2);
arcsNorm.plotLogGaussMap(.1,'k');
arcTimesArcNorm.plot('b','lineWidth',3);
% arcTimesArcNorm.plotMap(0,.1,'b');
scatter(real(arcNormSmpTimesArcNormSmp),imag(arcNormSmpTimesArcNormSmp),1,'g.');
fimplicit(@(x,y) (x-arcsNorm(1).Center).^2 + y.^2 - arcsNorm(1).Radius.^2,'k:');
fimplicit(@(x,y) (x-arcsNorm(2).Center).^2 + y.^2 - arcsNorm(2).Radius.^2,'k:');
fimplicit(@(x,y) (x-arcTimesArcNorm.Center).^2 + y.^2 - arcTimesArcNorm.Radius.^2,'b:');
if all(centers~=0)
    R1 = arcsNorm(1).Radius;
    R2 = arcsNorm(2).Radius;
    fimplicit(@(x,y) (x.^2+y.^2-2*x+(1-R1^2)*(1-R2^2)).^2 - ...
                      4*R1^2*R2^2*(x.^2+y.^2),'k:');
end

% Plot arcs in the log-domain
subplot(2,2,3);cla;hold on;axis equal;title('Log')
plot(imag(log(arcSmp)),real(log(arcSmp)),'k.');
plot(imag(log(arcTimesArcSmp)),real(log(arcTimesArcSmp)),'b.');
xlabel('Angle')
ylabel('Log-Range')
% Plot zero LGM values when one operand is zero-centered
if centers(1)==0
    pZeroLGM = arcs(2).findLGM(0);
    plot(imag(log(pZeroLGM)),real(log(pZeroLGM)),'ko');
end
if centers(2)==0
    pZeroLGM = arcs(1).findLGM(0);
    plot(imag(log(pZeroLGM)),real(log(pZeroLGM)),'ko');
end
