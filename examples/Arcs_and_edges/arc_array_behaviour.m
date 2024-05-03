clear
close all

%%  Generate random arcs

M = 4;
N = 3;

% Create arc array
centers = 3*complex(rand(M,N),rand(M,N));
radii = 2*rand(M,N)-1;
angInf = 2*pi*rand(M,N);
angSup = 2*pi*rand(M,N);
arcs = ciat.Arc(centers, radii, ciat.RealInterval(angInf,angSup));

% Plot
figure;
clf;hold on;axis equal
arcs.plot;

%% Plot properties

arcCenter = arcs.Center;
arcRadius = arcs.Radius;
arcMidAngle = arcs.ArcAngle.Midpoint;
arcStart = arcs.Startpoint;
arcStop = arcs.Endpoint;
startPoint = centers + radii .* exp(1j*angInf);
endPoint = centers + radii .* exp(1j*angSup);

radLine = [arcCenter(:) , arcCenter(:) + arcRadius(:) .* exp(1i*arcMidAngle(:))].';

plot(real(arcCenter),imag(arcCenter),'x')
plot(real(radLine),imag(radLine),'--')
plot(real(arcStart),imag(arcStart),'>')
plot(real(arcStop),imag(arcStop),'s')
plot(real(startPoint(:,:,1)),imag(startPoint(:,:,1)),'>')
plot(real(endPoint(:,:,1)),imag(endPoint(:,:,1)),'s')


%% Enclose them in rectangular intervals

arcReal = real(arcs);
arcImag = imag(arcs);
arcAbs = abs(arcs);
arcAng = angle(arcs);

arc2rec = ciat.RectangularInterval(arcReal,arcImag);
arc2pol = ciat.PolarInterval(arcAbs,arcAng);


arc2rec.plot
arc2pol.plot
