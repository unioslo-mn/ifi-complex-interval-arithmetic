clear
close all

%%  Generate random edges

M = 3;
N = 4;

% Create edge array
startpoints = complex(rand(M,N),rand(M,N));
endpoints = complex(rand(M,N),rand(M,N));
edges = ciat.Edge(startpoints,endpoints);

% Plot
figure;
clf;hold on;axis equal
edges.plot;

%% Plot properties

edgeMid = edges.Midpoint;

plot(real(edgeMid),imag(edgeMid),'x')

%% Enclose them in rectangular intervals

edgeReal = real(edges);
edgeImag = imag(edges);
edgeAbs = abs(edges);
edgeAng = angle(edges);

edge2rec = ciat.RectangularInterval(edgeReal,edgeImag);
edge2pol = ciat.PolarInterval(edgeAbs,edgeAng);


edge2rec.plot
edge2pol.plot
