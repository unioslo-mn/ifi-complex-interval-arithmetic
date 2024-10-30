clear
close all

%%  Generate a random arc

M = 3;
N = 3;
L = 2;

% Create arc array
centers = 3*complex(rand(M,N,L),rand(M,N,L));
radii = 2*rand(M,N,L)-1;
angInf = 2*pi*rand(M,N,L);
angSup = 2*pi*rand(M,N,L);
arcs1 = ciat.Arc(centers(:,:,1), radii(:,:,1), ...
                 ciat.RealInterval(angInf(:,:,1),angSup(:,:,1)));
arcs2 = ciat.Arc(centers(:,:,2), radii(:,:,2), ...
                 ciat.RealInterval(angInf(:,:,2),angSup(:,:,2)));
  

% Add arcs
arcPlus = arcs1 + arcs2;

%% 

% Normalize arcs
figure;clf;hold on;axis equal
arcs1.plot('k:')
arcs2.plot('k--')
arcPlus.plot('b')