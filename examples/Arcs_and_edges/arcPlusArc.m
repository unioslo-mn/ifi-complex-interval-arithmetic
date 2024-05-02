clear
% close all

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

% startPoint = centers + radii .* exp(1j*angInf);
% endPoint = centers + radii .* exp(1j*angSup);
% 
% figure;clf;hold on;axis equal
% arcs1.plot
% plot(real(startPoint(:,:,1)),imag(startPoint(:,:,1)),'>')
% plot(real(endPoint(:,:,1)),imag(endPoint(:,:,1)),'s')


%% 



%% 

% Normalize arcs
arcsNorm = arcs + (-arcs.Center);
arcPlusNorm = arcPlus + (-arcPlus.Center);

% figure;
% Plot arcs
subplot(1,2,1);cla;hold on;axis equal;title('Original')
arcs.plot('k','lineWidth',2);
arcPlus.plot('b','lineWidth',2);
% Plot normalized arcs
subplot(1,2,2);cla;hold on;axis equal;title('Normalized')
arcsNorm.plot('k','lineWidth',2)
arcsNorm.plotMap(0,.1,'k');
arcPlusNorm.plot('b','lineWidth',2)
arcPlusNorm.plotMap(0,.1,'b');
fimplicit(@(x,y) x.^2 + y.^2 - radii(1).^2,'k:')
fimplicit(@(x,y) x.^2 + y.^2 - radii(2).^2,'k:')
fimplicit(@(x,y) x.^2 + y.^2 - sum(radii).^2,'b:')


