
% Set polyarcular interval
a(1) = ciat.Arc(1-1i,1,-0.7*pi,-0.4*pi);
a(2) = ciat.Arc(2+1i,-0.7,0.3*pi,0.5*pi);
a(3) = ciat.Arc(0,0,0,0);
A = ciat.PolyarcularInterval(a);

% Sample interval
A_smp = A.sample(1e3);
a2smp = a(2).sample(100);

% Cast interval
G = ciat.PolygonalInterval(A,'tolerance',0.1);
R = ciat.RectangularInterval(A);
C = ciat.CircularInterval(A);
P = ciat.PolarInterval(A);


%% Plot polygonal interval

figure(1);clf;
subplot(1,2,1);cla;hold on;axis equal
set(gca, 'Position', [0,0.06,0.5,0.9]);

% True interval
A.plot('k','linewidth',1);
fill(real(A_smp),imag(A_smp),0.9*ones(1,3),'EdgeColor','none')
fill(real(a2smp),imag(a2smp),'w','EdgeColor','none')

% Polygon
G.plot('k','linewidth',2);
scatter(real(G.Points),imag(G.Points),50,'ko')


% Rectangle
R.plot('k-.','linewidth',2);
text(0,0.55,'[0,2]+[-2,0.7]i','Interpreter','latex')

% Vertex labels
n=1;xO=-0.05;yO=-0.15;
text(real(G.Points(1))+xO,imag(G.Points(1))+yO, '$P_8\!=\!P_1$','Interpreter','latex','HorizontalAlignment','center')
n=2;xO=0.15;yO=0.1;
text(real(G.Points(n))+xO,imag(G.Points(n))+yO, sprintf('$P_%i$',n),'Interpreter','latex','HorizontalAlignment','center')
n=3;xO=0.1;yO=0;
text(real(G.Points(n))+xO,imag(G.Points(n))+yO, sprintf('$P_%i$',n),'Interpreter','latex','HorizontalAlignment','center')
n=4;xO=-0.1;yO=0.1;
text(real(G.Points(n))+xO,imag(G.Points(n))+yO, sprintf('$P_%i$',n),'Interpreter','latex','HorizontalAlignment','center')
n=5;xO=-0.1;yO=0.1;
text(real(G.Points(n))+xO,imag(G.Points(n))+yO, sprintf('$P_%i$',n),'Interpreter','latex','HorizontalAlignment','center')
n=6;xO=-0.1;yO=-0.1;
text(real(G.Points(n))+xO,imag(G.Points(n))+yO, sprintf('$P_%i$',n),'Interpreter','latex','HorizontalAlignment','center')
n=7;xO=-0.1;yO=0;
text(real(G.Points(n))+xO,imag(G.Points(n))+yO, sprintf('$P_%i$',n),'Interpreter','latex','HorizontalAlignment','center')

% Edge labels
edgeMid = ciat.PolyarcularInterval(G).Edges.Midpoint;
n=1;xO=0.15;yO=-0.06;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',n-0.5),'Interpreter','latex','HorizontalAlignment','center')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex','HorizontalAlignment','center')
n=2;xO=0.15;yO=-0.2;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',n-0.5),'Interpreter','latex','HorizontalAlignment','center')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex','HorizontalAlignment','center')
n=3;xO=0.1;yO=0.25;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',n-0.5),'Interpreter','latex','HorizontalAlignment','center')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex','HorizontalAlignment','center')
n=4;xO=-0.2;yO=0.12;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',n-0.5),'Interpreter','latex','HorizontalAlignment','center')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex','HorizontalAlignment','center')
n=5;xO=-0.2;yO=-0.1;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',n-0.5),'Interpreter','latex','HorizontalAlignment','center')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex','HorizontalAlignment','center')
n=6;xO=-0.2;yO=-0.04;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',n-0.5),'Interpreter','latex','HorizontalAlignment','center')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex','HorizontalAlignment','center')

% Settings
grid on
xlim([-0.2,2.5]);
ylim([-2.3,1]);
xlabel('Real','Position',[-0.3,-2.35])
ylabel('Imag')
fontsize(25,'points')

% Vertex list
annotText = {[sprintf('$P_1=%0.1f,%0.1fi$',real(G.Points(1)),imag(G.Points(1)))],...
             [sprintf('$P_2=%0.1f,%0.1fi$',real(G.Points(2)),imag(G.Points(2)))],...
             [sprintf('$P_3=%0.1f,%0.1fi$',real(G.Points(3)),imag(G.Points(3)))],...
             [sprintf('$P_4=%0.1f,%0.1fi$',real(G.Points(4)),imag(G.Points(4)))],...
             [sprintf('$P_5=%0.1f,%0.1fi$',real(G.Points(5)),imag(G.Points(5)))],...
             [sprintf('$P_6=%0.1f,%0.1fi$',real(G.Points(6)),imag(G.Points(6)))],...
             [sprintf('$P_7=%0.1f,%0.1fi$',real(G.Points(7)),imag(G.Points(7)))]};
annotation('textbox',[.12,.22,.2,.3],'String',annotText, ...
           'BackgroundColor','none',...
           'FontSize',25,...
           'EdgeColor','none',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');


% Tigthess
annotText = {['$\tau(A^\mathcal{R})=' num2str(100*A.Area ./ R.Area,3),'\%$'],...
             ['$\tau(A^\mathcal{G})=' num2str(100*A.Area ./ G.Area,3),'\%$'],...
             };
annotation('textbox',[.17,.55,.1,.1],'String',annotText, ...
           'BackgroundColor','none',...
           'FontSize',25,...
           'EdgeColor','k',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');
%% Plot polyarcular interval

subplot(1,2,2);cla;hold on;axis equal
set(gca, 'Position', [0.5,0.06,0.5,0.9]);

% True interval
fill(real(A_smp),imag(A_smp),0.9*ones(1,3),'EdgeColor','none')
fill(real(a2smp),imag(a2smp),'w','EdgeColor','none')

% Polyarc
A.plot('k','linewidth',2);

% Arcs
a.plot('k-','linewidth',5);
scatter(real(a.Center),imag(a.Center),50,'ko')
scatter(real(a(3).Center),imag(a(3).Center),100,'ko','filled')
ciat.Edge(a.Center,a.Startpoint).plot('k--');
ciat.Edge(a.Center,a.Endpoint).plot('k--');
n=1;xO=-0.1;yO=-0.1;
text(real(a(n).Center)+xO,imag(a(n).Center)+yO,'$O_1=1-1i$','Interpreter','latex','HorizontalAlignment','Center')
text(real(a(n).Center)+xO,imag(a(n).Center)+yO-0.12,'$r_1=1$','Interpreter','latex','HorizontalAlignment','Center')
text(real(a(n).Center)+xO,imag(a(n).Center)+yO-0.24,'$\varphi_1=[-0.7,-0.4]\pi$','Interpreter','latex','HorizontalAlignment','Center')
n=2;xO=0.25;yO=-0.1;
text(real(a(n).Center)+xO,imag(a(n).Center)+yO,'$O_2=2+1i$','Interpreter','latex','HorizontalAlignment','Right')
text(real(a(n).Center)+xO,imag(a(n).Center)+yO-0.12,'$r_2=-0.7$','Interpreter','latex','HorizontalAlignment','Right')
text(real(a(n).Center)+xO,imag(a(n).Center)+yO-0.24,'$\varphi_2=[0.3,0.5]\pi$','Interpreter','latex','HorizontalAlignment','Right')
n=3;xO=0.12;yO=-0.2;
text(real(a(n).Center)+xO,imag(a(n).Center)+yO,'$O_3=0$','Interpreter','latex','HorizontalAlignment','Left')
text(real(a(n).Center)+xO,imag(a(n).Center)+yO-0.12,'$r_3=0$','Interpreter','latex','HorizontalAlignment','Left')
text(real(a(n).Center)+xO,imag(a(n).Center)+yO-0.24,'$\varphi_3=0$','Interpreter','latex','HorizontalAlignment','Left')

% Arc labels
arcMid = A.DefArcs.Midpoint;
n=1;xO=-0.3;yO=-0.1;
text(real(arcMid(n))+xO,imag(arcMid(n))+yO, sprintf('$t=%0.1f $',2*(n-1)+0.5),'Interpreter','latex')
text(real(arcMid(n))+xO,imag(arcMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex')
n=2;xO=-0.4;yO=-0.15;
text(real(arcMid(n))+xO,imag(arcMid(n))+yO, sprintf('$t=%0.1f $',2*(n-1)+0.5),'Interpreter','latex')
text(real(arcMid(n))+xO,imag(arcMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex')
n=3;xO=-0.1;yO=0.2;
text(real(arcMid(n))+xO,imag(arcMid(n))+yO, sprintf('$t=%0.1f $',2*(n-1)+0.5),'Interpreter','latex')
text(real(arcMid(n))+xO,imag(arcMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex')

% Edge labels
edgeMid = A.Edges.Midpoint;
n=1;xO=0.05;yO=0;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',2*n-0.5),'Interpreter','latex')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex')
n=2;xO=-0.15;yO=0.2;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',2*n-0.5),'Interpreter','latex')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex')
n=3;xO=0.1;yO=0.1;
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO, sprintf('$t=%0.1f $',2*n-0.5),'Interpreter','latex')
text(real(edgeMid(n))+xO,imag(edgeMid(n))+yO-0.1, sprintf('$n=%i $',n),'Interpreter','latex')

% Vertex labels
vertMid = A.Vertices.Center;
n=1;xO=0.1;yO=-0.05;
text(real(vertMid(n))+xO,imag(vertMid(n))+yO, '$P_4=P_5$','Interpreter','latex')
n=2;xO=-0.15;yO=0.1;
text(real(vertMid(n))+xO,imag(vertMid(n))+yO, '$P_2$','Interpreter','latex')
n=3;xO=-0.1;yO=-0.1;
text(real(vertMid(n))+xO,imag(vertMid(n))+yO, '$P_3$','Interpreter','latex')
n=4;xO=0.05;yO=0;
text(real(vertMid(n))+xO,imag(vertMid(n))+yO, '$P_6=P_1$','Interpreter','latex')
n=5;xO=-0.2;yO=-0.1;
text(real(vertMid(n))+xO,imag(vertMid(n))+yO, '$P_2$','Interpreter','latex')


% Circular
C.plot('k-.','linewidth',2);
text(0,0.7,'$1-0.7i+1.5e^{i[-\pi,\pi]}$','Interpreter','latex')

% Polar
P.plot('k-.','linewidth',2);
text(1.1,-2.2,'$[0,2.4]e^{i[-0.4,0.5]\pi}$','Interpreter','latex')

% Settings
grid on
xlim([-0.2,2.5]);
ylim([-2.3,1]);
xlabel('Real','Position',[-0.3,-2.35])
ylabel('Imag')
fontsize(25,'points')

% Tigthess
annotText = {['$\tau(A^\mathcal{P})=' num2str(100*A.Area ./ P.Area,3),'\%$'],...
             ['$\tau(A^\mathcal{C})=' num2str(100*A.Area ./ C.Area,3),'\%$'],...
             ['$\tau(A^\mathcal{A})=' num2str(100*A.Area ./ A.Area,3),'\%$'],...
             };
annotation('textbox',[.7,.55,.1,.1],'String',annotText, ...
           'BackgroundColor','none',...
           'FontSize',25,...
           'EdgeColor','k',...
           'HorizontalAlignment','center', ...
           'Interpreter','latex');