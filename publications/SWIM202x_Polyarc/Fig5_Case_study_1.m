clear
% close all

%% Calculate polar circular product

% Set intervals
A = ciat.CircularInterval(1+0.5i,0.10);
E = ciat.PolarInterval(1.6 , 2.0 , 0.05*pi, 0.08*pi);

% Wrap intervals
Ar = ciat.RectangularInterval(A);
Er = ciat.RectangularInterval(E);
Ec = ciat.CircularInterval(E);

% Calculate product intervals
EAg = ciat.PolygonalInterval(E,A,'tolerance',0.1);
EAa = ciat.PolyarcularInterval(E,A);
EAx = ciat.PolyarxInterval(E,A);
EAr = Ar * Er;
EAc = A * Ec;

% Calculate corner circles
Ap(4,1) = ciat.CircularInterval;
Ap(1) = A * E.Abs.inf*exp(1j*E.Angle.inf);
Ap(2) = A * E.Abs.sup*exp(1j*E.Angle.inf);
Ap(3) = A * E.Abs.inf*exp(1j*E.Angle.sup);
Ap(4) = A * E.Abs.sup*exp(1j*E.Angle.sup);

% Calculate rotated polar
Ep = E * A.Center;

% Sample operands and calculate Minkowski product
Esmp = E.sample(30);
Asmp = A.Center + A.Radius * exp(1j*linspace(-pi,pi,30));
EAgsmp = Esmp * Asmp;

%% Plot polar circular product

% Set parameters
cols = [200 165 0]/255;
fontSizeM = 25;
fontSizeL = 30;
cList = getColorList(3);

% Intialize figure
figure(1);clf;
subplot(1,2,1);hold on; axis equal

% Polar intervals
p1 = E.plot('k-','LineWidth',2,'DisplayName','Polar');

% Circular intervals
A.plot('k-','LineWidth',2);
% plot(real(A.Center),imag(A.Center),'+','color',cols(1,:))

% Wrapped intervals
Ar.plot('c-','LineWidth',2);
Er.plot('c-','LineWidth',2);
p2 = Ec.plot('-','LineWidth',2,'Color',cols(1,:),'DisplayName','Circular');

% Product intervals
p3 = EAg.plot('b-','LineWidth',2,'DisplayName','Convex polygon');
p4 = EAa.plot('g-','LineWidth',2,'DisplayName','Concave polyarc');
p5 = EAx.plot('r-','LineWidth',2,'DisplayName','Convex polyarc');
p6 = EAr.plot('c-','LineWidth',2,'DisplayName','Rectangular');
EAc.plot('-','LineWidth',2,'color',cols(1,:));

% Circular interval labels
text(real(A.Center),imag(A.Center), '$A_m$',...
    'fontsize',fontSizeL,'HorizontalAlignment','center','Interpreter','latex')

% Polar interval labels
pnt = E.Abs.mid * exp(1j*E.Angle.mid);
text(real(pnt),imag(pnt), '$E_m$', ...
    'fontsize',fontSizeL,'HorizontalAlignment','center','Interpreter','latex')

% Product interval labels
text(1.6,1.3, '$E_m A_m$', ...
    'fontsize',fontSizeL,'HorizontalAlignment','center','Interpreter','latex')

% Tightness values
annotText = {['$\tau_{EA}^\mathcal{R}=' num2str(100*EAa.Area ./ EAr.Area,3),'\%$'],...
             ['$\tau_{EA}^\mathcal{C}=' num2str(100*EAa.Area ./ EAc.Area,3),'\%$'],...
             ['$\tau_{EA}^\mathcal{G^*}=' num2str(100*EAa.Area ./ EAg.Area,3),'\%$'],...
             ['$\tau_{EA}^\mathcal{A^*}=' num2str(100*EAa.Area ./ EAx.Area,3),'\%$'],...
             ['$\tau_{EA}^\mathcal{A}=' num2str(100*EAa.Area ./ EAa.Area,3),'\%$']};
annotation('textbox',[.133,.35,.105,.21],'String',annotText, ...
           'FontSize',fontSizeM,...
           'BackgroundColor','w',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Misc
xlim([0.55,2.15])
yLim = ylim();
ylim(yLim+0.2)
legend([p2(1),p6(1),p3(1),p5(1),p4(1)],'location','northwest', ...
        'FontSize',fontSizeM)
xlabel('Real')
ylabel('Imag')

%% Calculate polar circular sum-product

% Set intervals
A(2) = ciat.CircularInterval(1-0.5i,0.10);
E(2) = ciat.PolarInterval(1.6 , 2.0 , -0.05*pi, -0.08*pi);

% Wrap intervals in various types
Er = ciat.RectangularInterval(E); 
Ec = ciat.CircularInterval(E);    
Eg = ciat.PolygonalInterval(E,'Tolerance',0.1); 
Ex = ciat.PolyarxInterval(E);
Ea = ciat.PolyarcularInterval(E);
Ar = ciat.RectangularInterval(A); 
Ac = A;    

% Multiply E and A
EAr = Er .* Ar; 
EAc = Ec .* Ac; 
EAg = ciat.PolygonalInterval(E,A,'Tolerance',0.1); 
EAx = ciat.PolyarxInterval(E,A); 
EAa = ciat.PolyarcularInterval(E,A); 

% Sum the EA intervals
Br = sum(EAr);
Bc = sum(EAc);
Bg = sum(EAg);
Bx = sum(EAx);
Ba = sum(EAa);

% Set parameters
Nsmp = 100;

% Sample E and A intervals
Es = E.sample(Nsmp);
As = A.sample(Nsmp);

% Multiply E and A intervals
EAs{1} = Es{1} * As{1}.';
EAs{2} = Es{2} * As{2}.';

% Sum EA intervals
Bs = EAs{1}(:) + EAs{2}(:).';

% Get convex hulls
idx = convhull(real(EAs{1}),imag(EAs{1}));
EAs{1} = EAs{1}(idx);
idx = convhull(real(EAs{2}),imag(EAs{2}));
EAs{2} = EAs{2}(idx);
% idx = convhull(real(Bs),imag(Bs));
% Bs = Bs(idx);



%% Plot polar circular sum-product
% figure(2);clf;
subplot(1,2,2);cla;hold on;axis equal
set(0,'DefaultLineLineWidth',2)

% Plot wraps of E intervals
p1 = E.plot('k-','DisplayName','Polar');
p2 = Er.plot('c-','DisplayName','Rectangular');
p3 = Ec.plot('-','color',cols(1,:),'DisplayName','Circular');

% Plot wraps of A intervals
Ar.plot('c-');
Ac.plot('-','color',cols(1,:));

% Plot origin
plot(0,0,'k+')

% Plot EA intervals
EAr.plot('c-');
EAc.plot('-','color',cols(1,:));
p4 = EAg.plot('b-','DisplayName','Convex polygon');
p5 = EAx.plot('r-','lineWidth',2,'DisplayName','Convex polyarc');
p6 = EAa.plot('g-','lineWidth',2,'DisplayName','Concave polyarc');

% Plot B interval
Br.plot('c-');
Bc.plot('-','color',cols(1,:));
Bg.plot('b');
Bx.plot('r','lineWidth',2);
Ba.plot('g-','lineWidth',2);

% Interval labels
text(real(Ec(1).Center),imag(Ec(1).Center),'$E_1$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Ec(2).Center),imag(Ec(2).Center),'$E_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(EAc(1).Center),imag(EAc(1).Center),'$E_1 A_1$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(EAc(2).Center),imag(EAc(2).Center),'$E_2 A_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Ac(1).Center),imag(Ac(1).Center),'$A_1$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Ac(2).Center),imag(Ac(2).Center),'$A_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Bc.Center),imag(Bc.Center)+0.15,'$B=$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')
text(real(Bc.Center),imag(Bc.Center)-0.15,'$E_1 A_1+E_2 A_2$', ...
    'HorizontalAlignment', 'center','Interpreter','latex')

% Tightness values
annotText = {['$\tau_B^\mathcal{R}=' num2str(100*Ba.Area ./ Br.Area,4),'\%$'],...
             ['$\tau_B^\mathcal{C}=' num2str(100*Ba.Area ./ Bc.Area,4),'\%$'],...
             ['$\tau_B^\mathcal{G^*}=' num2str(100*Ba.Area ./ Bg.Area,4),'\%$'],...
             ['$\tau_B^\mathcal{A^*}=' num2str(100*Ba.Area ./ Bx.Area,4),'\%$'],...
             ['$\tau_B^\mathcal{A}=' num2str(100*Ba.Area ./ Ba.Area,4),'\%$']};
annotation('textbox',[.79,.70,.105,.21],'String',annotText, ...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');


% Plot axes
ylim([-2.5,2.5])
xL = xlim;
yL = ylim;
plot(xL,[0,0],'color',0.7*ones(1,3),'linewidth',0.1)
plot([0,0],yL,'color',0.7*ones(1,3),'linewidth',0.1)

% Misc
fontsize(fontSizeM,'points')
legend([p2(1),p3(1),p4(1),p5(1),p6(1)],'location','southeast')
xlabel('Real')
ylabel('Imag')

annotation('textbox',[0.35,.83,.1,.1],'String','\textbf{a)}', ...
           'FontSize',50,...
           'BackgroundColor','none',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

annotation('textbox',[0.52,.83,.1,.1],'String','\textbf{b)}', ...
           'FontSize',50,...
           'BackgroundColor','none',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');