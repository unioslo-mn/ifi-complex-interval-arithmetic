clear
% close all

%% Set parameters

m = [6+12i ; -2-5i ; -3+10i];
d = ciat.RealInterval([11 8 5],[12 10 7])';
a = ciat.RealInterval([14 -147 63],[100 -75 150])'/180*pi;

%% Generate intervals

F_p = ciat.PolarInterval(d,a);
F_a = ciat.PolyarcularInterval(F_p);
iF_a = -F_a + m;
L_a = intersection(iF_a);
L_c = ciat.CircularInterval(L_a);
Loc = L_c.Center;
Msr = F_a + Loc;


%% Plot
figure(1);clf;hold on;axis equal

% Plot reference points
p1 = scatter(real(m),imag(m),100,'ks','filled','DisplayName','Landmark');
text(real(m(1))+0.3,imag(m(1))+0.5,sprintf('$P_%i$',1), ...
        'HorizontalAlignment','left','Interpreter','latex')
text(real(m(2))+0.2,imag(m(2)),sprintf('$P_%i$',2), ...
        'HorizontalAlignment','left','Interpreter','latex')
text(real(m(3))+0.2,imag(m(3)),sprintf('$P_%i$',3), ...
        'HorizontalAlignment','left','Interpreter','latex')


% Plot solution
p4 = L_a.plot('k-','linewidth',3,'DisplayName','Solution set (intersection)');
p3 = iF_a.plot('k--','DisplayName','Inverted measurement');


% Plot measurements
p5 = scatter(real(Loc),imag(Loc),100,'ko','filled','DisplayName','Center of solution set');
p2 = Msr.plot('k-','DisplayName','Measurement (relative to set center)');

% Interval labels
text(10,10,'$A_1$', 'HorizontalAlignment', 'center','Interpreter','latex')
text(-5,-2,'$A_2$', 'HorizontalAlignment', 'center','Interpreter','latex')
text(-3,12,'$A_3$', 'HorizontalAlignment', 'center','Interpreter','latex')
text(8,2,'$P_1-A_1$', 'HorizontalAlignment', 'left','Interpreter','latex')
text(4,3,'$P_2-A_2$', 'HorizontalAlignment', 'left','Interpreter','latex')
text(3,8,'$P_3-A_3$', 'HorizontalAlignment', 'center','Interpreter','latex')
text(-2,2,'$B$', 'HorizontalAlignment', 'center','Interpreter','latex')

fontsize(30,'points')
legend([p1(1),p2(1),p3(1),p4(1),p5(1)],'Location','northwest')
xlabel('Real (x)')
ylabel('Imag (y)')
xLim = xlim();
xlim(xLim-7)