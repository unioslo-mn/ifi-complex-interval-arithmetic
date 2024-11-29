





clear

% Parameters
    % Sample counts
nr = [10;20];
nR = [4;4];
na = [40;60];
nA = [7;13];
    % Interval dimensions
r = [0.10;... 
     0.65];
offs = [0.9+1.1i;...
        1.8+0.6i];

 
% Sample annular sector one
r1 = linspace(0,r(1),nr(1));
R1 = linspace(0,r(1),nR(1));
a1 = linspace(-pi,pi,na(1));
A1 = linspace(-pi,pi,nA(1));
int1 = [reshape(r1'*exp(1j*A1),[],1);...
        reshape(R1'*exp(1j*a1),[],1)] + offs(1);
Int1 = reshape(R1'*exp(1j*A1),[],1) + offs(1);
Bnd1 = (r(1)*exp(1j*a1) + offs(1)).';
Cnt1 = offs(1);

% Sample annular sector two
r2 = linspace(0,r(2),nr(2));
R2 = linspace(0,r(2),nR(2));
a2 = linspace(-pi,pi,na(2));
A2 = linspace(-pi,pi,nA(2));
int2 = [reshape(r2'*exp(1j*A2),[],1);...
        reshape(R2'*exp(1j*a2),[],1)] + offs(2);
Int2 = reshape(R2'*exp(1j*A2),[],1) + offs(2);
Bnd2 = (r(2)*exp(1j*a2) + offs(2)).';
Cnt2 = offs(2);

% Calculate sum
ISum = (Int1.' + Int2).';
iSum = (int1.' + int2).';
iSum1 = (int1.' + Int2).';
iSum2 = (int2.' + Int1).';
bSum1 = (Int2.' + Bnd1).';
bSum2 = (Int1.' + Bnd2).';
bSum = convhull(real(iSum(:)),imag(iSum(:)));
bSum = iSum(bSum);
pSum = convhull(real(iSum(:)),imag(iSum(:)));
pSum = iSum(pSum);
cSum = Cnt1 + Int2;

% Calculate product
IPrd = Int1 * Int2.';
iPrd = int1 * int2.';
iPrd1 = Int1 * int2.';
iPrd2 = Int2 * int1.';
bPrd1 = Int2 * Bnd1.';
bPrd2 = Int1 * Bnd2.';
bPrd = convhull(real(iPrd(:)),imag(iPrd(:)));
bPrd = iPrd(bPrd);
% gPrd = convhull(real(iPrd2(:)),imag(iPrd2(:)));
% gPrd = iPrd2(gPrd);
cPrd = Cnt1 * Int2;

% Calculate inverse
iInv = 1./int2;
IInv = 1./Int2;
bInv = convhull(real(iInv),imag(iInv));
bInv = iInv(bInv);
pInv = convhull(real(iInv),imag(iInv));
pInv = iInv(pInv);
zInv = 1/offs(2);
rInv = (r(2)*abs(zInv)^2)/( (r(2)*abs(zInv))^2 - 1 );
oInv = - abs(zInv) * cos(angle(zInv)) / ( (r(2)*abs(zInv))^2 - 1 ) - ...
    1j * abs(zInv) * sin(angle(zInv)) / ( (r(2)*abs(zInv))^2 - 1 );

% Interval arithmetic
sCnt = 100;
A_x = ciat.CircularInterval(offs(1),r(1));
B_x = ciat.CircularInterval(offs(2),r(2));
C_x = A_x + B_x;
D_x = A_x * B_x;
A_g = ciat.PolygonalInterval(A_x,'tolerance',0.01);
B_g = ciat.PolygonalInterval(B_x,'tolerance',0.01);
C_g = A_g + B_g;
D_g = A_g * B_g;
E_g = ciat.PolygonalInterval(1./B_g.Points);
A_a = ciat.PolyarcularInterval(A_x);
B_a = ciat.PolyarcularInterval(B_x);
C_a = A_a + B_a;
D_a = ciat.PolyarcularInterval(ciat.PolygonalInterval(bPrd));% D_a = A_a * B_a;
E_a = recip(B_a);
E_x = ciat.CircularInterval(E_a);


%% Plot

figure(3);clf;hold on;axis equal;grid on

set(gca, 'Position', [0.0,0.06,1,0.9]);

% Parameters
lineWidth = 3;
fontSize = 30;
dotSize = [5 15 20]*2;



% Interval 1
pA1 = plot(real(Bnd1),imag(Bnd1),'b-','LineWidth',lineWidth, ...
                                    'DisplayName','A interval');
pA2 = scatter(real(int1),imag(int1),dotSize(1),'bo','filled', ...
                                    'DisplayName','A samples');
pA3 = scatter(real(Cnt1),imag(Cnt1),dotSize(3),'bo',...
                                    'DisplayName','A center');

% Interval 2
pB1 = plot(real(Bnd2),imag(Bnd2),'r-','LineWidth',lineWidth, ...
                                    'DisplayName','B interval');
pB2 = scatter(real(int2),imag(int2),dotSize(1),'ro','filled', ...
                                    'DisplayName','B samples');
pB3 = scatter(real(Int2),imag(Int2),dotSize(3),'ro',...
                                    'DisplayName','B anchors');


% Sum
pS1 = scatter(real(cSum),imag(cSum),dotSize(3),'ms', ...
                                'DisplayName','A center + B anchors');
for idx = 1:size(bSum1,1)
    if mod(idx + fix((idx-1)/4),2)
       col = 'c';
    else
       col = 'm';
    end
    pS2 = scatter(real(iSum1(:,idx)),imag(iSum1(:,idx)),dotSize(1),[col,'s'], ...
                                'filled','DisplayName','A samples + B anchors');
    pS3 = plot(real(bSum1(idx,:)),imag(bSum1(idx,:)),[col,'-'],'LineWidth',lineWidth,...
                                 'DisplayName','A interval + B anchors');
end
C_a.plot('k-','LineWidth',lineWidth);
C_g.plot('k--','LineWidth',lineWidth);
C_x.plot('k:','LineWidth',lineWidth);

% Product
pP1 = scatter(real(cPrd),imag(cPrd),dotSize(3),'md', ...
                                'DisplayName','A center \times B anchors');
for idx = 1:size(bPrd1,1)
    if mod(idx + fix((idx-1)/4),2)
       col = 'm';
    else
       col = 'c';
    end
    pP2 = scatter(real(iPrd2(idx,:)),imag(iPrd2(idx,:)),dotSize(1),[col,'d'], ...
                                'filled','DisplayName','A samples \times B anchors');
    % scatter(real(IPrd(:,idx)),imag(IPrd(:,idx)),dotSize(2),[col,'o'],'filled');
    pP3 = plot(real(bPrd1(idx,:)),imag(bPrd1(idx,:)),[col,'-'],'LineWidth',lineWidth,...
                                 'DisplayName','A interval \times B anchors');
end
D_a.plot('k-','LineWidth',lineWidth);
D_g.plot('k--','LineWidth',lineWidth);
D_x.plot('k:','LineWidth',lineWidth);

% Reciprocal
pR1 = scatter(real(iInv),imag(iInv),dotSize(1),'rs','filled',...
                        'DisplayName','1 / B samples');
% scatter(real(IInv),imag(IInv),dotSize(3),'ro')
pI = E_a.plot('k-','LineWidth',lineWidth,'DisplayName','Result in polyarc');
pG = E_g.plot('k--','LineWidth',lineWidth,'DisplayName','Result in polygon');
pX = E_x.plot('k:','LineWidth',lineWidth,'DisplayName','Result in circle');

% Adjust x-limits
xLim = xlim();
xlim(xLim-0.5)


% Axes
xlabel('Real','FontSize',fontSize,'Position',[-3.2,-0.4,-1.26])
ylabel('Imag','FontSize',fontSize,'Position',[-3.5,-0.2,-1.27])
fontsize(30,'points')
xTicks = xticks;
xticks(xTicks(2:end))
yTicks = yticks;
yticks(yTicks(2:end))
xl = xlim();
yl = ylim();
line(xl,[0,0],'Color','black');
line([0,0],yl,'Color','black');


% Texts
text(0.65,1.08, '$A$','FontSize',fontSize,'horizontalAlignment','Center','Interpreter','latex')
text(1,0.5,     '$B$','FontSize',fontSize,'horizontalAlignment','Center','Interpreter','latex')
text(2.70,2.56, '$C\!=\!A+B$','FontSize',fontSize,'horizontalAlignment','center','Interpreter','latex')
text(2,3.2,     '$D\!=\!A\times B$','FontSize',fontSize,'horizontalAlignment','left','Interpreter','latex')
text(0.25,-0.19,'$E\!=\!1/B$','FontSize',fontSize,'horizontalAlignment','right','Interpreter','latex')

% Tightness values
annotText = {['$\tau(C^\mathcal{A})=' num2str(100*C_a.Area ./ C_a.Area,3),'\%$'],...
             ['$\tau(C^\mathcal{G})=' num2str(100*C_a.Area ./ C_g.Area,3),'\%$'],...
             ['$\tau(C^\mathcal{C})=' num2str(100*C_a.Area ./ C_x.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{A})=' num2str(100*D_a.Area ./ D_a.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{G})=' num2str(100*D_a.Area ./ D_g.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{C})=' num2str(100*D_a.Area ./ D_x.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{A})=' num2str(100*E_a.Area ./ E_a.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{G})=' num2str(100*E_a.Area ./ E_g.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{C})=' num2str(100*E_a.Area ./ E_x.Area,3),'\%$']};
annotation('textbox',[.77,.40,.15,.55],'String',annotText, ...
           'FontSize',fontSize,...
           'BackgroundColor','w',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Legend
legend([pA1(1),pA2(1),pA3(1), ...
        pB1(1),pB2(1),pB3, ...
        pS1(1),pS2(1),pS3(1), ...
        pP1(1),pP2(1),pP3(1),...
        pR1(1),pI(1),pG(1),pX(1)], ...
        'Location','northwest','FontSize',24)


annotation('textbox',[.18,.83,.1,.1],'String','\textbf{c)}', ...
           'FontSize',50,...
           'BackgroundColor','w',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');