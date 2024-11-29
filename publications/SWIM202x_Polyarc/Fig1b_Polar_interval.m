clear

% Parameters
    % Sample counts
nr = [10;30];
nR = [4;4];
na = [10;50];
nA = [4;4];
    % Interval dimensions
r = [0.73,0.86;... 
     1.21,1.76];
a  = [0.20,0.25;...
      0.04,0.21]*pi;

% Sample annular sector one
r1 = linspace(r(1,1),r(1,2),nr(1));
R1 = linspace(r(1,1),r(1,2),nR(1));
a1 = linspace(a(1,1),a(1,2),na(1));
A1 = linspace(a(1,1),a(1,2),nA(1));
int1 = [reshape(r1'*exp(1j*A1),[],1);...
        reshape(R1'*exp(1j*a1),[],1)];
Int1 = reshape(R1'*exp(1j*A1),[],1);
Bnd1 = [r(1,1)*exp(1j*a1), ...
        r(1,2)*exp(1j*flip(a1)),...
        r(1,1)*exp(1j*a(1,1))].';
Cnt1 = mean(r1)*exp(1j*mean(a1));

% Sample annular sector two
r2 = linspace(r(2,1),r(2,2),nr(2));
R2 = linspace(r(2,1),r(2,2),nR(2));
a2 = linspace(a(2,1),a(2,2),na(2));
A2 = linspace(a(2,1),a(2,2),nA(2));
int2 = [reshape(r2'*exp(1j*A2),[],1);...
        reshape(R2'*exp(1j*a2),[],1)];
Int2 = reshape(R2'*exp(1j*A2),[],1);
Bnd2 = reshape([r(2,1)*exp(1j*a2) ...
                r(2,2)*exp(1j*flip(a2)),...
                r(2,1)*exp(1j*a(2,1))],[],1);
Cnt2 = mean(r2)*exp(1j*mean(a2));

% Calculate sum
ISum = (Int1.' + Int2).';
iSum = (int1.' + int2).';
iSum1 = (int1.' + Int2).';
iSum2 = (int2.' + Int1).';
bSum1 = (Int2.' + Bnd1).';
bSum2 = (Int1.' + Bnd2).';
bSum = boundary(real(iSum(:)),imag(iSum(:)),0.3);
bSum = iSum(bSum);
pSum = convhull(real(iSum(:)),imag(iSum(:)));
pSum = iSum(pSum);
cSum = Cnt1 + Int2;

% Calculate sum bounding rectangle
rSum = [min(abs(bSum)),max(abs(bSum))];
aSum = linspace(min(angle(bSum)),max(angle(bSum)),na(2));
rSum = [rSum(1)*exp(1j*aSum), ...
        rSum(2)*exp(1j*flip(aSum)),...
        rSum(1)*exp(1j*aSum(1))].';

% Calculate product
IPrd = Int1 * Int2.';
iPrd1 = Int1 * int2.';
iPrd2 = Int2 * int1.';
bPrd1 = Int2 * Bnd1.';
bPrd2 = Int1 * Bnd2.';
bPrd = boundary(real(iPrd2(:)),imag(iPrd2(:)),0.5);
bPrd = iPrd2(bPrd);
gPrd = convhull(real(iPrd2(:)),imag(iPrd2(:)));
gPrd = iPrd2(gPrd);
cPrd = Cnt1 * Int2;

% Calculate inverse
iInv = 1./int2;
IInv = 1./Int2;
bInv = boundary(real(iInv),imag(iInv),0.7);
bInv = iInv(bInv);
pInv = convhull(real(iInv),imag(iInv));
pInv = iInv(pInv);


% Interval arithmetic
A_x = ciat.PolarInterval(r(1,1),r(1,2),a(1,1),a(1,2));
B_x = ciat.PolarInterval(r(2,1),r(2,2),a(2,1),a(2,2));
D_x = A_x * B_x;
A_g = ciat.PolygonalInterval(A_x,'tolerance',0.01);
B_g = ciat.PolygonalInterval(B_x,'tolerance',0.01);
C_g = A_g + B_g;
D_g = A_g * B_g;
A_a = ciat.PolyarcularInterval(A_x);
B_a = ciat.PolyarcularInterval(B_x);
C_a = A_a + B_a;
D_a = ciat.PolyarcularInterval(ciat.PolygonalInterval(bPrd));% D_a = A_a * B_a;
E_a = recip(B_a);
C_x = ciat.PolarInterval(C_a);
E_x = ciat.PolarInterval(E_a);
E_g = ciat.PolygonalInterval(E_a);

%% Plot

figure(2);clf;hold on;axis equal;grid on

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
    % PS2 = scatter(real(ISum(:,idx)),imag(ISum(:,idx)),dotSize(2),[col,'o'], ...
    %                             'filled','DisplayName','A center + B anchors');
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
pI = E_a.plot('k-','LineWidth',lineWidth,'DisplayName','Result in polyarc');
pG = E_g.plot('k--','LineWidth',lineWidth,'DisplayName','Result in polygon');
pX = E_x.plot('k:','LineWidth',lineWidth,'DisplayName','Result in polar');
% Adjust x-limits
xLim = xlim();
xlim(xLim-0.3)



% Axes
xlabel('Real','FontSize',fontSize,'Position',[-1.3,-0.52,-1.26])
ylabel('Imag','FontSize',fontSize,'Position',[-1.45,-0.35,-1.27])
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
text(0.64,0.43, '$A$','FontSize',fontSize,'horizontalAlignment','Left','Interpreter','latex')
text(1.33,0.56, '$B$','FontSize',fontSize,'horizontalAlignment','Left','Interpreter','latex')
text(2.4,1.11, '$C\!=\!A+B$','FontSize',fontSize,'horizontalAlignment','Left','Interpreter','latex')
text(0.70,1.40, '$D\!=\!A\times B$','FontSize',fontSize,'horizontalAlignment','left','Interpreter','latex')
text(0.33,-0.15,'$E\!=\!1/B$','FontSize',fontSize,'horizontalAlignment','center','Interpreter','latex')


% Tightness values
annotText = {['$\tau(C^\mathcal{A})=' num2str(100*C_a.Area ./ C_a.Area,3),'\%$'],...
             ['$\tau(C^\mathcal{G})=' num2str(100*C_a.Area ./ C_g.Area,3),'\%$'],...
             ['$\tau(C^\mathcal{P})=' num2str(100*C_a.Area ./ C_x.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{A})=' num2str(100*D_a.Area ./ D_a.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{G})=' num2str(100*D_a.Area ./ D_g.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{P})=' num2str(100*D_a.Area ./ D_x.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{A})=' num2str(100*E_a.Area ./ E_a.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{G})=' num2str(100*E_a.Area ./ E_g.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{P})=' num2str(100*E_a.Area ./ E_x.Area,3),'\%$']};
annotation('textbox',[.77,.07,.15,.55],'String',annotText, ...
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

annotation('textbox',[.18,.83,.1,.1],'String','\textbf{b)}', ...
           'FontSize',50,...
           'BackgroundColor','w',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');