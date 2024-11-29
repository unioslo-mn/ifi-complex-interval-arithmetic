%% Rectangular interval

clear

% Parameters
    % Sample counts
nx = [10;30];
nX = [4;4];
ny = [10;30];
nY = [4;4];
    % Interval dimensions
x = [1.15,0.95;... 
     1.35,2.10];
y  = [0.40,0.55;...
      0.4,0.9];
  
% Sample rectangle sector one
x1 = linspace(x(1,1),x(1,2),nx(1));
X1 = linspace(x(1,1),x(1,2),nX(1));
y1 = linspace(y(1,1),y(1,2),ny(1));
Y1 = linspace(y(1,1),y(1,2),nY(1));
int1 = [reshape(x1' + 1j*Y1,[],1);...
        reshape(X1' + 1j*y1,[],1)];
Int1 = reshape(X1' + 1j*Y1,[],1);
Bnd1 = [x1 + 1j*y1(1),...
        x1(end) + 1j*y1,...
        flip(x1) + 1j*y1(end),...
        x1(1) + flip(1j*y1)].';
Cnt1 = mean(x1)+1j*mean(y1);
            
% Sample rectangle sector two
x2 = linspace(x(2,1),x(2,2),nx(2));
X2 = linspace(x(2,1),x(2,2),nX(2));
y2 = linspace(y(2,1),y(2,2),ny(2));
Y2 = linspace(y(2,1),y(2,2),nY(2));
int2 = [reshape(x2' + 1j*Y2,[],1);...
        reshape(X2' + 1j*y2,[],1)];
Int2 = reshape(X2' + 1j*Y2,[],1);
Bnd2 = [x2 + 1j*y2(1),...
        x2(end) + 1j*y2,...
        flip(x2) + 1j*y2(end),...
        x2(1) + flip(1j*y2)].';

% Calculate sum
ISum = (Int1.' + Int2).';
iSum = (int1.' + int2).';
iSum1 = (int1.' + Int2).';
iSum2 = (int2.' + Int1).';
bSum1 = (Int2.' + Bnd1).';
bSum2 = (Int1.' + Bnd2).';
bSum = convhull(real(iSum(:)),imag(iSum(:)));
bSum = iSum(bSum);
cSum = Cnt1 + Int2;

% Calculate product
IPrd = Int1 * Int2.';
iPrd = int1 * int2.';
iPrd1 = Int1 * int2.';
iPrd2 = Int2 * int1.';
bPrd1 = Int2 * Bnd1.';
bPrd2 = Int1 * Bnd2.';
bPrd = boundary(real(iPrd(:)),imag(iPrd(:)),0.5);
bPrd = iPrd(bPrd);

% Calculate product bounding rectangle
xPrd = [min(real(bPrd)),max(real(bPrd))];
yPrd = [min(imag(bPrd)),max(imag(bPrd))];
rPrd = [xPrd + 1j*yPrd(1),...
        xPrd(end) + 1j*yPrd,...
        flip(xPrd) + 1j*yPrd(end),...
        xPrd(1) + flip(1j*yPrd)];
gPrd = convhull(real(iPrd),imag(iPrd));
gPrd = iPrd(gPrd);
cPrd = Cnt1 * Int2;

% Calculate inverse
iInv = 1./int2;
IInv = 1./Int2;
bInv = boundary(real(iInv),imag(iInv),0.5);
bInv = iInv(bInv);
pInv = convhull(real(iInv),imag(iInv));
pInv = iInv(pInv);

% Calculate inverse bounding rectangle
xInv = [min(real(bInv)),max(real(bInv))];
yInv = [min(imag(bInv)),max(imag(bInv))];
rInv = [xInv + 1j*yInv(1),...
        xInv(end) + 1j*yInv,...
        flip(xInv) + 1j*yInv(end),...
        xInv(1) + flip(1j*yInv)];


% Interval arithmetic
A_x = ciat.RectangularInterval(x(1,1),x(1,2),y(1,1),y(1,2));
B_x = ciat.RectangularInterval(x(2,1),x(2,2),y(2,1),y(2,2));
C_x = A_x + B_x;
A_g = ciat.PolygonalInterval(A_x,'tolerance',0.01);
B_g = ciat.PolygonalInterval(B_x,'tolerance',0.01);
C_g = A_g + B_g;
D_g = A_g * B_g;
A_a = ciat.PolyarcularInterval(A_x);
B_a = ciat.PolyarcularInterval(B_x);
C_a = A_a + B_a;
D_a = ciat.PolyarcularInterval(ciat.PolygonalInterval(bPrd));% D_a = A_a * B_a;
E_a = recip(B_a);% E_a = ciat.PolyarcularInterval(E_a.DefArcs([])    )
D_x = ciat.RectangularInterval(D_a);
E_x = ciat.RectangularInterval(E_a);
E_g = ciat.PolygonalInterval(E_a,'tolerance',0.01);

%% Plot

figure(1);clf;hold on;axis equal;grid on

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
% plot(real(bSum),imag(bSum),'k-','LineWidth',lineWidth);
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
% plot(real(bPrd),imag(bPrd),'k-','LineWidth',lineWidth)
% plot(real(rPrd),imag(rPrd),'k--','LineWidth',lineWidth)
% plot(real(gPrd),imag(gPrd),'k:','LineWidth',lineWidth)
D_a.plot('k-','LineWidth',lineWidth);
D_g.plot('k--','LineWidth',lineWidth);
D_x.plot('k:','LineWidth',lineWidth);

% Reciprocal
pR1 = scatter(real(iInv),imag(iInv),dotSize(1),'rs','filled',...
                        'DisplayName','1 / B samples');
% scatter(real(IInv),imag(IInv),dotSize(3),'ro')
% pI = plot(real(bInv),imag(bInv),'k-','LineWidth',lineWidth,...
%                             'DisplayName','Result bounds');
% pG = plot(real(rInv),imag(rInv),'k--','LineWidth',lineWidth,...
%                             'DisplayName','Result in convex polygon');
% pX = plot(real(pInv),imag(pInv),'k:','LineWidth',lineWidth,...
%                             'DisplayName','Result in rectangle');
pI = E_a.plot('k-','LineWidth',lineWidth,'DisplayName','Result in polyarc');
pG = E_g.plot('k--','LineWidth',lineWidth,'DisplayName','Result in polygon');
pX = E_x.plot('k:','LineWidth',lineWidth,'DisplayName','Result in rectangle');

% Adjust x-limits
xLim = xlim();
xlim(xLim-1)


% Axes
xlabel('Real','FontSize',fontSize,'Position',[-1.85,-0.36,-1.26])
ylabel('Imag','FontSize',fontSize,'Position',[-2.0,-0.26,-1.27])
xTicks = xticks;
xticks(xTicks(2:end))
xl = xlim();
yl = ylim();
line(xl,[0,0],'Color','black');
line([0,0],yl,'Color','black');


% Texts
text(1.02,0.61, '$A$','FontSize',fontSize,'horizontalAlignment','Center','Interpreter','latex')
text(1.69,0.66, '$B$','FontSize',fontSize,'horizontalAlignment','Center','Interpreter','latex')
text(2.70,1.50, '$C\!=\!A+B$','FontSize',fontSize,'horizontalAlignment','center','Interpreter','latex')
text(0.90,2.07, '$D\!=\!A\times B$','FontSize',fontSize,'horizontalAlignment','left','Interpreter','latex')
text(0.75,-0.2,'$E\!=\!1/B$','FontSize',fontSize,'horizontalAlignment','left','Interpreter','latex')

%Tightness values
annotText = {['$\tau(C^\mathcal{A})=' num2str(100*C_a.Area ./ C_a.Area,3),'\%$'],...
             ['$\tau(C^\mathcal{G})=' num2str(100*C_a.Area ./ C_g.Area,3),'\%$'],...
             ['$\tau(C^\mathcal{R})=' num2str(100*C_a.Area ./ C_x.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{A})=' num2str(100*D_a.Area ./ D_a.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{G})=' num2str(100*D_a.Area ./ D_g.Area,3),'\%$'],...
             ['$\tau(D^\mathcal{R})=' num2str(100*D_a.Area ./ D_x.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{A})=' num2str(100*E_a.Area ./ E_a.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{G})=' num2str(100*E_a.Area ./ E_g.Area,3),'\%$'],...
             ['$\tau(E^\mathcal{R})=' num2str(100*E_a.Area ./ E_x.Area,3),'\%$']};
annotation('textbox',[.30,.40,.15,.55],'String',annotText, ...
           'FontSize',fontSize,...
           'BackgroundColor','w',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Set font size
fontsize(30,'points')

annotation('textbox',[.18,.83,.1,.1],'String','\textbf{a)}', ...
           'FontSize',50,...
           'BackgroundColor','w',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

% Legend
legend([pA1(1),pA2(1),pA3(1), ...
        pB1(1),pB2(1),pB3, ...
        pS1(1),pS2(1),pS3(1), ...
        pP1(1),pP2(1),pP3(1),...
        pR1(1),pI(1),pG(1),pX(1)], ...
        'Location','northwest','FontSize',24)