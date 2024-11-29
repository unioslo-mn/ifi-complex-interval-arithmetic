clear
% close all; figure();



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
pPrd = convhull(real(iPrd),imag(iPrd));
pPrd = iPrd(pPrd);
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


%% Plot

subplot(1,1,1);cla;hold on
set(gca, 'Position', [0.04,0.06,0.95,0.9]);

% Parameters
lineWidth = 2;
fontSize = 18;
dotSize = [5 15 20];
dotColor = [0 0 1; 1 0 0];

title('Rectangular interval arithmetic properties','FontSize',fontSize)

    % Interval 1
scatter(real(int1),imag(int1),dotSize(1),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1)
scatter(real(Int1),imag(Int1),dotSize(2),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1)
plot(real(Bnd1),imag(Bnd1),'b-','LineWidth',lineWidth)
    % Interval 2
scatter(real(int2),imag(int2),dotSize(1),dotColor(2,:),'o','filled')%,'MarkerFaceAlpha',1)
scatter(real(Int2),imag(Int2),dotSize(3),dotColor(2,:),'o')%,'MarkerFaceAlpha',1)
plot(real(Bnd2),imag(Bnd2),'r-','LineWidth',lineWidth)

    % Sum
scatter(real(cSum),imag(cSum),dotSize(3),dotColor(2,:),'o');
for idx = 1:size(bSum1,1)
    p1 = scatter(real(iSum1(:,idx)),imag(iSum1(:,idx)),dotSize(1),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    p2 = scatter(real(ISum(:,idx)),imag(ISum(:,idx)),dotSize(2),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    lh = plot(real(bSum1(idx,:)),imag(bSum1(idx,:)),'-','LineWidth',lineWidth);
    if mod(idx + fix((idx-1)/4),2)
       lh.Color = [0,0,1,1];
       p1.MarkerFaceColor = [0,0.2,1];
       p2.MarkerFaceColor = [0,0.2,1];
    else
       lh.Color = [0,0.8,1,1];
       p1.MarkerFaceColor = [0,0.8,1];
       p2.MarkerFaceColor = [0,0.8,1];
    end
end
plot(real(bSum),imag(bSum),'k-','LineWidth',lineWidth)

    % Product
scatter(real(cPrd),imag(cPrd),dotSize(3),dotColor(2,:),'o');
for idx = 1:size(bPrd1,1)
    p1 = scatter(real(iPrd2(idx,:)),imag(iPrd2(idx,:)),dotSize(1),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    p2 = scatter(real(IPrd(:,idx)),imag(IPrd(:,idx)),dotSize(2),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    lh = plot(real(bPrd1(idx,:)),imag(bPrd1(idx,:)),'-','LineWidth',lineWidth);
    if mod(idx + fix((idx-1)/4),2)
       lh.Color = [0,0,1,1];
       p1.MarkerFaceColor = [0,0.2,1];
       p2.MarkerFaceColor = [0,0.2,1];
    else
       lh.Color = [0,0.8,1,1];
       p1.MarkerFaceColor = [0,0.8,1];
       p2.MarkerFaceColor = [0,0.8,1];
    end
end
plot(real(bPrd),imag(bPrd),'k-','LineWidth',lineWidth)
plot(real(rPrd),imag(rPrd),'k--','LineWidth',lineWidth)
plot(real(pPrd),imag(pPrd),'k:','LineWidth',lineWidth)

    % Reciprocal
scatter(real(iInv),imag(iInv),dotSize(1),dotColor(2,:),'o','filled')
scatter(real(IInv),imag(IInv),dotSize(3),dotColor(2,:),'o')
plot(real(bInv),imag(bInv),'k-','LineWidth',lineWidth)
plot(real(rInv),imag(rInv),'k--','LineWidth',lineWidth)
plot(real(pInv),imag(pInv),'k:','LineWidth',lineWidth)

    % Axes
axis equal;grid on
xlabel('Real','FontSize',fontSize)
ylabel('Imaginary','FontSize',fontSize)
xl = xlim();
yl = ylim();
line(xl,[0,0],'Color','black');
line([0,0],yl,'Color','black');


    % Texts
text(1.02,0.61,'A','FontSize',20)
text(1.69,0.66,'B','FontSize',20)
text(2.70,1.50,'A+B','FontSize',20)
text(0.90,2.07,'A×B','FontSize',20)
text(0.47,-0.05,'1/B','FontSize',20)





