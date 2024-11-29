clear
% close all;figure();
cla;hold on

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
pPrd = convhull(real(iPrd2(:)),imag(iPrd2(:)));
pPrd = iPrd2(pPrd);
cPrd = Cnt1 * Int2;

% Calculate inverse
iInv = 1./int2;
IInv = 1./Int2;
bInv = boundary(real(iInv),imag(iInv),0.7);
bInv = iInv(bInv);
pInv = convhull(real(iInv),imag(iInv));
pInv = iInv(pInv);

%% Plot

subplot(1,1,1);cla;hold on
set(gca, 'Position', [0.04,0.06,0.95,0.9]);

% Parameters
lineWidth = 2;
fontSize = 18;
dotSize = [5 15 20];
dotColor = [0 0 1; 1 0 0];

title('Polar interval arithmetic properties','FontSize',fontSize)


    % Interval 1
scatter(real(int1),imag(int1),dotSize(1),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1)
scatter(real(Int1),imag(Int1),dotSize(2),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1)
plot(real(Bnd1),imag(Bnd1),'b-','LineWidth',lineWidth)

    % Interval 2
scatter(real(int2),imag(int2),dotSize(1),dotColor(2,:),'o','filled')%,'MarkerFaceAlpha',0.5)
scatter(real(Int2),imag(Int2),dotSize(3),dotColor(2,:),'o')%,'MarkerFaceAlpha',0.5)
plot(real(Bnd2),imag(Bnd2),'r-','LineWidth',lineWidth)

        % Sum
scatter(real(cSum),imag(cSum),20,'ro');
for idx = 1:size(bSum1,1)
    p1 = scatter(real(iSum1(:,idx)),imag(iSum1(:,idx)),dotSize(1),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1);
    p2 = scatter(real(ISum(:,idx)),imag(ISum(:,idx)),dotSize(2),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1);
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
plot(real(pSum),imag(pSum),'k:','LineWidth',lineWidth)
plot(real(rSum),imag(rSum),'k--','LineWidth',lineWidth)

    % Product
scatter(real(cPrd),imag(cPrd),20,'ro');
for idx = 1:size(bPrd1,1)
    p1 = scatter(real(iPrd2(idx,:)),imag(iPrd2(idx,:)),dotSize(1),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1);
    p2 = scatter(real(IPrd(:,idx)),imag(IPrd(:,idx)),dotSize(2),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1);
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
plot(real(pPrd),imag(pPrd),'k:','LineWidth',lineWidth)

    % Inverse
scatter(real(iInv),imag(iInv),dotSize(1),dotColor(2,:),'o','filled')
scatter(real(IInv),imag(IInv),20,'ro')
plot(real(bInv),imag(bInv),'k-','LineWidth',lineWidth)
plot(real(pInv),imag(pInv),'k:','LineWidth',lineWidth)

    % Axes
axis equal;grid on
xlabel('Real','FontSize',fontSize)
ylabel('Imaginary','FontSize',fontSize)
xl = xlim();
yl = ylim();
line(xl,[0,0],'Color','black');
line([0,0],yl,'Color','black');
axis equal;grid on

    % Texts
text(0.64,0.43,'A','FontSize',fontSize)
text(1.33,0.56,'B','FontSize',fontSize)
text(1.92,1.11,'A+B','FontSize',fontSize)
text(0.70,1.40,'A×B','FontSize',fontSize)
text(0.38,-0.15,'1/B','FontSize',fontSize)