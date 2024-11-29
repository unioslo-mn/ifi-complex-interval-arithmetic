clear
% close all;figure();
cla;hold on

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

% Sample annular sector two
r2 = linspace(0,r(2),nr(2));
R2 = linspace(0,r(2),nR(2));
a2 = linspace(-pi,pi,na(2));
A2 = linspace(-pi,pi,nA(2));
int2 = [reshape(r2'*exp(1j*A2),[],1);...
        reshape(R2'*exp(1j*a2),[],1)] + offs(2);
Int2 = reshape(R2'*exp(1j*A2),[],1) + offs(2);
Bnd2 = (r(2)*exp(1j*a2) + offs(2)).';

% Calculate sum
ISum = (Int1.' + Int2).';
iSum = (int1.' + int2).';
iSum1 = (int1.' + Int2).';
iSum2 = (int2.' + Int1).';
bSum1 = (Int2.' + Bnd1).';
bSum2 = (Int1.' + Bnd2).';
bSum = convhull(real(iSum(:)),imag(iSum(:)));
bSum = iSum(bSum);

% Calculate product
IPrd = Int1 * Int2.';
iPrd = int1 * int2.';
iPrd1 = Int1 * int2.';
iPrd2 = Int2 * int1.';
bPrd1 = Int2 * Bnd1.';
bPrd2 = Int1 * Bnd2.';
bPrd = convhull(real(iPrd(:)),imag(iPrd(:)));
bPrd = iPrd(bPrd);

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


%% Plot

subplot(1,1,1);cla;hold on
set(gca, 'Position', [0.04,0.06,0.95,0.9]);

% Parameters
lineWidth = 2;
fontSize = 18;
dotSize = [5 15 20];
dotColor = [0 0 1; 1 0 0];

title('Circular interval arithmetic properties','FontSize',fontSize)

    % Interval 1
scatter(real(int1),imag(int1),dotSize(1),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1)
scatter(real(Int1),imag(Int1),dotSize(2),dotColor(1,:),'o','filled')%,'MarkerFaceAlpha',1)
plot(real(Bnd1),imag(Bnd1),'b-','LineWidth',2)

    % Interval 2
scatter(real(int2),imag(int2),dotSize(1),dotColor(2,:),'o','filled')%,'MarkerFaceAlpha',0.5)
scatter(real(Int2),imag(Int2),dotSize(2),dotColor(2,:),'o','filled')%,'MarkerFaceAlpha',0.5)
plot(real(Bnd2),imag(Bnd2),'r-','LineWidth',2)

    % Sum
%scatter(real(iSum2),imag(iSum2),2,dotColor(2,:),'o','filled')%,'MarkerFaceAlpha',0.5);
for idx = 1:size(bSum1,1)
    p1 = scatter(real(iSum1(:,idx)),imag(iSum1(:,idx)),dotSize(1),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    p2 = scatter(real(ISum(:,idx)),imag(ISum(:,idx)),dotSize(2),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    lh = plot(real(bSum1(idx,:)),imag(bSum1(idx,:)),'-','LineWidth',2);
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
plot(real(bSum),imag(bSum),'k-','LineWidth',2)

    % Product
%scatter(real(iPrd1)',imag(iPrd1)',2,dotColor(2,:),'o','filled')%,'MarkerFaceAlpha',0.5)
for idx = 1:size(bPrd1,1)
    p1 = scatter(real(iPrd2(idx,:)),imag(iPrd2(idx,:)),dotSize(1),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    p2 = scatter(real(IPrd(:,idx)),imag(IPrd(:,idx)),dotSize(2),dotColor(1,:),'o','filled');%,'MarkerFaceAlpha',1);
    lh = plot(real(bPrd1(idx,:)),imag(bPrd1(idx,:)),'-','LineWidth',2);
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
plot(real(bPrd),imag(bPrd),'k-','LineWidth',2)
fimplicit(@(x,y) (x-real(prod(offs))).^2 + ...
                 (y-imag(prod(offs))).^2 -...
                 (abs(offs(1))*r(2) + abs(offs(2))*r(1) + prod(r)).^2,...
                 'k--')
             
    % Inverse
scatter(real(iInv),imag(iInv),dotSize(1),dotColor(2,:),'o','filled')
scatter(real(IInv),imag(IInv),dotSize(3),dotColor(2,:),'o')
plot(real(bInv),imag(bInv),'k-','LineWidth',2)
plot(real(pInv),imag(pInv),'k:','LineWidth',2)
fimplicit(@(x,y) (x - real(oInv)).^2 + ...
                 (y - imag(oInv)).^2 -...
                 rInv.^2,...
                 'k--','LineWidth',2)

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
text(0.65,1.08,'A','FontSize',fontSize)
text(1.23,0.73,'B','FontSize',fontSize)
text(2.60,2.56,'A+B','FontSize',fontSize)
text(0.80,1.40,'A×B','FontSize',fontSize)
text(0.08,-0.19,'1/B','FontSize',fontSize)