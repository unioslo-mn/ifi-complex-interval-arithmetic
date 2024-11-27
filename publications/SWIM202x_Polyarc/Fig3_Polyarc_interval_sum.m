
clear

% Interval parameters
rad = [[1.6 ; 1.8] , [1.7 ; 1.9]];
ang = [[0.2 ; 0.4] , [-0.2 ; -0.1] ];

% Arc sampling parameter
sCnt = 10;
sAxis = linspace(0,1,sCnt)';

% Generate arcs
pla=cell(2,1);
for m = 1:2
for n = 1:2
    pla{m}(n).typ = 'arc';
    pla{m}(n).cnt = 0;
    pla{m}(n).rad = rad(mod(n,2)+1,m) * (2*(n==1)-1);
    pla{m}(n).nor = ang(:,m) + (n==2);
    pla{m}(n).len = abs(pla{m}(n).rad*pi*diff(pla{m}(n).nor));
    pla{m}(n).smp = pla{m}(n).rad ...
                    * exp(1j*pi*(sAxis*diff(pla{m}(n).nor) ...
                                  + pla{m}(n).nor(1)));
    if(n==2)
        pla{m}(n).smp = flip(pla{m}(n).smp);
        pla{m}(n).nor = flip(pla{m}(n).nor);
    end
end
pla{m} = [pla{m} pla{m}(1)];
end

% Generate segments
seg=cell(2,1);
for m = 1:2
for n=1:2
    % Arc
    seg{m}(4*(n-1)+2) = pla{m}(n);
    % Points
    for p=1:2
        seg{m}(4*(n-1)+2*(p-1)+1).typ = 'point';
        seg{m}(4*(n-1)+2*(p-1)+1).cnt = pla{m}(n).smp(1+(p==2)*(sCnt-1));
        seg{m}(4*(n-1)+2*(p-1)+1).smp = seg{m}(4*(n-1)+2*(p-1)+1).cnt ...
                                     * ones(sCnt,1);
        seg{m}(4*(n-1)+2*(p-1)+1).nor = pla{m}(n).nor(p) ...
                                      + sort([0; p-1.5]);
        seg{m}(4*(n-1)+2*(p-1)+1).len = 0;
        seg{m}(4*(n-1)+2*(p-1)+1).rad = 0;
    end
    % Edge
    seg{m}(4*(n-1)+4).typ = 'edge';
    seg{m}(4*(n-1)+4).cnt = [pla{m}(n).smp(end);pla{m}(n+1).smp(1)];
    seg{m}(4*(n-1)+4).smp = sAxis*diff(seg{m}(4*(n-1)+4).cnt) + ...
                            seg{m}(4*(n-1)+4).cnt(1);
    seg{m}(4*(n-1)+4).nor = seg{m}(4*(n-1)+3).nor(2) * [1;1];
    seg{m}(4*(n-1)+4).len = abs(diff(seg{m}(4*(n-1)+4).cnt));
    seg{m}(4*(n-1)+4).rad = 0;
end
end

% Split segments with pi crossing normal interval
for m = 1:2
    n = 1;
    while n<=length(seg{m})
        if (seg{m}(n).nor(1)<1 && seg{m}(n).nor(2)>1) || ...
           (seg{m}(n).nor(1)>1 && seg{m}(n).nor(2)<1)
            seg{m} = [seg{m}(1:n),...
                          seg{m}(n),...
                          seg{m}(n+1:end)];
            norTemp = seg{m}(n).nor;
            seg{m}(n).nor = [norTemp(2);1];
            seg{m}(n+1).nor = [1;norTemp(1)];
        end
        n=n+1;
    end
end

% Wrap and reorder segments
for m = 1:2
    for n = 1:length(seg{m})
        seg{m}(n).nor = wrapToPi(seg{m}(n).nor*pi)/pi;
        if abs(seg{m}(n).nor(1)) == 1
            seg{m}(n).nor(1) = sign(seg{m}(n).nor(2));
        elseif abs(seg{m}(n).nor(2)) == 1
            seg{m}(n).nor(2) = sign(seg{m}(n).nor(1));
        end
        if seg{m}(n).rad >= 0
            seg{m}(n).nor = sort(seg{m}(n).nor,'ascend');
        else
            seg{m}(n).nor = sort(seg{m}(n).nor,'descend');
        end
    end
    [~,segOrder]=sortrows([seg{m}.nor]',[1,2]);
    seg{m} = seg{m}(segOrder);
    segIdx = find([seg{m}.rad]<0);
    for k = 1:length(segIdx)
        if segIdx <= 2
            seg{m} = [seg{m}(segIdx(k)),...
                      seg{m}(segIdx(k)-1),...
                      seg{m}(segIdx(k)+1:length(seg{m}))];
        elseif segIdx < length(seg{m})
            seg{m} = [seg{m}(1:segIdx(k)-2),...
                      seg{m}(segIdx(k)),...
                      seg{m}(segIdx(k)-1),...
                      seg{m}(segIdx(k)+1:length(seg{m}))];
        else
            seg{m} = [seg{m}(1:segIdx(k)-2),...
                      seg{m}(segIdx(k)),...
                      seg{m}(segIdx(k)-1)];
        end
    end
end

% Generate Gauss map matching matrix
gauMatch = zeros(8);
gauAngle = zeros(16,4);
segCnt = 1;
for k = 1:length(seg{1})
for l = 1:length(seg{2})
    % Get segments
    seg1 = seg{1}(k);
    seg2 = seg{2}(l);

    % Check if there is a Gauss map match
        % Get normal intervals
    nor1 = wrapToPi(seg1.nor);
    nor2 = wrapToPi(seg2.nor);
       % Calculate intersection
    nor3 = [max([min(nor1),min(nor2)]) ; ...
            min([max(nor1),max(nor2)]) ];
        % Check if intersection is empty
    gauMatch(k,l) = nor3(1) <= nor3(2);

    % Add segments if there is a match
    if gauMatch(k,l) > 0
        if strcmp(seg1.typ,'point')
            seg{3}(segCnt) = seg2;
            seg{3}(segCnt).cnt = seg2.cnt + seg1.cnt;
            seg{3}(segCnt).smp = seg2.smp + seg1.cnt;
            seg{3}(segCnt).nor = nor3;
            segCnt = segCnt+1;
        elseif strcmp(seg2.typ,'point')
            seg{3}(segCnt) = seg1;
            seg{3}(segCnt).cnt = seg1.cnt + seg2.cnt;
            seg{3}(segCnt).smp = seg1.smp + seg2.cnt;
            seg{3}(segCnt).nor = nor3;
            segCnt = segCnt+1;
        elseif strcmp(seg1.typ,'edge')
            seg{3}(segCnt) = seg2;
            seg{3}(segCnt).cnt = seg2.cnt + seg1.cnt(1);
            seg{3}(segCnt).smp = seg2.smp + seg1.cnt(1);
            seg{3}(segCnt).nor = nor3;
            %
            seg{3}(segCnt+1) = seg2;
            seg{3}(segCnt+1).cnt = seg2.cnt + seg1.cnt(2);
            seg{3}(segCnt+1).smp = seg2.smp + seg1.cnt(2);
            seg{3}(segCnt+1).nor = nor3;
            segCnt = segCnt+2;
        elseif strcmp(seg2.typ,'edge')
            seg{3}(segCnt) = seg1;
            seg{3}(segCnt).cnt = seg1.cnt + seg2.cnt(1);
            seg{3}(segCnt).smp = seg1.smp + seg2.cnt(1);
            seg{3}(segCnt).nor = nor3;
            %
            seg{3}(segCnt+1) = seg1;
            seg{3}(segCnt+1).cnt = seg1.cnt + seg2.cnt(2);
            seg{3}(segCnt+1).smp = seg1.smp + seg2.cnt(2);
            seg{3}(segCnt+1).nor = nor3;
            segCnt = segCnt+2;
        else
            seg{3}(segCnt) = seg1;
            seg{3}(segCnt).cnt = seg1.cnt + seg2.cnt;
            seg{3}(segCnt).rad = seg1.rad + seg2.rad;
            seg{3}(segCnt).nor = nor3;
            seg{3}(segCnt).smp = seg{3}(segCnt).rad ...
                                   * exp(1j*pi*(sAxis*diff(seg{3}(segCnt).nor) ...
                                  + seg{3}(segCnt).nor(1)));
            segCnt = segCnt+1;
        end
    end
end
end

% Delete arcs with zero normal range
n = 1;
while n <= length(seg{3})
    if strcmp(seg{3}(n).typ,'arc') && diff(seg{3}(n).nor)<1e-6
        if n==1
            seg{3} = seg{3}(n+1:length(seg{3}));
        elseif n < length(seg{3})
            seg{3} = [seg{3}(1:n-1),...
                      seg{3}(n+1:length(seg{3}))];
        else
            seg{3} = seg{3}(1:n-1);
        end
    else
        n = n+1;
    end
end

% Wrap and reorder segments
for m = 3
    for n = 1:length(seg{m})
        seg{m}(n).nor = wrapToPi(seg{m}(n).nor*pi)/pi;
        if abs(seg{m}(n).nor(1)) == 1
            seg{m}(n).nor(1) = sign(seg{m}(n).nor(2));
        elseif abs(seg{m}(n).nor(2)) == 1
            seg{m}(n).nor(2) = sign(seg{m}(n).nor(1));
        end
        if seg{m}(n).rad >= 0
            seg{m}(n).nor = sort(seg{m}(n).nor,'ascend');
        else
            seg{m}(n).nor = sort(seg{m}(n).nor,'descend');
        end
    end
    [~,segOrder]=sortrows([seg{m}.nor]',[1,2]);
    seg{m} = seg{m}(segOrder);
    segIdx = find([seg{m}.rad]<0);
    for k = 1:length(segIdx)
        if segIdx <= 2
            seg{m} = [seg{m}(segIdx(k)),...
                      seg{m}(segIdx(k)-1),...
                      seg{m}(segIdx(k)+1:length(seg{m}))];
        elseif segIdx < length(seg{m})
            seg{m} = [seg{m}(1:segIdx(k)-2),...
                      seg{m}(segIdx(k)),...
                      seg{m}(segIdx(k)-1),...
                      seg{m}(segIdx(k)+1:length(seg{m}))];
        else
            seg{m} = [seg{m}(1:segIdx(k)-2),...
                      seg{m}(segIdx(k)),...
                      seg{m}(segIdx(k)-1)];
        end
    end
end

% Add interval boundaries
bnd1 = [seg{1}.smp];
bnd2 = [seg{2}.smp];
bnd3 = bnd1(:) + bnd2(:).';

% Plot parameters
% col = lines(100);
% col = hsv(100);
col = turbo(100);
gauRad = [0.2 , 0.5 , 0.8];
gauSca = 0.2 * [1,1,1];
norLen = 0.2;

%% Plot
figure(1);clf;hold on
% subplot(1,2,1);hold on;title('Complex interval sum')
fimplicit(@(x,y) x.^2+y.^2-1,'k-','LineWidth',1)
scatter(0,0,10,'+')
scatter(real(bnd3),imag(bnd3),1,[1 1 1]*0.8,'filled')
axis equal
% subplot(1,2,2);hold on;title('Gauss map')
fimplicit(@(x,y) x.^2+y.^2-((gauRad(1)+gauSca(1)+gauRad(2))/2)^2,'k-','LineWidth',1)
fimplicit(@(x,y) x.^2+y.^2-(gauRad(2)+gauSca(2)*1.1)^2,'k-','LineWidth',1)
%1fimplicit(@(x,y) x.^2+y.^2-(gauRad(3)+gauSca(3)*1.1)^2,'k-','LineWidth',1)
axis equal
for m = 1:3
    gauSca(m) = gauSca(m) / sum([seg{m}.len]);
for n = 1:length(seg{m})
    % Draw interval
    % subplot(1,2,1);
    c = round(round(mean(seg{m}(n).nor)/2+0.5,2)*100);
    if strcmp(seg{m}(n).typ,'point')
        scatter(real(seg{m}(n).cnt),imag(seg{m}(n).cnt),50,col(c,:),...
                'o','filled','MarkerFaceColor',col(c,:))
    else
        plot(real(seg{m}(n).smp),imag(seg{m}(n).smp),...
                '-','color',col(c,:),'LineWidth',2)
    end
    % Draw Gauss map
    gau = exp(1j*pi*(sAxis*diff(seg{m}(n).nor)+seg{m}(n).nor(1))) .* ...
                (sAxis * seg{m}(n).len * gauSca(m) + gauRad(m));
    plot(real(gau),imag(gau),...
        '-','color',col(c,:),'LineWidth',2) 
    norAng = linspace(seg{m}(n).nor(1),seg{m}(n).nor(2),sCnt)';
    quiver(real(seg{m}(n).smp),imag(seg{m}(n).smp),...
           norLen*cos(norAng*pi),norLen*sin(norAng*pi),...
           '-','color',col(c,:),'AutoScale','off');

    gauRad(m) = gauRad(m) + seg{m}(n).len * gauSca(m);
end
end

% Settings
xlabel('Real')
ylabel('Imag')
text(1,1.35,'A','HorizontalAlignment','center')
text(1.6,-.8,'B','HorizontalAlignment','center')
text(2.5,.6,'A+B','HorizontalAlignment','center')
text(0,0.1,'A','HorizontalAlignment','center')
text(0,0.5,'B','HorizontalAlignment','center')
text(0,.8,'A+B','HorizontalAlignment','center')
