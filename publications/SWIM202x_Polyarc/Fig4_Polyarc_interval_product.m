
clear

% Interval parameters
rad = [[0.3 ; 1.9] , [0.6 ; 2.8]];
ang = [ [-0.7 ; -0.3] , [0.6 ; 1.1] ];

% Arc sampling parameter
sCnt = 50;
sAxis = linspace(0,1,sCnt)';

% Plot parameters
col = lines(100);
gauRad = [0.2 , 0.5 , 0.8]*3;
gauSca = 0.8 * [1,1,1];
gauOffset = -5+6i;

% Generate polar intervals
A_p = ciat.PolarInterval(rad(1,1), rad(2,1), ang(1,1)*pi, ang(2,1)*pi);
B_p = ciat.PolarInterval(rad(1,2), rad(2,2), ang(1,2)*pi, ang(2,2)*pi);

% Cast to polyarc and tranlate to real 1 center
A_a = (ciat.PolyarcularInterval(A_p)-1);
B_a = (ciat.PolyarcularInterval(B_p)-1);


% Add interval boundaries
A_smp = A_a.sample(sCnt);
B_smp = B_a.sample(sCnt);
C_smp = A_smp(:) * B_smp(:).';
C_idx = boundary(real(C_smp(:)),imag(C_smp(:)),0.7);
C_bnd = C_smp(C_idx);

% Log-Gauss map intersections
% LGM_arc_arc = ~isnan(cap(A_a.Arcs.LogGaussMap,B_a.Arcs.LogGaussMap.'))
% LGM_arc_edge = ~isnan(cap([A_a.Arcs.LogGaussMap; A_a.Edges.LogGaussMap],...
%                           [B_a.Arcs.LogGaussMap; B_a.Edges.LogGaussMap].'))
% LGM_edge_edge = ~isnan(cap(A_a.Edges.LogGaussMap,B_a.Edges.LogGaussMap.'))

%% Plot
figure(2);clf;hold on;axis equal
A_a.plot('b');
B_a.plot('r');
% A_a.plotLogGaussMap(0.2,'b');
% B_a.plotLogGaussMap(0.2,'r');
scatter(real(C_smp),imag(C_smp),1,'g.');
plot(real(C_bnd),imag(C_bnd),'k-');

%%
% Generate arcs
plaOffs = -[1,1];
pla=cell(2,1);
for m = 1:2
for n = 1:2
    pla{m}(n).typ = 'arc';
    pla{m}(n).cnt = plaOffs(m);
    pla{m}(n).rad = rad(mod(n,2)+1,m) * (2*(n==1)-1);
    pla{m}(n).nor = ang(:,m) + (n==2);
    pla{m}(n).len = abs(pla{m}(n).rad*pi*diff(pla{m}(n).nor));
    pla{m}(n).smp = pla{m}(n).cnt + pla{m}(n).rad ...
                    * exp(1j*pi*(sAxis*diff(pla{m}(n).nor) ...
                                  + pla{m}(n).nor(1)));
    if(n==2)
        pla{m}(n).smp = flip(pla{m}(n).smp);
        pla{m}(n).nor = flip(pla{m}(n).nor);
    end
    pla{m}(n).log = [pla{m}(n).nor(1) - angle(pla{m}(n).smp(1))/pi;...
                    pla{m}(n).nor(2) - angle(pla{m}(n).smp(end))/pi  ];
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
        seg{m}(4*(n-1)+2*(p-1)+1).log = seg{m}(4*(n-1)+2*(p-1)+1).nor - ...
                                     angle(seg{m}(4*(n-1)+2*(p-1)+1).cnt)/pi;
    end
    % Edge
    seg{m}(4*(n-1)+4).typ = 'edge';
    seg{m}(4*(n-1)+4).cnt = [pla{m}(n).smp(end);pla{m}(n+1).smp(1)];
    seg{m}(4*(n-1)+4).smp = sAxis*diff(seg{m}(4*(n-1)+4).cnt) + ...
                            seg{m}(4*(n-1)+4).cnt(1);
    seg{m}(4*(n-1)+4).nor = seg{m}(4*(n-1)+3).nor(2) * [1;1];
    seg{m}(4*(n-1)+4).len = abs(diff(seg{m}(4*(n-1)+4).cnt));
    seg{m}(4*(n-1)+4).log = seg{m}(4*(n-1)+4).nor - ...
                         angle(seg{m}(4*(n-1)+4).smp([1,end]))/pi;
end
end

% Split segments with pi crossing log-normal interval
for m = 1:2
    n = 1;
    while n<=length(seg{m})
        if (seg{m}(n).log(1)<1 && seg{m}(n).log(2)>1) || ...
           (seg{m}(n).log(1)>1 && seg{m}(n).log(2)<1)
            seg{m} = [seg{m}(1:n),...
                          seg{m}(n),...
                          seg{m}(n+1:end)];
            logTemp = seg{m}(n).log;
            seg{m}(n).log = [logTemp(2);1];
            seg{m}(n+1).log = [1;logTemp(1)];
        end
        n=n+1;
    end
end

% Generate log-Gauss map match matrix
gauMatch = zeros(8);
gauAngle = zeros(16,4);
segCnt = 1;
for k = 1:length(seg{1})
for l = 1:length(seg{2})
    % Get segments
    seg1 = seg{1}(k);
    seg2 = seg{2}(l);

    % Check if there is a Gauss map match
        % Get logmal intervals
    log1 = wrapToPi(seg1.log);
    log2 = wrapToPi(seg2.log);
       % Calculate intersection
    log3 = [max([min(log1),min(log2)]) ; ...
            min([max(log1),max(log2)]) ];
        % Check if intersection is empty
    gauMatch(k,l) = log3(1) <= log3(2);

    % Add segments if there is a match
    if gauMatch(k,l) > 0
        if strcmp(seg1.typ,'point')
            seg{3}(segCnt) = seg2;
            seg{3}(segCnt).cnt = seg2.cnt * seg1.cnt;
            seg{3}(segCnt).smp = seg2.smp * seg1.cnt;
            seg{3}(segCnt).log = log3;
            segCnt = segCnt+1;
        elseif strcmp(seg2.typ,'point')
            seg{3}(segCnt) = seg1;
            seg{3}(segCnt).cnt = seg1.cnt * seg2.cnt;
            seg{3}(segCnt).smp = seg1.smp * seg2.cnt;
            seg{3}(segCnt).log = log3;
            segCnt = segCnt+1;
        elseif strcmp(seg1.typ,'edge')
            
        elseif strcmp(seg2.typ,'edge')
            
        else
            
        end
    end
end
end

% Add interval boundaries
bnd1 = [seg{1}.smp];
bnd2 = [seg{2}.smp];
bnd3 = bnd1(:) * bnd2(:).';



% figure(2);clf;hold on;axis equal
% fimplicit(@(x,y) x.^2+y.^2-1,'k:')
scatter(real(bnd3),imag(bnd3),'g.')

scatter(real(gauOffset),imag(gauOffset),10,'+')
fimplicit(@(x,y) (x-real(gauOffset)).^2+(y-imag(gauOffset)).^2 - ...
                        ((gauRad(1)+gauSca(1)+gauRad(2))/2)^2,'k:')
fimplicit(@(x,y) (x-real(gauOffset)).^2+(y-imag(gauOffset)).^2 - ...
                        (gauRad(2)+gauSca(2)*1.1)^2,'k:')
fimplicit(@(x,y) (x-real(gauOffset)).^2+(y-imag(gauOffset)).^2 - ...
                        (gauRad(3)+gauSca(3)*1.1)^2,'k:')
% axis equal
for m = 1:3
    gauSca(m) = gauSca(m) / sum([seg{m}.len]);
for n = 1:length(seg{m})
    % Draw interval
    if strcmp(seg{m}(n).typ,'point')
        scatter(real(seg{m}(n).cnt),imag(seg{m}(n).cnt),50,col(n,:),...
                'o','filled','MarkerFaceColor',col(n,:))
    else
        plot(real(seg{m}(n).smp),imag(seg{m}(n).smp),...
                '-','color',col(n,:),'LineWidth',2)
    end
    % Draw Gauss map
    gau = exp(1j*pi*(sAxis*diff(seg{m}(n).nor)+seg{m}(n).nor(1))) .* ...
                (sAxis * seg{m}(n).len * gauSca(m) + gauRad(m));
    loggau = exp(1j*pi*(sAxis*diff(seg{m}(n).log)+seg{m}(n).log(1))) .* ...
                    (sAxis * seg{m}(n).len * gauSca(m) + gauRad(m)) + ...
                    gauOffset;
    plot(real(loggau),imag(loggau),...
        '-','color',col(n,:),'LineWidth',2) 

    gauRad(m) = gauRad(m) + seg{m}(n).len * gauSca(m);
end
end
