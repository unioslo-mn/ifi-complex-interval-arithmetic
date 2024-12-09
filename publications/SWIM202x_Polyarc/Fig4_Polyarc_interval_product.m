
clear

% Interval parameters
% rad = [[0.3 ; 1.9] , [0.6 ; 2.8]];
% ang = [ [-0.7 ; -0.3]+0.35 , [0.6 ; 1.1]-0.3 ];

% Interval parameters
rad = [[0.3 ; 1.9] , [0.6 ; 2.8]];
ang = [[-0.6 ; -0.3] , [-0.9 ; -0.4] ];

% Arc sampling parameter
sCnt = 20;
sAxis = linspace(0,1,sCnt)';


% Generate polar intervals
A_p = ciat.PolarInterval(rad(1,1), rad(2,1), ang(1,1)*pi, ang(2,1)*pi);
B_p = ciat.PolarInterval(rad(1,2), rad(2,2), ang(1,2)*pi, ang(2,2)*pi);

% Cast to polyarc and trasnlate to real 1 center
A_a = (ciat.PolyarcularInterval(A_p)-1);
B_a = (ciat.PolyarcularInterval(B_p)-1);

% Multiply polyarcular intervals
C_a = A_a * B_a; 

% Multiply interval boundaries
A_smp = A_a.sample(sCnt);
B_smp = B_a.sample(sCnt);
C_smp = A_smp(:) * B_smp(:).';

% Multiply arcs and edges
ArcTimesArc     = A_a.Arcs .* B_a.Arcs.';           
ArcTimesEdge    = A_a.Arcs .* B_a.Edges.';         
EdgeTimesEdge   = A_a.Edges .* B_a.Edges.';       
VertTimesArc    = A_a.Vertices .* B_a.Arcs.';      
VertTimesEdge   = A_a.Vertices .* B_a.Edges.';   
ArcTimesVert    = A_a.Arcs .* B_a.Vertices.';      
EdgeTimesVert   = A_a.Edges .* B_a.Vertices.';
VertTimesVert   = A_a.Vertices .* B_a.Vertices.';
ArcTimesArc     = ArcTimesArc  (~isnan(ArcTimesArc   ));
ArcTimesEdge    = ArcTimesEdge (~isnan(ArcTimesEdge  ));
EdgeTimesEdge   = EdgeTimesEdge(~isnan(EdgeTimesEdge ));
VertTimesArc    = VertTimesArc (~isnan(VertTimesArc  ));
VertTimesEdge   = VertTimesEdge(~isnan(VertTimesEdge ));
ArcTimesVert    = ArcTimesVert (~isnan(ArcTimesVert  ));
EdgeTimesVert   = EdgeTimesVert(~isnan(EdgeTimesVert ));
VertTimesVert   = VertTimesVert(~isnan(VertTimesVert ));

% % Select segments that are not NaN
C_arcs = [ArcTimesArc;ArcTimesEdge;EdgeTimesEdge;...
          VertTimesArc;ArcTimesVert];
C_edges = [VertTimesEdge;EdgeTimesVert];
C_vertices = VertTimesVert;




%% Plot

% Plot parameters
col = lines(100);
gauRad = [0.2 , 0.5 , 0.8];
gauSca = 0.2 * [1,1,1];
norLen = 0.5;
gauOffset = 0;
lineWidth = 2;

% Initialize figure and subplots
figure(1);clf;
    % Log-Gauss map diagram subplot
subplot(2,3,1);hold on;axis equal
scatter(real(gauOffset),imag(gauOffset),10,'+')
fimplicit(@(x,y) (x-real(gauOffset)).^2+(y-imag(gauOffset)).^2 - ...
                        ((gauRad(1)+gauSca(1)+gauRad(2))/2)^2,'k:','linewidth',2)
fimplicit(@(x,y) (x-real(gauOffset)).^2+(y-imag(gauOffset)).^2 - ...
                        (gauRad(2)+gauSca(2)*1.1)^2,'k:','linewidth',2)
fimplicit(@(x,y) (x-real(gauOffset)).^2+(y-imag(gauOffset)).^2 - ...
                        (gauRad(3)+gauSca(3)*1.1)^2,'k:','linewidth',2)
text(0,0.2,'A','HorizontalAlignment','center')
text(0,0.55,'B','HorizontalAlignment','center')
text(0,0.9,'A \times B','HorizontalAlignment','center')
xlabel('Real')
ylabel('Imag')
fontsize(30,'points')
annotation('textbox',[0.07,.83,.1,.1],'String','\textbf{a)}', ...
           'FontSize',30,...
           'BackgroundColor','none',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

    % Complex plane subplot
subplot(2,3,[2,3,5,6]);hold on;axis equal
scatter(real(C_smp),imag(C_smp),1,'g.');
xlabel('Real')
ylabel('Imag')
text(-1.07,-1.2,'A','HorizontalAlignment','center')
text(-2.63,-1.5,'B','HorizontalAlignment','center')
text(-0.28,5,'A \times B','HorizontalAlignment','center')
fontsize(30,'points')
annotation('textbox',[0.35,.83,.1,.1],'String','\textbf{b)}', ...
           'FontSize',30,...
           'BackgroundColor','none',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');

    % Log-domain subplot
subplot(2,3,4);hold on;axis equal
text(4.5,-0.2,'A','HorizontalAlignment','center')
text(3.6,1,'B','HorizontalAlignment','center')
text(1.5,1,'A \times B','HorizontalAlignment','center')
xlabel('Arg')
ylabel('Log-Abs')
fontsize(30,'points')
annotation('textbox',[0.07,.35,.1,.1],'String','\textbf{c)}', ...
           'FontSize',30,...
           'BackgroundColor','none',...
           'EdgeColor','none',...
           'HorizontalAlignment','right', ...
           'Interpreter','latex');


% Plot segments
allSegments = {A_a.Boundary ; B_a.Boundary ; C_a.Boundary};
for m = 1:length(allSegments)
    % Calculate total boundary length of the interval
    bndLen = 0;
    for n = 1:length(allSegments{m})
        bndLen = bndLen + allSegments{m}{n}.Length;
    end
    gauSca(m) = gauSca(m) / bndLen;

    % Plot segments
    for n = 1:length(allSegments{m})
        seg = allSegments{m}{n};
        c = max(round(round(wrapToPi(seg.LogGaussMap.mid)/pi/2+0.5,2)*100),1);
        subplot(2,3,1);
        loggau = exp(1j * seg.LogGaussMap.sample(sCnt)).' .* ...
                    (sAxis * seg.Length * gauSca(m) + gauRad(m)) + ...
                    gauOffset;
        plot(real(loggau),imag(loggau),'-','color',col(c,:),'LineWidth',2) 
        gauRad(m) = gauRad(m) + seg.Length * gauSca(m);

        % Plot in the complex plane
        subplot(2,3,[2,3,5,6]);
        seg.plot('color',col(c,:),'linewidth',lineWidth);
        seg.plotLogGaussMap(norLen,'color',col(c,:),'arrowCount',10);

        % Plot in log-domain
        subplot(2,3,4);
        if isa(seg,'ciat.Arc') && seg.Radius == 0 
            isPoint = 1;
        else
            isPoint = 0;
        end

        if ~isPoint
            logSmp = log(seg.sample(sCnt));
            plot(wrapTo2Pi(imag(logSmp)),real(logSmp),...
                    '-','color',col(c,:),'LineWidth',2)
        else
            logPnt = log(seg.Center);
            scatter(wrapTo2Pi(imag(logPnt)),real(logPnt),50,col(n,:),...
                 'o','filled','MarkerFaceColor',col(n,:))
        end
    end
end

