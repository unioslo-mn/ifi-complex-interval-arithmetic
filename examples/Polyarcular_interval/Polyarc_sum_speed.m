clear
% close all

%% Generate two random polar intervals

% Define random polar intervals
N = 10;
absMin = rand(N,1);
absMax = absMin + rand(N,1)/10;
angMin = 2*pi*rand(N,1);
angMax = angMin + 2*pi*rand(N,1)/20;
pI = ciat.PolarInterval(absMin,absMax,angMin,angMax);

% Convert them to various other formats
pcI = ciat.RectangularInterval(pI);
pgI = ciat.PolygonalInterval(pI);
paI = ciat.PolyarcularInterval(pI);
pxI = paI.convexify;

% Sum rectangular and polyarcular as reference
tic
pcIsum = sum(pcI);
recTime = toc;
tic
paIsum = sum(paI);
arcTime = toc;

%% Sum polygonal intervals

% Prepare objects
obj1 = pgI(1);
obj2 = pgI(2);

% Measure time
tic
pgIsum = pgI(1);
for n = 2:N
    pgIsum = gonPlus(pgIsum,pgI(n));
end
gonTime = toc;

%% Sum polyarcular intervals

% Prepare arcs
arcs = cell(N,1);
for n = 1:N
    arcs{n} = [pxI(n).Arcs{:} ; pxI(n).Vertices{:}];
    [~,idx] = sort(arcs{n}.ArcAngle.Infimum);
    arcs{n} = arcs{n}(idx);
    if arcs{n}(1).ArcAngle.Infimum ~= -pi
        arcs{n} = [arcs{n}(end) ; arcs{n}];
        arcs{n}(1).ArcAngle.Infimum = -pi;
        arcs{n}(1).ArcAngle.Supremum = arcs{n}(1).ArcAngle.Supremum - 2*pi;
        arcs{n}(end).ArcAngle.Supremum= pi;
    end
end

% Measure time
tic
arcsSum = arcs{1};
for n = 2:N
    arcsSum = arxPlus(arcsSum,arcs{n});
end
arxTime = toc;

pxIsum = ciat.PolyarcularInterval(arcsSum);


%% Report

% Measure areas
recArea = pcIsum.Area;
gonArea = pgIsum.Area;
arcArea = paIsum.Area;
arxArea = pxIsum.Area;

% Print report
sprintf(['Rectangular interval: Area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polygonal interval area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n'...
         'Polyarcular (convex) area: %0.4f (tightness: %0.1f%%), Time: %0.1fms\n',...
         'Polyarcular (concave) area: %0.4f (tightness: %0.1f%%), Time: %0.1fms'], ...
         recArea, arcArea / recArea * 100, recTime*1e3,...
         gonArea, arcArea / gonArea * 100, gonTime*1e3,...
         arxArea, arcArea / arxArea * 100, arxTime*1e3,...
         arcArea, arcArea / arcArea * 100, arcTime*1e3)

%% Plot
% figure;clf
cla;hold on;axis equal
paI.plot('b');
paI.plotGaussMap(0.05,'c');
paIsum.plot('k','linewidth',2);
pgIsum.plot('c','linewidth',2);
pxIsum.plot('r','linewidth',2);



%% Stripped down polygonal sum algorithm

function r = gonPlus(obj1,obj2)

    v = obj1.Points{:};
    w = obj2.Points{:};
    
    % Handle exception when one of the inputs is a degenerate interval
    if length(v)==1 || length(w)==1 
        points = reshape(v + w.' ,[],1);
        r = ciat.PolygonalInterval(points);
        return
    end
    
    i = 1; j = 1; 
    eps10 = 10*eps; % needed so avoid skipping vertices due to numerical precision
    
    I = length(obj1.Points{:}) + 1; 
    v = [v; v(1); v(2)];
    v_arg = ciat.wrapTo2Pi( angle( v(2:end) - v(1:end-1) ));
    v_arg(end) = v_arg(end) + 2*pi; % otherwise it wraps around
    
    J = length(obj2.Points{:}) + 1; 
    w = [w; w(1); w(2)];
    w_arg = ciat.wrapTo2Pi( angle( w(2:end) - w(1:end-1) ));
    w_arg(end) = w_arg(end) + 2*pi; % otherwise it wraps around
    
    p = zeros( I + J, 1);
    n = 0;
    
    % v_arg(i) = angle(v_i+1 - v_i), which is why we also
    % repeat when i=I and j=J (as opposed to in DeBerg2008). 
    % The loop breaks in the following turn.
    while (i <= I) && (j <= J) % continue finding more points
        n = n+1;
        p(n) = v(i) + w(j); % add what we found
    
        if     v_arg(i) < w_arg(j) + eps10 
            i = i + 1;
        elseif v_arg(i) > w_arg(j) + eps10
            j = j + 1;
        else
            i = i + 1;
            j = j + 1;
        end            
    end
    
    r = ciat.PolygonalInterval(p(1:n-1));

end

%% Stripped down convex polyarcular plus algorithm

function r = arxPlus(arc1,arc2)

    % Extract parameters
    cen1 = arc1.Center;
    rad1 = arc1.Radius;
    ang1 = arc1.ArcAngle.Supremum;
    cen2 = arc2.Center;
    rad2 = arc2.Radius;
    ang2 = arc2.ArcAngle.Supremum;
    
    % Taken from the polygonal plus function
    N1 = size(ang1,1);
    N2 = size(ang2,1);
    N3 = N1 + N2;
    cen3 = zeros(N3,1);
    rad3 = zeros(N3,1);
    ang3 = zeros(N3,1);
    n1 = 1;
    n2 = 1;
    n3 = 0;
    eps10 = eps*10;
    while (n1 <= N1) && (n2 <= N2) % continue finding more points
        n3 = n3 + 1;
    
        % Sum arcs
        cen3(n3) = cen1(n1) + cen2(n2);
        rad3(n3) = rad1(n1) + rad2(n2);
        ang3(n3) = min(ang1(n1),ang2(n2));
        
        % Increment index
        if  ang1(n1) < ang2(n2) + eps10 
            n1 = n1 + 1;
        elseif ang1(n1) > ang2(n2) + eps10
            n2 = n2 + 1;
        else
            n1 = n1 + 1;
            n2 = n2 + 1;
        end            
    end
    cen3 = cen3(1:n3,:);
    rad3 = rad3(1:n3,:);
    ang3 = ang3(1:n3,:);
    
    
    % Calculate angle infimum and supremum values
    angInf = [-pi ; ang3(1:end-1)];
    angSup = ang3;
    
    % Generate polyarc
    r = ciat.Arc(cen3,rad3, ciat.RealInterval(angInf,angSup));
end
