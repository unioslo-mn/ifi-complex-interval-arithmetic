function r = quickSum(obj)
    
    % This function assumes that the input is a vertical array of convex 
    % polyarcular intervals, other input results unexpected behaviour

    % Prepare arcs
    N = length(obj);
    arcs = cell(N,1);
    for n = 1:N
        arcs{n} = [obj(n).Arcs{:} ; obj(n).Vertices{:}];
        [~,idx] = sort(arcs{n}.ArcAngle.Infimum);
        arcs{n} = arcs{n}(idx);
        if arcs{n}(1).ArcAngle.Infimum ~= -pi
            arcs{n} = [arcs{n}(end) ; arcs{n}];
            arcs{n}(1).ArcAngle.Infimum = -pi;
            arcs{n}(1).ArcAngle.Supremum = arcs{n}(1).ArcAngle.Supremum - 2*pi;
            arcs{n}(end).ArcAngle.Supremum= pi;
        end
    end
    
    % Sum arcs
    arcsSum = arcs{1};
    for n = 2:N
        arcsSum = quickPlus(arcsSum,arcs{n});
    end
    
    % Generate output interval
    r = ciat.PolyarcularInterval(arcsSum);

end

%% Function for quick convex addition

function r = quickPlus(arc1,arc2)

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