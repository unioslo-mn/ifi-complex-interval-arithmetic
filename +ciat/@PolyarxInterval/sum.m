function r = sum(obj)
    
    % This function assumes that the input is a vertical array of convex 
    % polyarcular intervals, other input results unexpected behaviour

    arx = obj.Arx;
    N = length(arx);

    % Sum arcs
    arxSum = arx{1};
    for n = 2:N
        arxSum = quickPlus(arxSum,arx{n});
    end
    
    % Generate output interval
    r = ciat.PolyarxInterval(arxSum);

end

%% Function for quick convex addition

function r = quickPlus(arx1,arx2)

    % Extract parameters
    cen1 = complex(arx1(:,1),arx1(:,2));
    rad1 = arx1(:,3);
    ang1 = arx1(:,4);
    cen2 = complex(arx2(:,1),arx2(:,2));
    rad2 = arx2(:,3);
    ang2 = arx2(:,4);
    
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

        if n3 > 1 && ang3(n3) <= ang3(n3-1) + eps10
            n3 = n3 - 1;
        end
    end
    cen3 = cen3(1:n3,:);
    rad3 = rad3(1:n3,:);
    ang3 = ang3(1:n3,:);
       
    % Generate polyarc
    r = [real(cen3),imag(cen3),rad3,ang3];
end