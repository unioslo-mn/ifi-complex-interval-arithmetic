function r = plusConvex(obj1,obj2)

    % Extract arcs including vertices
    arc1 = [obj1.Arcs{:} ; obj1.Vertices{:}];
    arc2 = [obj2.Arcs{:} ; obj2.Vertices{:}];

    % Split arcs with angle interval including pi
    arc1 = splitArc(arc1);
    arc2 = splitArc(arc2);

    % Sort arcs by angle
    [~,idx] = sort(arc1.ArcAngle.Infimum);
    arc1 = arc1(idx);
    [~,idx] = sort(arc2.ArcAngle.Infimum);
    arc2 = arc2(idx);

    % Taken from the polygonal plus function
    ang1 = arc1.ArcAngle.Supremum;
    ang2 = arc2.ArcAngle.Supremum;
    N1 = length(arc1);
    N2 = length(arc2);
    N3 = N1 + N2;
    arc3(N3,1) = ciat.Arc;
    n1 = 1;
    n2 = 1;
    n3 = 0;
    eps10 = eps*10;
    while (n1 <= N1) && (n2 <= N2) % continue finding more points
        n3 = n3 + 1;
        % arc3(n3) = arc1(n1) + arc2(n2); % add what we found
        arc3(n3) = ciat.Arc( arc1(n1).Center + arc2(n2).Center ,...
                             arc1(n1).Radius + arc2(n2).Radius ,...
                             cap( arc1(n1).ArcAngle, arc2(n2).ArcAngle));

        if  ang1(n1) < ang2(n2) + eps10 
            n1 = n1 + 1;
        elseif ang1(n1) > ang2(n2) + eps10
            n2 = n2 + 1;
        else
            n1 = n1 + 1;
            n2 = n2 + 1;
        end            
    end
    arc3 = arc3(1:n3);

    % If possible, merge elements at the pi crossing
    if arc3(1).Center==arc3(end).Center && arc3(1).Radius==arc3(end).Radius
        arc3(1).ArcAngle = ciat.RealInterval(arc3(end).ArcAngle.Infimum,...
                                      arc3(1).ArcAngle.Supremum + 2*pi);
        arc3 = arc3(1:n3-1);
    end

    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3);
    
end

%% Function for splitting arcs with angle interval crossing pi
% (this is different from the ciat.Arc splitAngle function)

function arcOut = splitArc(arcIn)
    % Find arcIn with angle interval including -pi or pi
    mask = (arcIn.ArcAngle.isin(-pi) & ...
            arcIn.ArcAngle.Infimum ~= -pi & ...
            arcIn.ArcAngle.Supremum ~= -pi) | ... 
           (arcIn.ArcAngle.isin(pi) & ...
            arcIn.ArcAngle.Infimum ~= pi & ...
            arcIn.ArcAngle.Supremum ~= pi);

    if any(mask,'all')
        % Create new arcIn with angle from pi to the supremum 
        newArc = arcIn(mask);
        newArc.ArcAngle = ciat.RealInterval(-pi * ones(length(newArc),1),...
                                        wrapToPi(newArc.ArcAngle.Supremum));
        
        % Modify the original arcIn to have angle from the infimum to pi
        arcIn(mask).ArcAngle = ciat.RealInterval(...
                                    wrapToPi(arcIn(mask).ArcAngle.Infimum),...
                                    pi * ones(sum(mask),1) );
        
        % Wrap the angles of all the other arcs
        arcIn(~mask).ArcAngle.Infimum = wrapToPi(arcIn(~mask).ArcAngle.Infimum);
        arcIn(~mask).ArcAngle.Supremum = wrapToPi(arcIn(~mask).ArcAngle.Supremum);
    
        % Add new arcIn to the array
        arcOut = [arcIn ; newArc];
    else
        arcOut = arcIn;
    end
end