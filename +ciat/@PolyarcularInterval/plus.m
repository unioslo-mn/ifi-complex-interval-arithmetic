function r = plus(obj1,obj2)

% Addition of polyarcular intervals (+ operator)
%
% This function creates the polyarcular intervals representing the sum
% of two sets of polyarcular intervals (see MATLAB plus function),
% it works on arrays of the same size or unit size along both dimensions.
% _________________________________________________________________________
% USAGE        
%   r = obj1 + obj2
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj1       : array of objects from the ciat.PolyarcularInterval class
%   obj2       : array of objects from the ciat.PolyarcularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________

    % Check input class
    mustBeA(obj1,["ciat.PolyarcularInterval","double"]);
    mustBeA(obj2,["ciat.PolyarcularInterval","double"]);
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % If the first object is a double make that the second
    if isa(obj1, 'double')
        objTemp = obj1;
        obj1 = obj2;
        obj2 = objTemp;
    end
            
    % Loop throught the arrays
    r(M,N) = ciat.PolyarcularInterval;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            if isa(obj2,'double')
                r(M,N) = obj1(m1,n1);
                r.ArcStorage{M,N}.Center = r.ArcStorage{M,N}.Center + ...
                                            obj2(m2,n2);
            else
                if obj1(m1,n1).isconvex && obj2(m2,n2).isconvex
                    r(M,N) = addConvex( obj1(m1,n1) , obj2(m2,n2) );
                else
                    r(M,N) = addConcave( obj1(m1,n1) , obj2(m2,n2) );
                end
            end
        end
    end
end

%% Function for adding two concave polyarcs

function r = addConcave(obj1,obj2)
    
    % Extract curve segments by type
    %   - extract arcs including vertices
    arc1 = [obj1.Arcs{:} ; obj1.Vertices{:}];
    arc2 = [obj2.Arcs{:} ; obj2.Vertices{:}];
    %   - extract edges with non-zero length
    edge1 = obj1.Edges{:};
    edge2 = obj2.Edges{:};    
    
    % Add arcs and vertices
    arcPlusArc = arc1 + arc2.';
    arcPlusEdge = arc1 + edge2.';
    edgePlusArc = edge1 + arc2.';
    
    % Extract valid segments
    arc3 = arcPlusArc(~isnan(arcPlusArc));
    edge3 = [ arcPlusEdge(~isnan(arcPlusEdge)) ; ...
              edgePlusArc(~isnan(edgePlusArc)) ];

    % Extract non-vertex segments
    arc3 = arc3(abs(arc3.Length)>10*eps);
    edge3 = edge3(abs(edge3.Length)>10*eps);

    % Split segments
    [arc3,edge3] = ciat.PolyarcularInterval.splitSegments(arc3,edge3);
    arc3 = arc3(abs(arc3.Length)>10*eps); % This should be unnecessary
    edge3 = edge3(abs(edge3.Length)>10*eps);
    
    % Trim segments
    arc3 = ciat.PolyarcularInterval.trimSegments(arc3,edge3);
    
    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3);

end


%% Function for adding two convex polyarcs

function r = addConvex(obj1,obj2)

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
        arc3(n3) = arc1(n1) + arc2(n2); % add what we found

        if  ang1(n1) < ang2(n2) + eps10 
            n1 = n1 + 1;
        elseif ang1(n1) > ang2(n2) + eps10
            n2 = n2 + 1;
        else
            n1 = n1 + 1;
            n2 = n2 + 1;
        end            
    end

    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3(1:n3));
    
end


%% Function for splitting arcs with angle interval crossing pi

function arcOut = splitArc(arcIn)
    % Find arcIn with angle interval including -pi or pi
    mask = arcIn.ArcAngle.isin(-pi) | arcIn.ArcAngle.isin(pi);

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
end
