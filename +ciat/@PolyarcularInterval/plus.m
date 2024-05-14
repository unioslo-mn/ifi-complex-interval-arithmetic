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
    
    % % Turn scalars to degenerate intervals
    % if isa(obj1, 'double')
    %     obj1 = ciat.PolyarcularInterval(obj1);
    % end
    % if isa(obj2, 'double')
    %     obj2 = ciat.PolyarcularInterval(obj2);
    % end 

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
                r(M,N).ArcStorage.Center = r(M,N).ArcStorage.Center + ...
                                            obj2(m2,n2);
            else
                r(M,N) = addConcaveArcs( obj1(m1,n1) , obj2(m2,n2) );
            end
        end
    end
end

%% Function for adding two concave polyarcs

function r = addConcaveArcs(obj1,obj2)
    
    % Extract curve segments by type
    %   - extract arcs including vertices
    arc1 = [obj1.Arcs ; obj1.Vertices];
    arc2 = [obj2.Arcs ; obj2.Vertices];
    %   - extract edges with non-zero length
    edge1 = obj1.Edges;
    edge2 = obj2.Edges;    
    
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
    
    % Create ordered lists of segments (arcs and vertices)
    seg1 = orderSegments(obj1);
    seg2 = orderSegments(obj2);aIsum.plot('k-','linewidth',2)

    % Create match matrix
    N = length(seg1);
    M = length(seg2);
    segMatch = zeros(N,M);
    segMatchIntv(N,M) = ciat.RealInterval;
    for n = 1:N
        for m = 1:M
            segCap = intersection([seg1(n).GaussMap,seg2(m).GaussMap]);
            segMatch(n,m) = ~isempty(segCap);
            if segMatch(n,m)
                segMatchIntv(n,m) = segCap;
            end
        end
    end
    seg1.GaussMap;

    % Generate result arcs
    arcs = [];
    currGauss = -pi;
    for n=1:N
        while sum(segMatch(n,:))
            infGauss = [segMatchIntv(n,:).Infimum];
            m = find(infGauss==currGauss);
            arcs = [arcs; ciat.Arc(seg1(n).Center + seg2(m).Center,...
                                   seg1(n).Radius + seg2(m).Radius,...
                                   segMatchIntv(n,m)) ];
            currGauss = segMatchIntv(n,m).Supremum;
            segMatch(n,m) = 0;
        end
    end


    r = ciat.PolyarcularInterval(arcs);


end

