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
    
    % Turn scalars to degenerate intervals
    if isa(obj1, 'double')
        obj1 = ciat.PolyarcularInterval(obj1);
    end
    if isa(obj2, 'double')
        obj2 = ciat.PolyarcularInterval(obj2);
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
            r(M,N) = add( obj1(m1,n1) , obj2(m2,n2) );
        end
    end
end

%% Function for adding two concave polyarcs

function r = add(obj1,obj2)
    
    % Extract curve segments by type
    %   - extract arcs with non-zero radius
    arc1 = obj1.Arcs(abs([obj1.Arcs.Radius])>0);
    arc2 = obj2.Arcs(abs([obj2.Arcs.Radius])>0);
    %   - extract edges with non-zero length
    edge1 = obj1.Edges([obj1.Edges.Length]>0);
    edge2 = obj2.Edges([obj2.Edges.Length]>0);    
    %   - extract all vertices including zero arcs and zero edges
    vert1 = [obj1.Vertices ; ...
             obj1.Arcs(abs([obj1.Arcs.Radius])==0) ; ...
             obj1.Edges([obj1.Edges.Length]==0)];
    vert2 = [obj2.Vertices ; ...
             obj2.Arcs(abs([obj2.Arcs.Radius])==0) ; ...
             obj2.Edges([obj2.Edges.Length]==0)];


    % Add arcs and vertices
    


    % Generate result arcs and vertices
    K = sum(mapMatchMat,'all');
    seg3(K,1) = ciat.Arc;
    k = 0;
    for n=1:N
        for m = find(mapMatchMat(n,:))
            k = k + 1;
            seg3(k).Center = seg1(n).Center + seg2(m).Center; 
            seg3(k).Radius = seg1(n).Radius + seg2(m).Radius;
            seg3(k).Angles = mapMatchVal(n,m);
        end
    end








    r = ciat.PolyarcularInterval(seg3);


end

%% Function for adding two convex polyarcs

function r = addConvex(obj1,obj2)
    
    % Create ordered lists of segments (arcs and vertices)
    seg1 = orderSegments(obj1);
    seg2 = orderSegments(obj2);

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


%% Function for creating an ordered set of arcs and vertices

function seg = orderSegments(obj)
    N = obj.ArcCount;
    M = sum([obj.Arcs.Radius]~=0);
    if N > 1 && M > 0
        seg(length(obj.Arcs)+length(obj.Vertices),1) = ciat.Arc;
        segIdx = 1;
        vertIdx = 1;
        for arcIdx = 1:length(obj.Arcs)
            if obj.Arcs(arcIdx).Radius ~= 0 
                seg(segIdx) = obj.Vertices(2*vertIdx-1);
                seg(segIdx+1) = obj.Arcs(arcIdx);
                seg(segIdx+2) = obj.Vertices(2*vertIdx);
                vertIdx = vertIdx + 1;
                segIdx = segIdx + 3;
            else
                seg(segIdx) = obj.Arcs(arcIdx);
                segIdx = segIdx + 1;
            end
        end
    else
        seg = obj.Arcs;
    end

    % Find and split segments with two Gauss intervals
    n = 1;
    while n<=length(seg)
        if length(seg(n).GaussMap)>1
            if n<length(seg)
                seg = [seg(1:n) ; seg(n) ; seg(n+1:end)];
            else
                seg = [seg(1:n) ; seg(n)];
            end
            seg(n+1).Angles = seg(n).GaussMap(2);
            seg(n).Angles = seg(n).GaussMap(1);
        end
        n=n+1;
    end

    % Shift order so the first segment is at -pi
    [~,n0] = min(inf([seg.GaussMap]));
    seg = circshift(seg,1-n0);



end