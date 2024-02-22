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
    r(M,N) = ciat.PolygonalInterval;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            r(M,N) = add( obj1(m1,n1) , obj2(m2,n2) );
            r(M,N).Points = r(M,N).Boundary;
        end
    end
end

%% Function for adding two polygons

function r = add(obj1,obj2)
    
    % Create ordered lists of segments (arcs and vertices)
    seg1 = orderSegments(obj1);
    seg2 = orderSegments(obj2);

    % Create match matrix
    N = length(seg1);
    M = length(seg2);
    matchMat = zeros(N,M);
    for n = 1:N
        for m = n:M
            matchMat(n,m) = isempty(intersection([seg1(n).GaussMap,...
                                                  seg2(m).GaussMap]));
        end
    end
    seg1.GaussMap;

    % % Gauss map match segments
    % arcs = [];
    % idx1 = 1;
    % idx2 = 1;
    % while idx1 <= length(seg1) && idx2 <= length(seg2)
    %     % Intersect Gauss maps
    %     match = intersection([seg1(idx1).GaussMap,...
    %                           seg2(idx2).GaussMap]);
    % 
    % 
    %     if ~isempty(match)
    % 
    %     else
    %         idx2 = idx2 + 1;
    %     end
    % 
    % end



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
            seg(segIdx) = obj.Arcs(arcIdx);
            if obj.Arcs(arcIdx).Radius ~= 0 
                seg(segIdx+1) = obj.Vertices(2*vertIdx-1);
                seg(segIdx+2) = obj.Vertices(2*vertIdx);
                vertIdx = vertIdx + 1;
                segIdx = segIdx + 3;
            else
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
end