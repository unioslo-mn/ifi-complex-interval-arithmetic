function r = times(obj1,obj2)
   
% Element-wise multiplication of polyarcular intervals (.* operator)
%
% This function creates the polyarcular intervals representing the 
% element-wise product of two sets of polyarcular intervals 
% (see MATLAB times function).
% _________________________________________________________________________
% USAGE        
%   r = obj1 .* obj2
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
    
   
    % Loop throught the arrays
    r(M,N) = ciat.PolyarcularInterval;
    for m = 1:M
        for n = 1:N
            if isa(obj1,'double')
                r(m,n) = ciat.PolyarcularInterval(obj2.DefArcs .* obj1);
            elseif isa(obj2,'double')
                r(m,n) = ciat.PolyarcularInterval(obj1.DefArcs .* obj2);
            else
            % Calculate product
            r(M,N) = multiply( obj1(m,n) , obj2(m,n) );
            end
        end
    end
end

%%
function r = multiply(obj1,obj2)
    
      % Extract curve segments by type
    arc1 = obj1.Arcs;
    arc2 = obj2.Arcs;
    vert1 = obj1.Vertices;
    vert2 = obj2.Vertices;
    edge1 = obj1.Edges;
    edge2 = obj2.Edges; 

    % Multiply segments
    segments = cell(4,1);
    segments{1} = arc1 * arc2.';
    segments{2} = arc1 * edge2.';
    segments{3} = edge1 * arc2.';
    segments{4} = edge1 * edge2.';
    segments{5} = arc1 * vert2.';
    segments{6} = vert1 * arc2.';
    segments{7} = edge1 * vert2.';
    segments{8} = vert1 * edge2.';

    % Separate segments to arcs and edges
    arc3 = ciat.Arc;
    edge3 = ciat.Edge;
    for idx = 1:length(segments)
        switch class(segments{idx})
            case 'ciat.Arc'
                arc3 = [arc3;segments{idx}(~isnan(segments{idx}))];
            case 'ciat.Edge'
                edge3 = [edge3;segments{idx}(~isnan(segments{idx}))];
            case 'cell'
                arc3 = [arc3;segments{idx}{1}(~isnan(segments{idx}{1}))];
                edge3 = [edge3;segments{idx}{2}(~isnan(segments{idx}{2}))];
        end
    end

    % Split and trim segments
    [arc3,edge3] = ciat.PolyarcularInterval.splitSegments(arc3,edge3);
    arc3 = ciat.PolyarcularInterval.trimSegments(arc3,edge3,...
                                                'attempts',15,...
                                                'tolerance',1e-6);
    
    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3);
    r = mergeArcs(r);
end