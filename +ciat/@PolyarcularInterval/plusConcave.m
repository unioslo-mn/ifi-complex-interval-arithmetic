function r = plusConcave(obj1,obj2)
    
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
    arc3 = arc3(abs(arc3.Length)>100*eps);
    edge3 = edge3(abs(edge3.Length)>100*eps);

    % Split and trim segments
    [arc3,edge3] = ciat.PolyarcularInterval.splitSegments(arc3,edge3);
    arc3 = ciat.PolyarcularInterval.trimSegments(arc3,edge3,1);
    
    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3);
    r = joinSegments(r);

end



