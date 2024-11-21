function r = plusConcave(obj1,obj2)
    
    % Extract curve segments by type
    %   - extract arcs including vertices
    arc1 = [obj1.Arcs ; obj1.Vertices];
    arc2 = [obj2.Arcs ; obj2.Vertices];
    %   - extract edges with non-zero length
    edge1 = obj1.Edges;
    edge2 = obj2.Edges;    
    
    % Add arcs
    arc3 = ciat.Arc;
    if ~isempty(arc1) && ~isempty(arc2)
        arc3 = arc1 + arc2.';
        arc3 = arc3(~isnan(arc3));
        arc3 = arc3(abs(arc3.Length)>100*eps);
    end

    % Add arcs and edges
    edge3 = ciat.Edge;
    if ~isempty(arc1) && ~isempty(edge2)
        arcPlusEdge = arc1 + edge2.';
        edge3 = [edge3 ; arcPlusEdge(~isnan(arcPlusEdge))];
    end
    if ~isempty(edge1) && ~isempty(arc2)
        edgePlusArc =  edge1 + arc2.';
        edge3 = [edge3 ; edgePlusArc(~isnan(edgePlusArc))];
    end
    edge3 = edge3(abs(edge3.Length)>100*eps);

    % Split and trim segments
    [arc3,edge3] = ciat.PolyarcularInterval.splitSegments(arc3,edge3);
    arc3 = ciat.PolyarcularInterval.trimSegments(arc3,edge3);
    
    % Generate polyarc
    r = ciat.PolyarcularInterval(arc3);
    r = mergeArcs(r);

end



