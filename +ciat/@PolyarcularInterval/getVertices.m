% function [vertices,convexity] = getVertices(arcs,edges)
function [vertices,varargout] = getVertices(obj)

    % Extract elements
    arcs = obj.Arcs{:};
    edges = obj.Edges{:};
    arcStart = arcs.Startpoint;
    arcEnd = arcs.Endpoint;
    arcRadius = arcs.Radius;
    arcGauss = arcs.GaussMap;
    arcStartAngle = (arcRadius>0) .* arcGauss.Infimum + ...
                    (arcRadius<0) .* arcGauss.Supremum;
    arcEndAngle = (arcRadius>0) .* arcGauss.Supremum + ...
                  (arcRadius<0) .* arcGauss.Infimum;
    edgeStart = edges.Startpoint;
    edgeEnd = edges.Endpoint;
    edgeAngle = edges.GaussMap.Infimum;

    % Find vertices ...
    center = [];
    angStart = [] ;
    angEnd = [];
        % at the intersection of arc pairs
    [M,N] = find(abs(arcEnd - arcStart.')<10*eps);
    center = [center ; arcEnd(M)];
    angStart = [angStart ; arcEndAngle(M)];
    angEnd = [angEnd ; arcStartAngle(N)];
        % at the intersection of edge pairs
    [M,N] = find(abs(edgeEnd - edgeStart.')<10*eps);
    center = [center ; edgeEnd(M)];
    angStart =  [angStart ; edgeAngle(M)];
    angEnd = [angEnd ; edgeAngle(N)];
        % at the intersection of arc-edge pairs
    [M,N] = find(abs(arcEnd - edgeStart.')<10*eps);
    center = [center ; arcEnd(M)];
    angStart = [angStart ; arcEndAngle(M)];
    angEnd = [angEnd ; edgeAngle(N)];
        % at the intersection of edge-arc pairs
    [M,N] = find(abs(edgeEnd - arcStart.')<10*eps);
    angStart = [angStart ; edgeAngle(M)];
    angEnd = [angEnd ;arcStartAngle(N)];
    center = [center ; edgeEnd(M)];

    % Set the angular interval of the vertices and check convexity
    angDif = wrapToPi(angEnd - angStart);
    angInf = (angDif>=0) .* angStart + (angDif<0) .* angEnd;
    angSup = angInf + sign(angDif) .* angDif;
    
    % Return vertices
    vertices = ciat.Arc( center, zeros(length(center),1),...
                          ciat.RealInterval(angInf,angSup) );

    % Return convexity if requested
    if nargout > 1
        varargout = {angDif >= 0};
    end

end