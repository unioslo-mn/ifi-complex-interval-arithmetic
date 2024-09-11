function r = timesPolar(xObj,pObj)

    % Extract arcs and vertices
    xArcs = [xObj.Arcs{:} ; xObj.Vertices{:}];
    pxObj = ciat.PolyarxInterval(pObj);
    pArc = pxObj.Arcs{:};
    pVert = pxObj.Vertices{:};

    % Multiply with the polar vertices
    rArcs = ciat.Arc;
    for n=1:4
        gMask = ~isnan(cap(xArcs.LogGaussMap,pVert(n).LogGaussMap.'));
        rArcs = [rArcs ; xArcs(gMask) .* pVert(n).Center] ;
    end
        
    % Find outermost point on the radial axis
    xArcMax = xArcs(xArcs.LogGaussMap.isin(0));
    if xArcMax.Radius > 0
        if xArcMax.Center == 0
            pArcMax = ciat.Arc( 0 , pArc.Radius*xArcMax.Radius , ...
                                    pArc.ArcAngle + xArcMax.ArcAngle);
        else
            pArcMax = ciat.Arc .* (abs(xArcMax.Center) + xArcMax.Radius) * ...
                                   exp(1i*angle(xArcMax.Center));
        end
    else
        pArcMax = pArc .* xArcMax.Center;
    end
    rArcs = [rArcs ; pArcMax ];

    % Generate output shape
    r = ciat.PolyarxInterval(ciat.PolyarcularInterval(rArcs));

end