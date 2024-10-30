function r = timesPolar(xObj,pObj)

    % Extract arcs and vertices
    xArc = [xObj.Arcs{:} ; xObj.Vertices{:}];
    paObj = ciat.PolyarcularInterval(pObj);
    pArc = paObj.Arcs{:}(2);
    pVert = paObj.Vertices{:};
    
    % Split arcs at the log Gauss map quadrant limits (-pi/2,0,pi/2,pi)
    xArc = splitArcsAtLogGauss(xArc);
    
    % Multiply with the polar vertices
    rArcs = ciat.Arc;
    for n=1:4
        gMask = ~isnan(cap(xArc.LogGaussMap,pVert(n).LogGaussMap.'));
        rArcs = [rArcs ; xArc(gMask) .* pVert(n).Center] ;
    end
        
    % Find outermost point on the radial axis
    xArcMax = xArc(xArc.LogGaussMap.isin(0));
    if xArcMax.Radius > 0
        if xArcMax.Center == 0
            pArcMax = ciat.Arc( 0 , pArc.Radius*xArcMax.Radius , ...
                                    pArc.ArcAngle + xArcMax.ArcAngle);
        else
            pArcMax = pArc .* (abs(xArcMax.Center) + xArcMax.Radius) * ...
                                   exp(1i*angle(xArcMax.Center));
        end
    else
        pArcMax = pArc .* xArcMax.Center;
    end
    rArcs = [rArcs ; pArcMax ];

    % Generate output shape
    r = ciat.PolyarxInterval(ciat.PolyarcularInterval(rArcs));

end


%% Utility functions

function xArc = splitArcsAtLogGauss(xArc)

    % Split at max radius
    N = length(xArc);
    rMask = xArc.Radius > 0;
    gMask = any(xArc.LogGaussMap.isinside([-2*pi,0,2*pi]),2);
    idx = find(rMask & gMask);
    if ~isempty(idx)
        c = xArc(idx).Center;
        r = xArc(idx).Radius;
        alpha = angle(xArc(idx).Center);
        pnt = c + r * exp(1i*alpha);
        if N>1
            xArc = [xArc(setdiff(1:N,idx)) ; xArc(idx).split(pnt)];
        else
            xArc = xArc(idx).split(pnt);
        end
    else
        pnt = xArc(rMask & gMask).Center;
    end
        % Split at min radius
    N = length(xArc);
    rMask = xArc.Radius > 0;
    gMask = any(xArc.LogGaussMap.isinside([-pi,pi]),2);
    idx = find(rMask & gMask);
    if ~isempty(idx)
        c = xArc(idx).Center;
        r = xArc(idx).Radius;
        alpha = angle(xArc(idx).Center);
        pnt = c - r * exp(1i*alpha);
        if N>1
            xArc = [xArc(setdiff(1:N,idx)) ; xArc(idx).split(pnt)];
        else
            xArc = xArc(idx).split(pnt);
        end
    else
        pnt = xArc(rMask & gMask).Center;
    end
        % Split at max angle
    N = length(xArc);
    rMask = xArc.Radius > 0;
    gMask = any(xArc.LogGaussMap.isinside([-3*pi/2,pi/2]),2); 
    idx = find(rMask & gMask);
    if ~isempty(idx)
        c = xArc(idx).Center;
        R = abs(c);
        alpha = angle(c);
        r = xArc(idx).Radius;
        pnt = sqrt(R^2-r^2) * exp(1j*(alpha + asin(r/R)));
        if N>1
            xArc = [xArc(setdiff(1:N,idx)) ; xArc(idx).split(pnt)];
        else
            xArc = xArc(idx).split(pnt);
        end
    else
        pnt = xArc(rMask & gMask).Center;
    end
        % Split at min angle
    N = length(xArc);
    rMask = xArc.Radius > 0;
    gMask = any(xArc.LogGaussMap.isinside([-pi/2,3*pi/2]),2); 
    idx = find(rMask & gMask);
    if ~isempty(idx)
        R = abs(xArc(idx).Center);
        r = xArc(idx).Radius;
        alpha = angle(xArc(idx).Center);
        pnt = sqrt(R^2-r^2) * exp(1j*(alpha - asin(r/R)));
        if N>1
            xArc = [xArc(setdiff(1:N,idx)) ; xArc(idx).split(pnt)];
        else
            xArc = xArc(idx).split(pnt);
        end
    else
        pnt = xArc(rMask & gMask).Center;
    end

end



