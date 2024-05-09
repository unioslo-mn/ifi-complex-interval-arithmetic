function arcOut = trimSegments(arcIn,edgeIn)

mustBeA(arcIn,'ciat.Arc')
mustBeA(edgeIn,'ciat.Edge')

% Initialize
K = length(arcIn) + length(edgeIn);
arcOut(K,1) = ciat.Arc;
startPoints = [arcIn.Startpoint ; edgeIn.Startpoint];
midGauss = [arcIn.GaussMap.Midpoint ; edgeIn.GaussMap.Midpoint];

% Find starting object
[~,idx] = min([arcIn.Real.Infimum ; edgeIn.Real.Infimum]);
startIdx = idx;

% Follow the boundary
k = 1;
edgeFlag = false;
while k <= K
    if idx <= length(arcIn)
        seg = arcIn(idx);
        arcOut(k) = arcIn(idx);
        edgeFlag = false;
    else
        seg = edgeIn(idx-length(arcIn));
        if edgeFlag
            arcOut(k) = ciat.Arc(seg.Startpoint,0,0);
        end
        edgeFlag = true;
    end

    idx = find( abs(seg.Endpoint - startPoints) < 10*eps );

    if idx == startIdx
        break
    end

    if isempty(idx)
        error('Boundary incontinuity')
    end
    if length(idx) > 1
        [~,minIdx] = min(midGauss(idx));
        idx = idx(minIdx);
    end
    k = k+1;
end

if k > K
    error('Boundary is not closed')
end

if edgeFlag
    arcOut(k) = ciat.Arc(seg.Startpoint,0,0);
end
arcOut = arcOut(~isnan(arcOut));

end