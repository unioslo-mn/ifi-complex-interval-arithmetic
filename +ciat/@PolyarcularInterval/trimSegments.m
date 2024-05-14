function arcOut = trimSegments(arcIn,edgeIn)

mustBeA(arcIn,'ciat.Arc')
mustBeA(edgeIn,'ciat.Edge')

% Initialize
K = length(arcIn) + length(edgeIn);
arcOut(K,1) = ciat.Arc;
startPoints = [arcIn.Startpoint ; edgeIn.Startpoint];
midGauss = [arcIn.GaussMap.Midpoint ; edgeIn.GaussMap.Midpoint];
startGauss = [arcIn.GaussMap.Infimum .* (arcIn.Radius > 0) + ...
              arcIn.GaussMap.Supremum .* (arcIn.Radius < 0) ; ... 
              edgeIn.GaussMap.Midpoint];
endGauss = [arcIn.GaussMap.Supremum .* (arcIn.Radius > 0) + ...
              arcIn.GaussMap.Infimum .* (arcIn.Radius < 0) ; ... 
              edgeIn.GaussMap.Midpoint];

% Find starting object
[~,idx] = min([arcIn.Real.Infimum ; edgeIn.Real.Infimum]);
startIdx = idx;

% Follow the boundary
k = 1;
edgeFlag = false;
while k==1 || (k<=K && idx~=startIdx)

    % Select segment of the given index and store defining arc if available
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

    % Find next segment
    prevIdx = idx;
    idx = find( abs(seg.Endpoint - startPoints) < 10*eps );
    if isempty(idx)
        error('Boundary incontinuity')
    end
    if length(idx) > 1
        [~,minIdx] = min(wrapToPi(startGauss(idx)-endGauss(prevIdx)));
        idx = idx(minIdx);
    end

    % Increment counter
    k = k+1;
end

% Check overrun condition
if k > K && idx ~= startIdx
    error('Boundary is not closed')
end

if edgeFlag && idx>length(arcIn)
    arcOut(k) = ciat.Arc(seg.Endpoint,0,0);
end

arcOut = arcOut(~isnan(arcOut));

end