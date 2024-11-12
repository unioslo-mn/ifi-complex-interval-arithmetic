function arcOut = trimSegments(arcIn,edgeIn,recurCount)

mustBeA(arcIn,{'ciat.Arc','double'})
mustBeA(edgeIn,{'ciat.Edge','double'})

% Initialize
K = length(arcIn) + length(edgeIn);
arcOut(K,1) = ciat.Arc;
startPoints = [arcIn.Startpoint ; edgeIn.Startpoint];
startGauss = [arcIn.GaussMap.Infimum .* (arcIn.Radius > 0) + ...
              arcIn.GaussMap.Supremum .* (arcIn.Radius < 0) ; ... 
              edgeIn.GaussMap.Midpoint];
endGauss = [arcIn.GaussMap.Supremum .* (arcIn.Radius > 0) + ...
              arcIn.GaussMap.Infimum .* (arcIn.Radius < 0) ; ... 
              edgeIn.GaussMap.Midpoint];
segLength = [arcIn.Length ; edgeIn.Length];
segRealInf = [inf(real(arcIn)) ; inf(real(edgeIn))];

% Find starting object
idx = find(segRealInf == min(segRealInf),1);
startIdx = idx;
if idx <= length(arcIn)
    edgeFlag = false;
else
    edgeFlag = true;
end

% Follow the boundary
k = 1;
prevIdx = [];
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
    prevIdx = [prevIdx;idx];
    idx = find( abs(seg.Endpoint - startPoints) < 100*eps );
    if isempty(idx)
        if recurCount < 5
            % Remove element and try trimming again (recursive hell)
            if prevIdx(end) <= length(arcIn)
                arcIn = arcIn((setdiff(1:end,prevIdx(end))));
            else
                edgeIn = edgeIn((setdiff(1:end,prevIdx(end)-length(arcIn))));
            end
            % warning('Boundary incontinuity, removing last segment and trying again')
            arcOut = ciat.PolyarcularInterval.trimSegments(arcIn,edgeIn,...
                                                           recurCount+1);
            return
        else
            warning('Boundary incontinuity after 5 attempts, returning NaN')
            arcOut(1,1) = ciat.Arc;
            return
        end
    end
    if length(idx) > 1
        if any(idx == startIdx)
            idx = startIdx;
        else
            % Select the segment with the smaller initial Gauss angle
            diffGauss = ciat.wrapToPi(startGauss(idx)-endGauss(prevIdx(end)));
            idx = idx(diffGauss == min(diffGauss));
        end
        if length(idx) > 1
            % Select the segment with the smallest Gauss angle derivative
            diffGauss = ciat.wrapToPi(endGauss(idx)-startGauss(idx))./segLength(idx);
            minIdx = find(diffGauss == min(diffGauss),1);
            idx = idx(minIdx);
        end
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