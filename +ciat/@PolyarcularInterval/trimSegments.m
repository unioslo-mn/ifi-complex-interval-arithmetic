function arcOut = trimSegments(arcIn,edgeIn,optional)


arguments
    arcIn            (:,:)  {mustBeA(arcIn,{'ciat.Arc','double'})}  = []
    edgeIn           (:,:)  {mustBeA(edgeIn,{'ciat.Edge','double'})}  = []
    optional.inner   (1,1)  {mustBeNumericOrLogical} = false
    optional.attempts (1,1) {mustBeNumeric} = 5
    optional.tolerance (1,1) {mustBeNumeric} = 100*eps
end 

% Remove dead-end segments



for iTry = 1:optional.attempts

    % Initialize
    K = length(arcIn) + length(edgeIn);
    arcOut(K,1) = ciat.Arc;
    arcIdx = zeros(K,1);
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
    while k==1 || (k<=K && idx~=startIdx && all(idx~=prevIdx))
    
        % Select segment of the given index and store defining arc if available
        if idx <= length(arcIn)
            seg = arcIn(idx);
            arcOut(k) = arcIn(idx);
            arcIdx(k) = idx;
            edgeFlag = false;
        else
            seg = edgeIn(idx-length(arcIn));
            if edgeFlag
                arcOut(k) = ciat.Arc(seg.Startpoint,0,0);
                arcIdx(k) = idx;
            end
            edgeFlag = true;
        end
    
        % Find next segment
        prevIdx = [prevIdx;idx];
        idx = find( abs(seg.Endpoint - startPoints) < optional.tolerance );
        if isempty(idx)
            if iTry < optional.attempts
                % Remove element and continue
                if prevIdx(end) <= length(arcIn)
                    arcIn = arcIn((setdiff(1:end,prevIdx(end))));
                else
                    edgeIn = edgeIn((setdiff(1:end,prevIdx(end)-length(arcIn))));
                end    
                break
            else
                warning(sprintf('Boundary incontinuity after %i attempts, returning NaN',optional.attempts))
                arcOut(1,1) = ciat.Arc;
                return
            end
            else
                
        end
        if length(idx) > 1
            if any(idx == startIdx)
                idx = startIdx;
            else
                % Select the segment with the smaller initial Gauss angle
                diffGauss = ciat.wrapToPi(startGauss(idx)-endGauss(prevIdx(end)));
                if optional.inner
                    idx = idx(diffGauss == max(diffGauss));
                else
                    idx = idx(diffGauss == min(diffGauss));
                end
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
end

% For inner boundary trimming, remove the first few segments
if any(idx==prevIdx)
    arcOut = arcOut(find(arcIdx==idx):end);
end

% Check overrun condition
if k > K && idx ~= startIdx && all(idx~=prevIdx)
    error('Boundary is not closed')
end

% Add last vertex
if edgeFlag && idx>length(arcIn)
    arcOut(k) = ciat.Arc(seg.Endpoint,0,0);
end

arcOut = arcOut(~isnan(arcOut));

end