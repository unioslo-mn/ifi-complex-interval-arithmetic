function [arcOut,edgeOut] = splitSegments(arcIn,edgeIn)

    % Split segments at intersection points
    arcBox = ciat.RectangularInterval(arcIn);
    edgeBox = ciat.RectangularInterval(edgeIn);
        % Split arcs
    arcOut = [];
    for k = 1:length(arcIn)
        if length(arcIn)>1
            capBox = cap([arcBox(setdiff(1:end,k));edgeBox],arcBox(k));
        else
            capBox = cap(edgeBox,arcBox(k));
        end
        boxMask = ~isnan(capBox) & capBox.Area>0;
        if any(boxMask,'all')
            boxIdx = find(boxMask);
            splitPoint = [];
            for l = 1:length(boxIdx)
                if boxIdx(l) < k
                    seg = arcIn(boxIdx(l));
                elseif boxIdx(l) > length(arcIn)-1
                    seg = edgeIn(boxIdx(l)-length(arcIn)+1);
                else
                    seg = arcIn(boxIdx(l)+1);
                end
                splitPoint = [splitPoint ; cap(arcIn(k),seg)];
            end
            arcOut = [arcOut ; arcIn(k).split(splitPoint)];
        else
            arcOut = [arcOut;arcIn(k)];
        end
    end
    
    % Split edges
    edgeOut = [];
    for k = 1:length(edgeIn)
        if length(edgeIn)>1
            capBox = cap([arcBox;edgeBox(setdiff(1:end,k))],edgeBox(k));
        else
            capBox = cap(arcBox,edgeBox(k));
        end
        boxMask = ~isnan(capBox) & capBox.Area>0;
        if any(boxMask,'all')
            boxIdx = find(boxMask);
            splitPoint = [];
            for l = 1:length(boxIdx)
                if boxIdx(l) <= length(arcIn)
                    seg = arcIn(boxIdx(l));
                elseif boxIdx(l) - length(arcIn) > k
                    seg = edgeIn(boxIdx(l)-length(arcIn)+1);
                else
                    seg = edgeIn(boxIdx(l)-length(arcIn));
                end
                splitPoint = [splitPoint ; cap(edgeIn(k),seg)];
            end

            edgeOut = [edgeOut ; edgeIn(k).split(splitPoint)];
        else
            edgeOut = [edgeOut;edgeIn(k)];
        end
    end

    % This should not be necessary
    arcOut = arcOut(abs(arcOut.Length)>10*eps); 
    edgeOut = edgeOut(abs(edgeOut.Length)>10*eps);

end