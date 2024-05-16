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
            splitPoint = splitPoint(~isnan(splitPoint));
            if ~isempty(splitPoint)
                arcCenter = arcIn(k).Center;
                arcRadius = arcIn(k).Radius;
                arcAngInf = arcIn(k).ArcAngle.Infimum;
                arcAngSup = arcIn(k).ArcAngle.Supremum;
                splitAngle = angle(splitPoint-arcCenter) + pi*(arcRadius<0);
                % splitAngle = sort(wrapToPi(splitAngle));
                splitAngle = wrapTo2Pi(splitAngle)-2*pi;
                splitAngle = splitAngle + 2*pi*(splitAngle<arcAngInf);
                splitAngle = sort(splitAngle);
                splitAngle = [arcAngInf ; splitAngle ; arcAngSup ];
                splitAngle = uniquetol(splitAngle,10*eps);
                for l = 1:length(splitAngle)-1
                    arcOut = [arcOut ; ...
                                ciat.Arc(arcCenter,arcRadius,...
                                         ciat.RealInterval(splitAngle(l), ...
                                                           splitAngle(l+1)))];
                end
            else
                arcOut = [arcOut;arcIn(k)];
            end
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
                elseif boxIdx(l) > k
                    seg = edgeIn(boxIdx(l)-length(arcIn)+1);
                else
                    seg = edgeIn(boxIdx(l)-length(arcIn));
                end
                splitPoint = [splitPoint ; cap(edgeIn(k),seg)];
            end
            splitPoint = splitPoint(~isnan(splitPoint));
            if ~isempty(splitPoint)
                p1 = edgeIn(k).Startpoint;
                p2 = edgeIn(k).Endpoint;
                [~,sortIdx] = sort(abs(splitPoint-p1));
                splitPoint = [p1 ; splitPoint(sortIdx); p2];
                [~,uniqueIdx,~] = uniquetol(abs(splitPoint-p1),10*eps);
                splitPoint = splitPoint(uniqueIdx);
                for l = 1:length(splitPoint)-1
                    edgeOut = [edgeOut ; ...
                                ciat.Edge(splitPoint(l),splitPoint(l+1))];
                end
            else
                edgeOut = [edgeOut;edgeIn(k)];
            end
        else
            edgeOut = [edgeOut;edgeIn(k)];
        end
    end

end