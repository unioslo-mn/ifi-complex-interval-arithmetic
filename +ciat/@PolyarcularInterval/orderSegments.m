function seg = orderSegments(obj)
    N = obj.ArcCount;
    M = sum([obj.Arcs.Radius]~=0);
    if N > 1 && M > 0
        seg(length(obj.Arcs)+length(obj.Vertices),1) = ciat.Arc;
        segIdx = 1;
        vertIdx = 1;
        for arcIdx = 1:length(obj.Arcs)
            if obj.Arcs(arcIdx).Radius ~= 0 
                seg(segIdx) = obj.Vertices(2*vertIdx-1);
                seg(segIdx+1) = obj.Arcs(arcIdx);
                seg(segIdx+2) = obj.Vertices(2*vertIdx);
                vertIdx = vertIdx + 1;
                segIdx = segIdx + 3;
            else
                seg(segIdx) = obj.Arcs(arcIdx);
                segIdx = segIdx + 1;
            end
        end
    else
        seg = obj.Arcs;
    end

    % Find and split segments with two Gauss intervals
    n = 1;
    while n<=length(seg)
        if length(seg(n).GaussMap)>1
            if n<length(seg)
                seg = [seg(1:n) ; seg(n) ; seg(n+1:end)];
            else
                seg = [seg(1:n) ; seg(n)];
            end
            seg(n+1).Angles = seg(n).GaussMap(2);
            seg(n).Angles = seg(n).GaussMap(1);
        end
        n=n+1;
    end

    % Shift order so the first segment is at -pi
    [~,n0] = min(inf([seg.GaussMap]));
    seg = circshift(seg,1-n0);
end