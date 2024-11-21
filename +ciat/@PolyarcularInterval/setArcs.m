function arcs = setArcs(arcs)
    N = length(arcs);
    if N > 1
        for n=1:N
            % Fix the angles property of zero radius arcs
            if arcs(n).Radius == 0
                % set current, previous and next arc
                arcCurr = arcs(n);
                if n > 1
                    arcPrev = arcs(n-1);
                else
                    arcPrev = arcs(N);
                end
                if n < N
                    arcNext = arcs(n+1);
                else
                    arcNext = arcs(1);
                end
                
                angMin = angle(arcCurr.Center - ...
                               arcPrev.Endpoint) - pi/2;
                angMax = angle(arcNext.Startpoint - ...
                               arcCurr.Center) - pi/2;
                if angMin > angMax
                    angMax = angMax + 2*pi;
                end
                arcs(n).ArcAngle = ciat.RealInterval(angMin,angMax);
            end
        end
    % else
        % arcs.ArcAngle = ciat.RealInterval(-pi,pi);
    end
end