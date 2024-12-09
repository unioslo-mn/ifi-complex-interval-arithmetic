function r = times(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.Edge","double","ciat.Arc"]);
    mustBeA(obj2,["ciat.Edge","double","ciat.Arc"]);
    
    % Get input sizes and check if they can be combined
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);

    % If one is arc, forward to arc times edge
    if isa(obj2,'ciat.Arc')
        r = obj2 .* obj1;
        return
    end

    % Loop throught the arrays
    arcOut(M,N) = ciat.Arc;
    edgeOut(M,N) = ciat.Edge;
    for m = 1:M
        for n = 1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);

            % One of the operands is a double
            if isa(obj1,"ciat.Edge") && isa(obj2,"double")
                segment= timesDouble(obj1(m1,n1),obj2(m2,n2));
            elseif isa(obj1,"double") && isa(obj2,"ciat.Edge")
                segment= timesDouble(obj2(m2,n2),obj1(m1,n1));
            else
                segment= edgeTimesEdge(obj1(m1,n1),obj2(m2,n2));
            end

            % Assign the result segment to the array of the appropritate type
            if ~isnan(segment)
                if isa(segment,'ciat.Arc')
                    arcOut(m,n) = segment;
                else
                    edgeOut(m,n) = segment;
                end
            end
        end
    end

    if all(isnan(arcOut),'all')
        r = edgeOut;
    elseif all(isnan(edgeOut),'all')
        r = arcOut;
    else
        r = {arcOut , edgeOut};
    end
end

%% Multiplication with a double

function r = timesDouble(edge,x)
    r = ciat.Edge(edge.Startpoint * x,...
                  edge.Endpoint * x );
end


%% Multiplication with an edge

function r = edgeTimesEdge(edge1,edge2)
    
    % Calculate log-Gauss map intersection
    LGMcap = ciat.Arc.capGaussMap(edge1.LogGaussMap,edge2.LogGaussMap);

    if ~isnan(LGMcap)
        if edge1.ZeroCrossing == 0 && edge2.ZeroCrossing == 0
            % This product is an edge between the product of the endpoints.
            r = ciat.Edge(edge1.Startpoint .* edge2.Startpoint,...
                          edge1.Endpoint .* edge2.Endpoint);
        elseif edge1.ZeroCrossing == 0
            % The zero-crossing edge's LGM is +/- pi, so the product is
            % this edge multiplied by the point in the other edge with the
            % corresponding LGM value
            if edge2.LogGaussMap.isin(edge1.LogGaussMap.mid)
                r = edge1 .* edge2.findLGM(edge1.LogGaussMap.mid);
            else
                r = NaN;
            end
        elseif edge2.ZeroCrossing == 0
            % The zero-crossing edge's LGM is +/- pi, so the product is
            % this edge multiplied by the point in the other edge with the
            % corresponding LGM value
            if edge1.LogGaussMap.isin(edge2.LogGaussMap.mid)
                r = edge2 .* edge1.findLGM(edge2.LogGaussMap.mid);
            else
                r = NaN;
            end
        else
            % The result is an arc fitted to an ellipse or hyperbole
            % segment
            r = fitArcToEdgeTimesEdge(edge1,edge2);
        end
    else
        r = NaN;
    end
end


%% Fit arc to the product of two edges
function r = fitArcToEdgeTimesEdge(edge1,edge2)
    % Extract radii and curve parameters
            S = edge1.CurveParameter;
            T = edge2.CurveParameter;
    
            % Define envelope function and touching condition
            H_st = @(s,t) 1-s.*t + 1i*(s+t);
            j_st = @(s,t) s - t;
    
            % Find endpoints of envelope segment 
            sJinf = fsolve(@(s) j_st(s,T.inf), S.mid,optimset('Display','off'));
            sJsup = fsolve(@(s) j_st(s,T.sup), S.mid,optimset('Display','off'));
            tJinf = fsolve(@(t) j_st(S.inf,t), T.mid,optimset('Display','off'));
            tJsup = fsolve(@(t) j_st(S.sup,t), T.mid,optimset('Display','off'));
            pHst = [];
            if S.isin(sJinf); pHst = [pHst ; [sJinf,T.inf]];end
            if S.isin(sJsup); pHst = [pHst ; [sJsup,T.sup]];end
            if T.isin(tJinf); pHst = [pHst ; [S.inf,tJinf]];end
            if T.isin(tJsup); pHst = [pHst ; [S.sup,tJsup]];end
            if ~all(size(pHst) == [2,2])
                r = ciat.Arc;
                return
            end
            pHxy = H_st(pHst(:,1),pHst(:,2));
    
            % Sample then envelope
            figure;
            J_smp = fimplicit(j_st,[S.inf S.sup T.inf T.sup]);
            J_smp = [J_smp.XData;J_smp.YData]';
            close
            H_smp = H_st(J_smp(:,1),J_smp(:,2));
            if isempty(H_smp)
                r = ciat.Arc;
                return
            end
    
            % Fit arc to the envelope
            cCnt = @(u) u * exp(1i*(angle(diff(pHxy))+pi/2)) + mean(pHxy);
            cRad = @(u) abs(cCnt(u)-pHxy(1));
            cAng1 = @(u) angle(pHxy(1)-cCnt(u));
            cAng2 = @(u) angle(pHxy(2)-cCnt(u));
            cFit = @(u) abs(min(abs(H_smp-cCnt(u))-cRad(u)));
            uSol = fminsearch(cFit,1);

            % Create concave arc
            arcCenter = cCnt(uSol);
            arcRadius = - cRad(uSol);
            arcAngle = ciat.RealInterval(cAng1(uSol), cAng2(uSol)) - pi;
            
            % Flip the arc-angle interval if necessary
            cAng0 = angle(H_smp(ceil(length(H_smp)/2)) - arcCenter) - pi ;
            if ~any(arcAngle.isin(cAng0+[-2,0,2]*pi))
                arcAngle = ciat.RealInterval(cAng2(uSol)-2*pi,cAng1(uSol));
            end

            % Assign arc and reverse normalization
            r = ciat.Arc(arcCenter , arcRadius , arcAngle);
            r = r * (1 / edge1.NormFactor / edge2.NormFactor);

            
            % 
            % r = ciat.Arc(cCnt(uSol),cRad(uSol),...
            %         ciat.RealInterval(cAng1(uSol),cAng2(uSol)));
            % 
            % % Flip the arc-angle interval if necessary
            % cAng0 = angle(H_smp(ceil(length(H_smp)/2)) - r.Center);
            % if ~r.ArcAngle.isin(cAng0)
            %     r = ciat.Arc(cCnt(uSol),cRad(uSol),...
            %         ciat.RealInterval(cAng2(uSol)-2*pi,cAng1(uSol)));
            % end
            % 
            % % Reverse normalization
            % r = r * (1 / edge1.NormFactor / edge2.NormFactor);
end
