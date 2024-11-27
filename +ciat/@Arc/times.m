function r = times(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.Arc","ciat.Edge","double"]);
    mustBeA(obj2,["ciat.Arc","ciat.Edge","double"]);
    
    % Get input sizes and check if they can be combined
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);

    % Loop throught the arrays
    % r(M,N) = ciat.Arc;
    r = cell(M,N);
    for m = 1:M
        for n = 1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);

            % One of the operands is a double
            if isa(obj1,"ciat.Arc") && isa(obj2,"double")
                r{m,n} = timesDouble(obj1(m1,n1),obj2(m2,n2));
            elseif isa(obj1,"double") && isa(obj2,"ciat.Arc")
                r{m,n} = timesDouble(obj2(m2,n2),obj1(n1,m1));

            % One of the operands is an edge
            elseif isa(obj1,"ciat.Arc") && isa(obj2,"ciat.Edge")
                r{m,n} = timesEdge(obj1(m1,n1),obj2(m2,n2));
            elseif isa(obj1,"ciat.Edge") && isa(obj2,"ciat.Arc")
                r{m,n} = timesEdge(obj2(m2,n2),obj1(n1,m1));

            % Both operands are arcs
            elseif isa(obj1,"ciat.Arc") && isa(obj2,"ciat.Arc")
                r{m,n} = timesArc(obj1(m1,n1),obj2(m2,n2));
            end
        end
    end

    
end

%% Multiplication with a double

function r = timesDouble(arc,x)
    r = ciat.Arc(arc.Center .* x , ...
                 arc.Radius .* abs(x) , ...
                 arc.ArcAngle + angle(x));
end

%% Multiplication with an edge

function r = timesEdge(arc,edge)
    
    arcCenter = arc.Center;
    edgeZero = edge.ZeroCrossing;

    if ~isnan(cap(arc.LogGaussMap,edge.LogGaussMap))
        if arcCenter == 0 && edgeZero == 0
            % This product is NaN, because the arc log-Gauss map is 0, and
            % the edge log-Gauss map is +/-pi/2, so there is no match.
            r = NaN;
        elseif arcCenter == 0
            % The arc log-Gauss map is zero, so the output is the arc 
            % multiplied by the point on the edge with 0 log-Gauss map
            if edge.LogGaussMap.isin(0)
                r = arc .* edge.findLGM(0);
            else
                r = NaN;
            end
        elseif abs(edgeZero) < 100*eps
            % The edge log-Gauss map is +/-pi, so the output is the edge
            % multiplied by the point on the arc with +/-pi log-Gauss map
            if any(arc.LogGaussMap.isin([-pi/2,pi/2]))
                LGM = edge.LogGaussMap;
                if LGM.inf == pi/2
                    r = edge .* arc.findLGM(pi/2);
                elseif LGM.sup == -pi/2
                    r = edge .* arc.findLGM(-pi/2);
                else
                    pArc = [arc.findLGM(-pi/2) ; arc.findLGM(pi/2)];
                    r = edge .* pArc(~isnan(pArc));
                end
            else
                r = NaN;
            end
        else
            % The result is an arc fitted to an ellipse or hyperbole
            % segment
            
            % Extract radii and curve parameters
            R = arc.Radius * abs(arc.NormFactor);
            S = arc.CurveParameter;
            T = edge.CurveParameter;
    
            % Define envelope function and touching condition
            H_st = @(s,t) (1 + R*cos(s) - R*t.*sin(s)) + ...
                      1i* (t + R*t.*cos(s) + R*sin(s));
            j_st = @(s,t) R + cos(s) + sin(s)./t;
    
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
    
            % Fit arc to the envelope
            cCnt = @(u) u * exp(1i*(angle(diff(pHxy))+pi/2)) + mean(pHxy);
            cRad = @(u) abs(cCnt(u)-pHxy(1));
            % if abs(R)<1
            %     cAng0 = angle(diff(pHxy))+pi/2;
            % else
            %     cAng0 = angle(diff(pHxy))-pi/2;
            % end
            cAng1 = @(u) angle(pHxy(1)-cCnt(u));
            cAng2 = @(u) angle(pHxy(2)-cCnt(u));
            cFit = @(u) abs(min(abs(H_smp-cCnt(u))-cRad(u)));
            uSol = fminsearch(cFit,1);
            r = ciat.Arc(cCnt(uSol),cRad(uSol),...
                    ciat.RealInterval(cAng1(uSol),cAng2(uSol)));
            
            % Flip the arc-angle interval if necessary
            % if ~r.ArcAngle.isin(cAng0)
            cAng0 = angle(H_smp(ceil(length(H_smp)/2)) - r.Center);
            if ~r.ArcAngle.isin(cAng0)
                r = ciat.Arc(cCnt(uSol),cRad(uSol),...
                    ciat.RealInterval(cAng2(uSol)-2*pi,cAng1(uSol)));
            end
    
            % Reverse normalization
            r = r * (1 / arc.NormFactor / edge.NormFactor);
            
        end
    else
        r = NaN;
    end
end


%% Multiplication with an arc

function r = timesArc(arc1,arc2)
    arc1Center = arc1.Center;
    arc2Center = arc2.Center;

    if ~isnan(cap(arc1.LogGaussMap,arc2.LogGaussMap))
        if arc1Center == 0 && arc2Center == 0 
            r = ciat.Arc(0,arc1.Radius * arc2.Radius,...
                        arc1.GaussMap + arc2.GaussMap);
        elseif arc1Center == 0 
            % The log-Gauss map of a zero-centered arc is 0, therefore 
            % in this case the product is the zero-centered arc
            % multplied by the point(s) of the non-zero centered arc
            % that has a log-Gauss map of zero
            r = arc1 .* arc2.findLGM(0);
        elseif arc2Center == 0 
            r = arc2 .* arc1.findLGM(0);
        else
            % The result is an arc fitted to a Cartesian oval segment 
    
            % Extract radii and curve parameters
            R1 = arc1.Radius * abs(arc1.NormFactor);
            R2 = arc2.Radius * abs(arc2.NormFactor);
            S = arc1.CurveParameter;
            T = arc2.CurveParameter;
    
            % Define envelope function and touching condition
            H_st = @(s,t) (1 + R1*cos(s) + R2*cos(t) + R1*R2*cos(s + t)) + ...
                  1i*(R1*sin(s) + R2*sin(t) + R1*R2*sin(s + t));
            j_st = @(s,t) sin(s-t)-R1*sin(t)+R2*sin(s);
    
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
    
            % Fit arc to the envelope
            cCnt = @(u) u * exp(1i*(angle(diff(pHxy))+pi/2)) + mean(pHxy);
            cRad = @(u) abs(cCnt(u)-pHxy(1));
            cAng0 = wrapToPi(angle(diff(pHxy))-pi/2);
            cAng1 = @(u) angle(pHxy(1)-cCnt(u));
            cAng2 = @(u) angle(pHxy(2)-cCnt(u));
            cFit = @(u) abs(max(abs(H_smp-cCnt(u))-cRad(u)));
            uSol = fminsearch(cFit,1);
            r = ciat.Arc(cCnt(uSol),cRad(uSol),...
                    ciat.RealInterval(cAng1(uSol),cAng2(uSol)));
    
            % Flip the arc-angle interval if necessary
            if ~r.ArcAngle.isin(cAng0)
                r = ciat.Arc(cCnt(uSol),cRad(uSol),...
                    ciat.RealInterval(cAng2(uSol)-2*pi,cAng1(uSol)));
            end
    
            % Reverse normalization
            r = r * (1 / arc1.NormFactor / arc2.NormFactor);
    
        end    
    end
end

