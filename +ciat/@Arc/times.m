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
    r(M,N) = ciat.Arc;
    for m = 1:M
        for n = 1:N
            % One of the operands is a double
            if isa(obj1,"ciat.Arc") && isa(obj2,"double")
                r = timesDouble(obj1,obj2);
            elseif isa(obj1,"double") && isa(obj2,"ciat.Arc")
                r = timesDouble(obj2,obj1);

            % One of the operands is an edge
            elseif isa(obj1,"ciat.Arc") && isa(obj2,"ciat.Edge")
                r = timesEdge(obj1,obj2);
            elseif isa(obj1,"ciat.Edge") && isa(obj2,"ciat.Arc")
                r = timesEdge(obj2,obj1);

            % Both operands are arcs
            elseif isa(obj1,"ciat.Arc") && isa(obj2,"ciat.Arc")
                r = timesArc(obj1,obj2);
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

    

    if ~isnan(logGaussMap)
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
        elseif edgeZero==0
            % The edge log-Gauss map is +/-pi, so the output is the edge
            % multiplied by the point on the arc with +/-pi log-Gauss map
            if any(arc.LogGaussMap.isin([-pi/2,pi/2]))
                LGM = edge.LogGaussMap;
                if LGM.inf == pi/2
                    r = edge .* arc.findLGM(pi/2);
                elseif LGM.sup == -pi/2
                    r = edge .* arc.findLGM(-pi/2);
                else
                    r = edge .* [arc.findLGM(-pi/2) ; ...
                                 arc.findLGM(pi/2)];
                end
            else
                r = NaN;
            end
        else
            % The result is an ellipse or hyperbole section between the 
            % log-Gauss map intersection bounds
            warning('This function is not yet implemented')
            % Get log-Gauss-map intersection
            LGM = ciat.Arc.capGaussMap(arc.LogGaussMap,edge.LogGaussMap);
            
            
            
        end
    else
        r = NaN;
    end
end


%% Multiplication with an arc

function r = timesArc(arc1,arc2)
    warning('This function is not yet implemented')
   
    arc1Center = arc1.Center;
    arc2Center = arc2.Center;

    if arc1Center == 0 && arc2Center == 0 
        r = ciat.Arc(0,arc1.Radius + arc2.Radius,...
                    cap(arc.GaussMap,arc2.GaussMap));
    elseif arc1Center == 0 


    elseif arc2Center == 0 

    else

    end        
end

