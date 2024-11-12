function r = isin(obj,x)
    
    [Mo,No] = size(obj);
    [Mx,Nx] = size(x);
    
    if Mo*No == 1
        r = inPolyarc(obj,x);
    elseif Mx*Nx == 1
        r(Mo,No) = false;
        for m = 1:Mo
            for n = 1:No
                r(m,n) = inPolyarc(obj(m,n), x);
            end
        end
    elseif No*Mx == 1
        r(Mo,Nx) = false;
        for m = 1:Mo
            r(m,:) = inPolyarc(obj(m), x);
        end
    elseif Mo*Nx == 1
        r(Mx,No) = false;
        for n = 1:No
            r(:,n) = inPolyarc(obj(n),x.').';
        end
    else
        error('Arrays have incompatible sizes for this operation.')
    end
end

%% Function for checking if points are inside a single polyarc

function r = inPolyarc(obj,x)

    % Make sure x is a horizontal vector
    [M,N] = size(x);
    x = x(:).';

    % Check first if the values are inside the inclusive rectangle
    r = ciat.RectangularInterval(obj).isin(x);

    % If there are, then check it more closely
    if any(r)
        % Extract arcs
        arcs = [obj.Arcs ; obj.Vertices];
    
        % Check if the point is inside the polygon
        vertexPoly = [arcs.Startpoint , arcs.Endpoint].';
        vertexPoly = [unique(vertexPoly(:),'stable') ; arcs.Startpoint(1)];
        inVertexPoly = inpolygon(real(x),imag(x),real(vertexPoly),imag(vertexPoly));
    
        % Check if the point is inside the arcs
        convexArcs = arcs(arcs.Radius>0);
        concaveArcs = arcs(arcs.Radius<0);
        inConvexArcs = any(convexArcs.isin(x),'all');
        inConcaveArcs = any(concaveArcs.isin(x),'all');
    
        % Combine conditions
        r = (inVertexPoly & ~inConcaveArcs) | inConvexArcs;
        r = reshape(r,M,N);
    end
end



