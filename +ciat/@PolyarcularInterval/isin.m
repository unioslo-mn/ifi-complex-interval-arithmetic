function r = isin(obj,x)
    
    [Mo,No] = size(obj);
    [Mx,Nx] = size(x);
    
    if Mo*No == 1
        r = inPolyarc(obj.ArcStorage{:},x);
    elseif Mx*Nx == 1
        r(Mo,No) = false;
        for m = 1:Mo
            for n = 1:No
                r(m,n) = inPolyarc(obj.ArcStorage{m,n},x);
            end
        end
    elseif No*Mx == 1
        r(Mo,Nx) = false;
        for m = 1:Mo
            r(m,:) = inPolyarc(obj.ArcStorage{m},x);
        end
    elseif Mo*Nx == 1
        r(Mx,No) = false;
        for n = 1:No
            r(:,n) = inPolyarc(obj.ArcStorage{n},x.').';
        end
    else
        error('Arrays have incompatible sizes for this operation.')
    end
end

%% Function for checking if points are inside a single polyarc

function r = inPolyarc(arcs,x)
    % Make sure x is a horizontal vector
    [M,N] = size(x);
    x = x(:).';

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