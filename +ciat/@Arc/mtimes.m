function r = mtimes(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.Arc","ciat.Edge","double"]);
    mustBeA(obj2,["ciat.Arc","ciat.Edge","double"]);
    
    % Get input sizes and check if they can be combined
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(N1 == M2 || (M1 == 1 && N1 == 1) || (M2 == 1 && N2 == 1) )
    if N1 == M2
        M = M1;
        N = N2;
    else
       M = max(M1,M2);
       N = max(N1,N2);
    end
    
    % Loop through the arrays
    r(M,N) = ciat.Arc;
    if N1 == M2 && N1>1
        for m = 1:M
            for n = 1:N
                % Calculate product using the times function
                r(m,n) = sum( obj1(m,:) .* obj2(:,n).' );
            end
        end
    else
        % Calculate product using the times function
        r = obj1 .* obj2;
    end
end