function r = mtimes(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.PolyarxInterval","double"]);
    mustBeA(obj2,["ciat.PolyarxInterval","double"]);
    
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
    
    
    % Loop throught the arrays
    r(M,N) = ciat.PolyarxInterval;
    for m = 1:M
        for n = 1:N
           if isa(obj1, 'ciat.PolyarxInterval') && isa(obj2,'double')
                r(M,N) = timesDouble(obj1(m,n) , obj2(m,n));
            elseif isa(obj1, 'double') && isa(obj2,'ciat.PolyarxInterval')
                r(M,N) = timesDouble(obj2(m,n) , obj1(m,n));
            end
        end
    end
end
