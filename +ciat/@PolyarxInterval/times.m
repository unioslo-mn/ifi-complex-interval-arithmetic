function r = times(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.PolyarxInterval","ciat.PolarInterval","double"]);
    mustBeA(obj2,["ciat.PolyarxInterval","ciat.PolarInterval","double"]);
    
    % Get input sizes and check if they can be combined
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    
    % Loop throught the arrays
    r(M,N) = ciat.PolyarxInterval;
    for m = 1:M
        for n = 1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);

            if isa(obj1, 'ciat.PolyarxInterval') && isa(obj2,'double')
                r(M,N) = timesDouble(obj1(m1,n1) , obj2(m2,n2));
            elseif isa(obj1, 'double') && isa(obj2,'ciat.PolyarxInterval')
                r(M,N) = timesDouble(obj2(m2,n2) , obj1(m1,n1));
            elseif isa(obj1, 'ciat.PolyarxInterval') && ...
                   isa(obj2,'ciat.PolarInterval')
                r(M,N) = timesPolar(obj1(m1,n1) , obj2(m2,n2));
            end
        end
    end
end

