function r = times(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.Edge","double"]);
    mustBeA(obj2,["ciat.Edge","double"]);
    
    % Get input sizes and check if they can be combined
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);

    % Loop throught the arrays
    r(M,N) = ciat.Edge;
    for m = 1:M
        for n = 1:N
            % One of the operands is a double
            if isa(obj1,"ciat.Edge") && isa(obj2,"double")
                r = timesDouble(obj1(m,n),obj2(m,n));
            elseif isa(obj1,"double") && isa(obj2,"ciat.Edge")
                r = timesDouble(obj2(m,n),obj1(m,n));
            else
                error('invalid input')
            end
        end
    end

    
end

%% Multiplication with a double

function r = timesDouble(edge,x)
    r = ciat.Edge(edge.Startpoint * x,...
                  edge.Endpoint * x );
end

