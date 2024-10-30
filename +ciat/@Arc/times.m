function r = times(obj1,obj2)
   
    % Check input class
    mustBeA(obj1,["ciat.Arc","double"]);
    mustBeA(obj2,["ciat.Arc","double"]);
    
    % Get input sizes and check if they can be combined
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    if isa(obj1,"ciat.Arc") && isa(obj2,"double")
        r = ciat.Arc(obj1.Center .* obj2 , ...
                     obj1.Radius .* abs(obj2) , ...
                     obj1.ArcAngle + angle(obj2));
    elseif isa(obj1,"double") && isa(obj2,"ciat.Arc")
        r = ciat.Arc(obj2.Center .* obj1 , ...
                     obj2.Radius .* abs(obj1) , ...
                     obj2.ArcAngle + angle(obj1));
    else
        error('Invalid input, only arc times double allowed.')
    end
end
