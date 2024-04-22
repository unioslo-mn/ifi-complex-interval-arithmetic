function r = plus(obj1, obj2)

% Check input class
    mustBeA(obj1,"ciat.ProbaGrid");
    mustBeA(obj2,"ciat.ProbaGrid");
    
    % Get input sizes and check if they can be added
    [M1,N1] = size(obj1);
    [M2,N2] = size(obj2);
    assert(M1 == M2 || M1 == 1 || M2 == 1)
    assert(N1 == N2 || N1 == 1 || N2 == 1)
    M = max([M1,M2]);
    N = max([N1,N2]);
    
    % Loop throught the arrays
    r(M,N) = ciat.ProbaGrid;
    for m=1:M
        for n=1:N
            % Calculate indexes
            m1 = min(m,M1);
            n1 = min(n,N1);
            m2 = min(m,M2);
            n2 = min(n,N2);
            
            % Calculate sum
            r(m,n) = add( obj1(m1,n1) , obj2(m2,n2) );
        end
    end
end

%% Function for adding two probability grids
function obj = add(obj1, obj2)
    % Compute the distribution of the sum of the two random
    % variables. This is done by convolving the two pdfs.

    new_nx = max(obj1.nx, obj2.nx);
    new_ny = max(obj1.ny, obj2.ny);

    % Interpolate the grid with the smallest step size
    [pdf1, pdf2, new_dx, new_dy] = match_step_size(obj1, obj2);

    % Convolve the two pdfs
    pdf = conv2(pdf1, pdf2, "full") .* new_dx .* new_dy;

    % Create a new ProbaGrid object
    new_x = linspace(obj1.x(1) + obj2.x(1), obj1.x(end) + obj2.x(end), length(pdf(1,:)));
    new_y = linspace(obj1.y(1) + obj2.y(1), obj1.y(end) + obj2.y(end), length(pdf(:,1)));
    obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
    % obj = obj.normalize();  

    % Keep the same grid size
    new_x = linspace(obj.x(1), obj.x(end), new_nx);
    new_y = linspace(obj.y(1), obj.y(end), new_ny);
    obj = obj.adjust(new_x, new_y);
end

%% Match step size
% Function that takes two ProbaGrid objects and returns their two pdfs,
% where the one with the smallest step size is interpolated on the grid of
% the obj2
function [pdf1, pdf2, new_dx, new_dy] = match_step_size(obj1, obj2)
    new_dx = max(obj1.dx, obj2.dx);
    new_dy = max(obj1.dy, obj2.dy);

    pdf1 = interp_pdf(obj1, new_dx, new_dy);
    pdf2 = interp_pdf(obj2, new_dx, new_dy);
end
