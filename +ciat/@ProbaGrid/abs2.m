function obj = abs2(obj)
    % Compute the distribution of the squared magnitude of the
    % random variable.
    % For a variable Z = X + iY, this is |Z|^2 = X^2 + Y^2
    % Thus this is a real random variable.

    % The new probability grid is only 1D as opposed to the other
    % operations

    new_nx = max(obj.nx, obj.ny);

    % First compute the new grid
    abs2_points = abs(obj.x + 1i * obj.y.').^2;
    new_x = linspace(min(abs2_points(:)), max(abs2_points(:)), new_nx);
    new_y = 0;

    % Compute the marginal distributions
    f_X = real(obj).Pdf;
    f_Y = imag(obj).Pdf;

    pdf = zeros(new_nx, 1);
    for i = 1:length(new_x)
        z = new_x(i);
        if z == 0
            continue;
        end
        % Compute the integral
        % Define functions to integrate
        f_Xf = @(u) interp1(obj.x, f_X, u, 'linear', 0);
        f_Yf = @(u) interp1(obj.y, f_Y, u, 'linear', 0);

        f = @(u) ( (f_Xf(u) + f_Xf(-u)) .* ...
                        (f_Yf(sqrt(z-u.^2)) + f_Yf(-sqrt(z-u.^2))) ...
                    + (f_Yf(u) + f_Yf(-u)) .* ...
                        (f_Xf(sqrt(z-u.^2)) + f_Xf(-sqrt(z-u.^2))) ...
                 ) ./ sqrt(z-u.^2);

        pdf(i) = integral(f, 0, sqrt(z/2))/2;
    end

    % Create a new ProbaGrid object
    obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);

    % Normalize the PDF
    obj = obj.normalize();
end