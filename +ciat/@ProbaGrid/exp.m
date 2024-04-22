function obj = exp(obj)
    % % Compute the distribution of the exp of the random variable.

    % First compute the new grid
    exp_points = exp(obj.x + 1i * obj.y.');
    new_x = linspace(min(real(exp_points(:))), max(real(exp_points(:))), obj.nx);
    new_y = linspace(min(imag(exp_points(:))), max(imag(exp_points(:))), obj.ny);
    [new_X, new_Y] = meshgrid(new_x, new_y);

    % Note that this would need to be summed for each multiple of 2*pi for the actual result
    % However, for the case of the product (which is the only one we need), the result only
    % needs three terms, as the other are zero.
    pdf = (interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)),         imag(log(new_X + 1i * new_Y)), "linear", 0) + ...
           interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)+2*1j*pi), imag(log(new_X + 1i * new_Y)+2*1j*pi), "linear", 0) + ...
           interp2(obj.x, obj.y, obj.Pdf, real(log(new_X + 1i * new_Y)-2*1j*pi), imag(log(new_X + 1i * new_Y)-2*1j*pi), "linear", 0)) ...
            ./ abs(new_X + 1i * new_Y);

    % Create a new ProbaGrid object
    obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
    % obj = obj.normalize();
end