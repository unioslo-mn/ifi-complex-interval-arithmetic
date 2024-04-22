function obj = log(obj)
    % Compute the distribution of the log of the random variable.

    % First compute the new grid
    log_points = log(obj.x + 1i * obj.y.');
    new_x = linspace(min(real(log_points(:))), max(real(log_points(:))), obj.nx);
    new_y = linspace(min(imag(log_points(:))), max(imag(log_points(:))), obj.ny);

    exp_points = exp(new_x + 1i * new_y.');

    pdf = interp2(obj.x, obj.y, obj.Pdf, real(exp_points), imag(exp_points), "linear", 0) .* abs(exp_points);

    % Create a new ProbaGrid object
    obj = ciat.ProbaGrid.from_pdf(pdf, new_x, new_y);
    % obj = obj.normalize();
end