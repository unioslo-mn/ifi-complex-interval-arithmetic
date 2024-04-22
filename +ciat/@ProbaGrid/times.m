function obj = times(obj, other)
    % Compute the distribution of the product of the two random
    % variables.

    new_nx = max(obj.nx, other.nx);
    new_ny = max(obj.ny, other.ny);

    obj = exp(log(obj) + log(other));
    % obj = obj.normalize();

    % Keep the same grid size
    new_x = linspace(obj.x(1), obj.x(end), new_nx);
    new_y = linspace(obj.y(1), obj.y(end), new_ny);
    obj = obj.adjust(new_x, new_y);
end