function plotAreas(obj, interval, n_areas, varargin)
    % Separate the pdf into n_areas areas, of the same shape
    % as the original interval but smaller and centered at the
    % baricenter of the original interval
    % For each area, plot the probability to be on that zone
    % as a heatmap
    % Thus the sum of the probabilities of all the areas is 1,
    % and the graphs contains n_areas colored zones.

    % Compute the center of the interval
    center = interval.Real.Midpoint + 1i * interval.Imag.Midpoint;

    interval_centered = interval - center;

    [X,Y] = meshgrid(obj.x, obj.y);

    probas = zeros(n_areas, 1);
    probas_grid = zeros(size(obj.Pdf));
    old_mask = zeros(size(obj.Pdf));

    t_values = linspace(0, 1, n_areas+1);
    hold on
    for i = 1:n_areas
        % Compute the new interval
        new_interval = interval_centered * t_values(i+1) + center;

        % Compute the mask of the new interval
        mask = new_interval.ininterval(X + 1i * Y);
        mask = mask & ~old_mask;

        % Save the mask for the next iteration to avoid overlapping
        old_mask = mask | old_mask;

        % Compute the probability of being in the new interval
        proba_area = sum(obj.Pdf(mask), "all") * obj.dx * obj.dy;

        probas(i) = proba_area;

        % Fill the probability grid with the computed probability
        probas_grid(mask) = proba_area;
    end

    % Plot the heatmap
    imagesc(obj.x, obj.y, probas_grid, 'AlphaData', interval.ininterval(X + 1i * Y), varargin{:})
    set(gca,'YDir','normal')
    % colormap turbo
    axis equal
    colorbar

    % Plot all the intervals
    for i = 1:n_areas
        % Compute the new interval
        new_interval = interval_centered * t_values(i+1) + center;

        % Plot the new interval
        % new_interval.plot("Color", [0 0 0], "LineWidth", 2, 'HandleVisibility','off');
        
        % Plot the interval in red if it's the first with cumulative probability greater than 0.5
        if sum(probas(1:i)) > 0.5 && sum(probas(1:i-1)) < 0.5
            new_interval.plot("Color", [1 0 0], "LineWidth", 4, 'DisplayName', "50% probability");
        end

        % Same in blue for 0.9
        if sum(probas(1:i)) > 0.9 && sum(probas(1:i-1)) < 0.9
            new_interval.plot("Color", [0 0 1], "LineWidth", 4, 'DisplayName', "90% probability");
        end

        % Same in green for 0.99
        if sum(probas(1:i)) > 0.99 && sum(probas(1:i-1)) < 0.99
            new_interval.plot("Color", [0 1 0], "LineWidth", 4, 'DisplayName', "99% probability");
        end 
    end
end