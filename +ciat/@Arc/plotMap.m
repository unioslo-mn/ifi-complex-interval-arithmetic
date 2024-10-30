function h = plotMap(obj, logMap, arrowSize, varargin)
    tf = ishold;
    if tf == false 
        clf
    end
    hold on
    h = [];

    % Plot normal vectors
    for n = 1:length(obj(:))
        % Extract variables
        arc = obj(n);
        if logMap == 0
            map = arc.GaussMap;
        else
            map = arc.LogGaussMap;
        end

        % Set arrow positions
        p = zeros(1,3);
        p(1) = arc.Startpoint;
        p(2) = arc.Midpoint;
        p(3) = arc.Endpoint;

        % Set vector lengths
        if length(map) == 1
            a = zeros(1,3);
        else
            a = zeros(1,6);
            p = reshape(repmat(p,1,2),1,[]);
        end

        % Set arrow angle and length
        for m = 1:length(map)
            a(1+3*(m-1)) = arrowSize * exp(1i*map(m).Infimum);
            a(2+3*(m-1)) = arrowSize * exp(1i*map(m).Midpoint);
            a(3+3*(m-1)) = arrowSize * exp(1i*map(m).Supremum);
        end

        % Flip vectors if arc is concave 
        if obj(n).Radius < 0
            a = flip(a);
        end

        % Plot arrows
        h=[h; ...
            quiver(real(p),imag(p),real(a),imag(a),...
                   varargin{:},'AutoScale','off')];
    end

    if tf == false 
        hold off
    end
end