function h = plotMap(obj,logMap,arrowSize,varargin)
    tf = ishold;
    if tf == false 
        clf
    end
    hold on
    h = [];

    % Plot normal vectors
    for n = 1:length(obj(:))
        % Extact variables
        edge = obj(n);
        if logMap == 0
            map = edge.GaussMap;
        else
            map = edge.LogGaussMap;
        end

        % Set arrow positions
        p(1) = edge.Startpoint;
        p(2) = edge.Midpoint;
        p(3) = edge.Endpoint;

        % Set vector lengths
        if length(map) == 1
            a = zeros(1,3);
        else
            a = zeros(1,6);
            p = repmat(p,1,2);
        end

        % Set arrow angle and length
        for m = 1:length(map)
            a(1+3*(m-1)) = arrowSize * exp(1i*map(m).Infimum);
            a(2+3*(m-1)) = arrowSize * exp(1i*map(m).Midpoint);
            a(3+3*(m-1)) = arrowSize * exp(1i*map(m).Supremum);
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