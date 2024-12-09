function h = plotMap(obj, logMap, arrowSize, varargin)


    % Extract arrowcount variable
    arrowCount = 3;
    if ~isempty(varargin)
        idx = 1;
        while idx <= length(varargin)
            if strcmp(varargin{idx},'arrowCount')
                if length(varargin) > idx
                    arrowCount = varargin{idx+1};
                    varargin = varargin(setdiff(1:length(varargin),[idx idx+1]));
                else
                    varargin = varargin(1:end-1);
                end
            end
            idx = idx+1;
        end
    end


    tf = ishold;
    if tf == false 
        clf
    end
    hold on
    h = [];

    % Plot normal vectors
    for n = 1:length(obj(:))
        % Extract variables
        edge = obj(n);
        if logMap == 0
            map = edge.GaussMap;
        else
            map = edge.LogGaussMap;
        end

        % Set arrow positions and length
        p = edge.sample(arrowCount);
        a = arrowSize * exp(1i*map.sample(arrowCount).');

        % Plot arrows
        h=[h; ...
            quiver(real(p),imag(p),real(a),imag(a),...
                   varargin{:},'AutoScale','off')];
    end

    if tf == false 
        hold off
    end
end