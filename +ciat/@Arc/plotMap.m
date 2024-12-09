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
        arc = obj(n);
        if logMap == 0
            map = arc.GaussMap;
        else
            map = arc.LogGaussMap;
        end

        % Set arrow positions and length
        p = arc.sample(arrowCount);
        a = arrowSize * exp(1i*map.sample(arrowCount).');
        
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