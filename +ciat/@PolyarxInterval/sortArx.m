function arx = sortArx(arx)
    if isa(arx,'cell')
        arx = arx{:};
    end
    if arx(end,4) ~= pi
        % [~,idx] = sort(arx(:,4));
        % arx = arx(idx,:);
        [~,idx] = min(arx(:,4));
        arx = circshift(arx,1-idx);
        if arx(end,4) ~= pi
            arx = [ arx ; arx(1,:) ];
            arx(end,4) = pi;
        end
    end
end