function r = sum(obj,varargin)
% Sum of polyarcular intervals
%
% This function creates the polyarcular interval representing the 
% sum of a set of polyarcular intervals
% _________________________________________________________________________
% USAGE        
%   r = sum(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.PolyarcularInterval class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   polyInt = sum([ciat.PolyarcularInterval([0,1,1i]), ...
%                    ciat.PolyarcularInterval([0,-1,-1i])]);
% _________________________________________________________________________

    [M,N] = size(obj);
    allConvex = all(obj.isconvex,'all');

    % Check for special case
    if allConvex && M == 1
        r = obj.quickSum;
        return
    end

    if ( M == 1 && N == 1)
        r = obj;
    elseif M == 1
        r = sum2(obj,allConvex);
    elseif N == 1
        r = sum1(obj,allConvex);
    else
        if size(varargin) == 0
            r = sum1(obj,allConvex);
        else
            if varargin{1} == 1
                r = sum1(obj,allConvex);
            elseif varargin{1} == 2
                r = sum2(obj,allConvex);
            elseif strcmp(varargin{1},'all')
                r = sumAll(obj,allConvex);
            else
                error('Parameter two is invalid.')
            end
        end
    end    
end

%% Function for summing along dimension 1
function r = sum1(obj,allConvex)
    [M,N] = size(obj);
    r = obj(1,:);
    for n = 1:N
        for m = 2:M
            % r(n) = r(n) + obj(m, n);
            r(n) = addPolyarc(r(n),obj(m, n),allConvex);
        end
    end
end

%% Function for summing along dimension 2
function r = sum2(obj,allConvex)
    [M,N] = size(obj);
    r = obj(:,1);
    for m = 1:M
        for n = 2:N
            % r(m) = r(m) + obj(m, n);
            r(m) = addPolyarc(r(m),obj(m, n),allConvex);
        end
    end
end

%% Function for summing along all dimension
function r = sumAll(obj,allConvex)
    [M,N] = size(obj);
    r = 0;
    for m = 1:M
        for n = 1:N
            % r = r + obj(m, n);
            r = addPolyarc(r,obj(m, n),allConvex);
        end
    end
end

%% Function for adding two polyarcs
function r = addPolyarc(a,b,allConvex)
    if allConvex
        r = ciat.PolyarcularInterval.plusConvex(a,b);
    else
        r = ciat.PolyarcularInterval.plusConcave(a,b);
    end
end