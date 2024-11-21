function r = sum(obj,varargin)
% Sum of arcs
%
% This function creates the edge representing the 
% sum of a set of edge intervals
% _________________________________________________________________________
% USAGE        
%   r = sum(obj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   obj       : array of objects from the ciat.Edge class
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
% _________________________________________________________________________

    [M,N] = size(obj);
    if ( M == 1 && N == 1)
        r = obj;
    elseif ( M == 1 || N == 1 )
        r = 0;
        for n = 1:max(M,N)
            r = r + obj(n);
        end
    else
        if size(varargin) == 0
            r = sum1(obj);
        else
            if varargin{1} == 1
                r = sum1(obj);
            elseif varargin{1} == 2
                r = sum2(obj);
            elseif strcmp(varargin{1},'all')
                r = sumAll(obj);
            else
                error('Parameter two is invalid.')
            end
        end
    end    
end

%% Function for summing along dimension 1
function r = sum1(obj)
    [M,N] = size(obj);
    r = obj(1,:);
    for n = 1:N
        for m = 2:M
            r(n) = r(n) + obj(m, n);
        end
    end
end

%% Function for summing along dimension 2
function r = sum2(obj)
    [M,N] = size(obj);
    r = obj(:,1);
    for m = 1:M
        for n = 2:N
            r(m) = r(m) + obj(m, n);
        end
    end
end

%% Function for summing along all dimension
function r = sumAll(obj)
    [M,N] = size(obj);
    r = 0;
    for m = 1:M
        for n = 1:N
            r = r + obj(m, n);
        end
    end
end