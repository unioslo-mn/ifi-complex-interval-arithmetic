function outObj = cast(inObj)

% Cast complex intervals of other types to circular interval type
%
% This function takes one or more complex intervals of another type
% and creates the smallest inclusive interval(s) of the circular
% interval type.
% _________________________________________________________________________
% USAGE        
%   outObj = ciat.CircularInterval.cast(inObj)
% _________________________________________________________________________
% NECESSARY ARGUMENT
%   inObj       : object of one of the following types:
%                   - ciat.RectangularInterval
%                   - ciat.PolarInterval
%                   - ciat.PolygonalInterval               
% _________________________________________________________________________
% OPTIONS
% _________________________________________________________________________
% EXAMPLES
%   circInt = ciat.CircularInterval(ciat.RectangularInterval(1,3,2,4));
% _________________________________________________________________________
%
% Copyright (C) 2023 H. Arnestad and G. Gereb, BSD-3
% If you use this software, please cite it as in CITATION.cff
% Project: Beampattern Interval Analysis 
% Website: doi.org/10.5281/zenodo.6856232
% Contact: haavaarn@uio.no, gaborge@uio.no
% (More information in README.md and LICENSE.md.)
% _________________________________________________________________________


    [M,N] = size(inObj);

    switch class(inObj)
        case 'double'
            outCenter = inObj;
            outRadius = zeros(M,N);
        case 'ciat.RectangularInterval'
            inReal = inObj.Real;
            inImag = inObj.Imag;
            outCenter = complex( [inReal.Midpoint] , [inImag.Midpoint] );
            outRadius = sqrt( ([inReal.Width]/2).^2 + ([inImag.Width]/2).^2 );
            
        case 'ciat.PolarInterval'
            % Extract parameters
            inAbs = [inObj.Abs];
            inAngle = [inObj.Angle];
            angDiff = [inAngle.Supremum] - [inAngle.Infimum];
            angAvg = ([inAngle.Supremum] + [inAngle.Infimum])/2;
            rU = inAbs.Supremum;
            rL = inAbs.Infimum;
            outCenterOpt1 = rU .* cos(angDiff/2);
            outCenterOpt2 = (rU + rL) ./ (2*cos(angDiff/2));
            outer_corner = rU * exp( 1j* (angAvg + angDiff/2) );
            
            % Calculate outCenter and outRadius
            outCenter = min(outCenterOpt1, outCenterOpt2);
            outCenter(angDiff >= pi) = 0;
            outCenter = outCenter .* exp( 1j * angAvg );
            outRadius = abs(outer_corner - outCenter);
            
        case 'ciat.PolygonalInterval'
            outCenter = zeros(M,N);
            outRadius = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    inReal = real(inObj(m,n).Points);
                    inImag = imag(inObj(m,n).Points);
                    % [c,R] = minboundcircle(inReal,inImag);
                    [R,c,~] = ExactMinBoundCircle([inReal, inImag]);
                    Rcorr = max(abs(complex(inReal-c(1),inImag-c(2))) - R);
                    outCenter(m,n) = c(1) + 1i*c(2);
                    outRadius(m,n) = R + Rcorr;
                end
            end
        case 'ciat.PolyarcularInterval' % Temporary solution
            outCenter = zeros(M,N);
            outRadius = zeros(M,N);
            
            for m = 1:M
                for n = 1:N
                    points = inObj(m,n).sample(10);
                    inReal = real(points);
                    inImag = imag(points);
                    % [c,R] = minboundcircle(inReal,inImag);
                    [R,c,~] = ExactMinBoundCircle([inReal(:),inImag(:)]);
                    outCenter(m,n) = c(1) + 1i*c(2);
                    outRadius(m,n) = R+10*eps;
                end
            end
        case 'ciat.PolyarxInterval' % Temporary solution
            outCenter = zeros(M,N);
            outRadius = zeros(M,N);

            for m = 1:M
                for n = 1:N
                    points = inObj(m,n).sample(10);
                    inReal = real(points);
                    inImag = imag(points);
                    % [c,R] = minboundcircle(inReal,inImag);
                    [R,c,~] = ExactMinBoundCircle([inReal,inImag]);
                    outCenter(m,n) = c(1) + 1i*c(2);
                    outRadius(m,n) = R+10*eps;
                end
            end
         case 'ciat.Arc' % Temporary solution
            outCenter = zeros(M,N);
            outRadius = zeros(M,N);
            for m = 1:M
                for n = 1:N
                    points = inObj(m,n).sample(10);
                    inReal = real(points);
                    inImag = imag(points);
                    [c,R] = minboundcircle(inReal,inImag);
                    outCenter(m,n) = c(1) + 1i*c(2);
                    outRadius(m,n) = R;
                end
            end
        otherwise
            error('Invalid input type')
    end
	outObj = ciat.CircularInterval(outCenter,outRadius);       
    outObj = reshape(outObj,M,N);
end

%% Utility functions

function [center,radius] = minboundcircle(x,y,hullflag)
% minboundcircle: Compute the minimum radius enclosing circle of a set of (x,y) pairs
% usage: [center,radius] = minboundcircle(x,y,hullflag)
%
% arguments: (input)
%  x,y - vectors of points, describing points in the plane as
%        (x,y) pairs. x and y must be the same size. If x and y
%        are arrays, they will be unrolled to vectors.
%
%  hullflag - boolean flag - allows the user to disable the
%        call to convhulln. This will allow older releases of
%        matlab to use this code, with a possible time penalty.
%        It also allows minboundellipse to call this code
%        efficiently.
% 
%        hullflag = false --> do not use the convex hull
%        hullflag = true  --> use the convex hull for speed
%
%        default: true
%
%
% arguments: (output)
%  center - 1x2 vector, contains the (x,y) coordinates of the
%        center of the minimum radius enclosing circle
%
%  radius - scalar - denotes the radius of the minimum
%        enclosing circle
%
%
% Example usage:
%   x = randn(50000,1);
%   y = randn(50000,1);
%   tic,[c,r] = minboundcircle(x,y);toc
%
%   Elapsed time is 0.171178 seconds.
%
%   c: [-0.2223 0.070526]
%   r: 4.6358
%
%
% See also: minboundrect
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 1/10/07
% default for hullflag
if (nargin<3) || isempty(hullflag)
  hullflag = true;
elseif ~islogical(hullflag) && ~ismember(hullflag,[0 1])
  error 'hullflag must be true or false if provided'
end
% preprocess data
x=x(:);
y=y(:);
% not many error checks to worry about
n = length(x);
if n~=length(y)
  error 'x and y must be the same sizes'
end
% start out with the convex hull of the points to
% reduce the problem dramatically. Note that any
% points in the interior of the convex hull are
% never needed.
if hullflag && (n>3)
  edges = convhulln([x,y]);
  % list of the unique points on the convex hull itself
  % convhulln returns them as edges
  edges = unique(edges(:));
  % exclude those points inside the hull as not relevant
  x = x(edges);
  y = y(edges);
    
end
% now we must find the enclosing circle of those that
% remain.
n = length(x);
% special case small numbers of points. If we trip any
% of these cases, then we are done, so return.
switch n
  case 0
    % empty begets empty
    center = [];
    radius = [];
    return
  case 1
    % with one point, the center has radius zero
    center = [x,y];
    radius = 0;
    return
  case 2
    % only two points. center is at the midpoint
    center = [mean(x),mean(y)];
    radius = norm([x(1),y(1)] - center);
    return
  case 3
    % exactly 3 points
    [center,radius] = enc3(x,y);
    return
end
% more than 3 points.
% Use an active set strategy.
aset = 1:3; % arbitrary, but quite adequate
iset = 4:n;
% pick a tolerance
tol = 10*eps*(max(abs(mean(x) - x)) + max(abs(mean(y) - y)));
% Keep a list of old sets as tried to trap any cycles. we don't need to
% retain a huge list of sets, but only a few of the best ones. Any cycle
% must hit one of these sets. Really, I could have used a smaller list,
% but this is a small enough size that who cares? Almost always we will
% never even fill up this list anyway.
old.sets = NaN(10,3);
old.rads = inf(10,1);
old.centers = NaN(10,2);
flag = true;
while flag
  % have we seen this set before? If so, then we have entered a cycle
  aset = sort(aset);
  if ismember(aset,old.sets,'rows')
    % we have seen it before, so trap out
    center = old.centers(1,:);
    radius = old.radius(1);
    
    % just reset flag then continue, and the while loop will terminate
    flag = false;
    continue
  end
  
  % get the enclosing circle for the current set
  [center,radius] = enc3(x(aset),y(aset));
  
  % is this better than something from the retained sets?
  if radius < old.rads(end)
    old.sets(end,:) = sort(aset);
    old.rads(end) = radius;
    old.centers(end,:) = center;
        
    % sort them in increasing order of the circle radii
    [old.rads,tags] = sort(old.rads,'ascend');
    old.sets = old.sets(tags,:);
    old.centers = old.centers(tags,:);
  end
  
  % are all the inactive set points inside the circle?
  r = sqrt((x(iset) - center(1)).^2 + (y(iset) - center(2)).^2);
  [rmax,k] = max(r);
  if (rmax - radius) <= tol
    % the active set enclosing circle also enclosed
    % all of the inactive points.
    flag = false;
  else
    % it must be true that we can replace one member of aset
    % with iset(k). Which one?
    s1 = [aset([2 3]),iset(k)];
    [c1,r1] = enc3(x(s1),y(s1));
    if (norm(c1 - [x(aset(1)),y(aset(1))]) <= r1)
      center = c1;
      radius = r1;
      
      % update the active/inactive sets
      swap = aset(1);
      aset = [iset(k),aset([2 3])];
      iset(k) = swap;
      
      % bounce out to the while loop
      continue
    end
    s1 = [aset([1 3]),iset(k)];
    [c1,r1] = enc3(x(s1),y(s1));
    if (norm(c1 - [x(aset(2)),y(aset(2))]) <= r1)
      center = c1;
      radius = r1;
      
      % update the active/inactive sets
      swap = aset(2);
      aset = [iset(k),aset([1 3])];
      iset(k) = swap;
      
      % bounce out to the while loop
      continue
    end
    s1 = [aset([1 2]),iset(k)];
    [c1,r1] = enc3(x(s1),y(s1));
    if (norm(c1 - [x(aset(3)),y(aset(3))]) <= r1)
      center = c1;
      radius = r1;
      
      % update the active/inactive sets
      swap = aset(3);
      aset = [iset(k),aset([1 2])];
      iset(k) = swap;
      
      % bounce out to the while loop
      continue
    end
    
    % if we get through to this point, then something went wrong.
    % Active set problem. Increase tol, then try again.
    tol = 2*tol;
    
  end
  
end
end
% =======================================
%  begin subfunction
% =======================================
function [center,radius] = enc3(X,Y)
% minimum radius enclosing circle for exactly 3 points
%
% x, y are 3x1 vectors
% convert to complex
xy = X + sqrt(-1)*Y;
% just in case the points are collinear or nearly so, get
% the interpoint distances, and test the farthest pair
% to see if they work.
Dij = @(XY,i,j) abs(XY(i) - XY(j));
D12 = Dij(xy,1,2);
D13 = Dij(xy,1,3);
D23 = Dij(xy,2,3);
% Find the most distant pair. Test if their circumcircle
% also encloses the third point.
if (D12>=D13) && (D12>=D23)
  center = (xy(1) + xy(2))/2;
  radius = D12/2;
  if abs(center - xy(3)) <= radius
    center = [real(center),imag(center)];
    return
  end
elseif (D13>=D12) && (D13>=D23)
  center = (xy(1) + xy(3))/2;
  radius = D13/2;
  if abs(center - xy(2)) <= radius
    center = [real(center),imag(center)];
    return
  end
elseif (D23>=D12) && (D23>=D13)
  center = (xy(2) + xy(3))/2;
  radius = D23/2;
  if abs(center - xy(1)) <= radius
    center = [real(center),imag(center)];
    return
  end
end
% if we drop down to here, then the points cannot
% be collinear, so the resulting 2x2 linear system
% of equations will not be singular.
A = 2*[X(2)-X(1), Y(2)-Y(1); X(3)-X(1), Y(3)-Y(1)];
rhs = [X(2)^2 - X(1)^2 + Y(2)^2 - Y(1)^2; ...
       X(3)^2 - X(1)^2 + Y(3)^2 - Y(1)^2];
     
center = (A\rhs)';
radius = norm(center - [X(1),Y(1)]);
end

%% Alternative smallest circle algorithm

function [R,C,Xb] = ExactMinBoundCircle(X)
% Compute exact minimum bounding circle of a 2D point cloud using 
% Welzl's algorithm [1]. 
%
% INPUT
%   - X     : M-by-2 list of point coordinates, where M is the total
%             number of points.
%
% OUTPUT
%   - R     : radius of the minimum bounding circle of X.
%   - C     : 1-by-2 vector specifying centroid co-ordinates of the minimum
%             bounding circle of X.
%   - Xb    : subset of X, listing K-by-2 list of point co-ordinates from 
%             which R and C were computed. See function
%             'FitCircle2Points' for more info.
%
% REREFERENCES:
% [1] Welzl, E. (1991), 'Smallest enclosing disks (balls and ellipsoids)',
%     Lecture Notes in Computer Science, Vol. 555, pp. 359-370
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%
if nargin<1 || isempty(X) || ~isnumeric(X) || ~ismatrix(X) 
    error('Unrecognized format for 1st input argument (X)')
elseif size(X,2)~=2
    error('This function only works for 2D data')
elseif any(~isfinite(X(:)))
    error('Point data contains NaN and/or Inf entries. Remove them and try again.')    
end
% Get minimum bounding circle
% -------------------------------------------------------------------------
if size(X,1)>2
    
    % Check points for collinearity
    Xo = mean(X,1);
    dX = bsxfun(@minus,X,Xo);
    [U,D] = svd((dX'*dX)/size(X,1),0);
    D = diag(D);
    if D(2)<1E-15
        dx = dX*U(:,1);
        [dx_min,id_min] = min(dx);
        [dx_max,id_max] = max(dx);
        R = (dx_max - dx_min)/2;
        C = U(:,1)'*(dx_max + dx_min)/2 + Xo;
        Xb = X([id_min; id_max],:);
        return
    end
    
    % Get convex hull of the point set
    F = convhull(X);
    F = unique(F(:));
    X = X(F,:);
    
end
try
    % Remove duplicates
    X = uniquetol(X,eps,'ByRows',true);
catch
    % older version of Matlab; 'uniquetol' is unavailable
end
if size(X,1)<3
    [R,C] = FitCircle2Points(X); 
    Xb = X;   
    return
end
% Randomly permute the point set
idx = randperm(size(X,1));
X = X(idx(:),:);
if size(X,1)<1E3
    try
        
        % Center and radius of the circle
        [R,C] = B_MinCircle(X,[]);
        
        % Co-ordinates of the points used to compute parameters of the 
        % minimum bounding circle
        D = sum(bsxfun(@minus,X,C).^2,2);
        [D,idx] = sort(abs(D-R^2));
        Xb = X(idx(1:4),:);
        D = D(1:4);
        Xb = Xb(D<1E-6,:);
        [~,idx] = sort(Xb(:,1));
        Xb = Xb(idx,:);
        return
    catch
    end
end
    
% If we got to this point, then either size(X,1)>=1E3 or recursion depth 
% limit was reached. So need to break-up point-set into smaller sets and 
% then recombine the results.
M = size(X,1);
dM = max(min(floor(M/4),300),3);
res = mod(M,dM);
n = ceil(M/dM);  
idx = dM*ones(1,n);
if res>0
    idx(end) = res;
end
 
if res<0.25*dM && res>0
    idx(n-1) = idx(n-1)+idx(n);
    idx(n) = [];
    n = n - 1;
end
X = mat2cell(X,idx,2);
Xb = [];
for i = 1:n
    
    % Center and radius of the circle
    [R,C,Xi] = B_MinCircle([Xb;X{i}],[]);    
    
    % 40 points closest to the circle
    if i<1
        D = abs(sum(bsxfun(@minus,Xi,C).^2,2)-R^2);
    else
        D = abs(sqrt(sum(bsxfun(@minus,Xi,C).^2,2))-R);
    end
    [D,idx] = sort(D);
    Xb = Xi(idx(1:min(40,numel(D))),:);
    
end
D = D(1:3);
Xb = Xb(D/R*100<1E-3,:);
[~,idx] = sort(Xb(:,1));
Xb = Xb(idx,:);
    function [R,C,P] = B_MinCircle(P,B)
        
    if size(B,1)==3 || isempty(P)
        [R,C] = FitCircle2Points(B); % fit circle to boundary points
        return
    end
    
    % Remove the last (i.e., end) point, p, from the list
    P_new = P;
    P_new(end,:) = [];
    p = P(end,:);
        
    % Check if p is on or inside the bounding circle. If not, it must be
    % part of the new boundary.
    [R,C,P_new] = B_MinCircle(P_new,B); 
    if isnan(R) || isinf(R) || R<=eps
        chk = true;
    else
        chk = norm(p-C)>(R+eps);
    end
    
    if chk
        B = [p; B];
        [R,C] = B_MinCircle(P_new,B);
        P = [p; P_new];
    end
        
    end
end

function [R,C]=FitCircle2Points(X)
    % Fit a circle to a set of 2 or at most 3 points in 3D space. Note that
    % point configurations with 3 collinear points do not have well-defined 
    % solutions (i.e., they lie on circles with infinite radius).
    %
    % INPUT:
    %   - X     : M-by-2 array of point coordinates, where M<=3.
    %
    % OUTPUT:
    %   - R     : radius of the circle. R=Inf when the circle is undefined, as 
    %             specified above.
    %   - C     : coordinates of the circle centroid. C=nan(1,2) when the 
    %             circle is undefined, as specified above.
    %
    % AUTHOR: Anton Semechko (a.semechko@gmail.com)
    %
    N=size(X,1);
    if N>3
        error('Input must a N-by-2 array of point coordinates, with N<=3')
    end
    % Empty set
    if isempty(X)
        C=nan(1,2);
        R=nan; 
        return
    end
    % A single point
    if N==1
        C=X;
        R=0;
        return
    end
    % Line segment
    if N==2
        C=mean(X,1);
        R=norm(X(2,:)-X(1,:))/2;
        return
    end
    % Remove duplicate vertices, if there are any
    D=bsxfun(@minus,permute(X,[1 3 2]),permute(X,[3 1 2]));
    D=sqrt(sum(D.^2,3));
    D(1:(N+1):end)=Inf;
    chk=D<=1E-12;
    if sum(chk(:))>0
        for i=1:(N-1)
            if size(X,1)<=i, break; end
            idx=chk(i,:);
            idx(1:i)=false;
            idx=find(idx);
            chk(idx,:)=[];
            chk(:,idx)=[];
            X(idx,:)=[];
        end
        [R,C]=FitCircle2Points(X);
        return
    end
    % Three unique, though possibly collinear points
    tol=1E-2; % collinearity threshold (in degrees)
       
    % Check for collinearity
    D12=X(2,:)-X(1,:); D12=D12/norm(D12);
    D13=X(3,:)-X(1,:); D13=D13/norm(D13);
    chk=abs(D12*D13(:));
    chk(chk>1)=1;
    if acos(chk)*(180/pi)<tol
        R=inf;
        C=nan(1,2);
        return
    end
    % Circle centroid
    A=2*bsxfun(@minus,X(2:3,:),X(1,:));
    b=sum(bsxfun(@minus,X(2:3,:).^2,X(1,:).^2),2);
    C=(A\b)';
    % Circle radius
    R=sqrt(sum((X(1,:)-C).^2,2));

end