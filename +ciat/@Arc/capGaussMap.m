function output = capGaussMap(input1, input2)
    
    % Extract sizes
    [M1,N1] = size(input1);
    [M2,N2] = size(input2);

    % Split angles
    input1 = splitAngle(input1);
    input2 = splitAngle(input2);
    L1 = size(input1,3) > 1;
    L2 = size(input2,3) > 1;

    % If the inputs are two arrays of different orientation
    % form matrices
    if (M1 == M2) && (N1 == N2)
        M = size(input1,1);
        N = size(input1,2);
    elseif (N1 == 1) && (M2 == 1)
        input1 = repmat(input1,1,N2);
        input2 = repmat(input2,M1,1);
        M = M1;
        N = N2;
    elseif (M1 == 1) && (N2 == 1)
        input1 = repmat(input1,M2,1);
        input2 = repmat(input2,1,N1);
        M = M2;
        N = N1;
    else
        error('Incorrect input size for operation.')
    end
  
    % Intersect angles
    primaryCap = cap(input1(:,:,1),input2(:,:,1));
    if ~L1 && L2
        secondaryCap = cap(input1(:,:,1),input2(:,:,2));
    elseif L1 && ~L2
        secondaryCap = cap(input1(:,:,2),input2(:,:,1));
    elseif L1 && L2
        % Merge the adjacent intersections 
        infPositive = primaryCap.Infimum > 0;
        tempCap = cap(input1(:,:,2),input2(:,:,2));
        primaryCap = union( primaryCap - infPositive*2*pi, ...
                            tempCap + ~infPositive*2*pi);

        % There can be only one secondary intersection
        tempCap = cat( 3 , cap(input1(:,:,1),input2(:,:,2)) , ...
                           cap(input1(:,:,2),input2(:,:,1)) );
        secondaryCap(M,N) = ciat.RealInterval;
        secondaryCap.Infimum = max(tempCap.Infimum,[],3);
        secondaryCap.Supremum = max(tempCap.Supremum,[],3);

    end

    % Join primary and secondary caps
    if ~L1 && ~L2
        output = primaryCap;
    else
        moveMask = isnan(primaryCap) & ~isnan(secondaryCap);
        if any(moveMask,'all')
            emptyInterval(1) = ciat.RealInterval;
            primaryCap(moveMask) = secondaryCap(moveMask);
            secondaryCap(moveMask) = emptyInterval;
        end
        
        if all(isnan(secondaryCap),'all')
            output = primaryCap;
        else
            output = cat(3 , primaryCap , secondaryCap);
            lastwarn('Surplus angle intersections found!')
        end
    end
end


%% Function for splitting the angles including the -pi or pi value

function output = splitAngle(input)

    % Check if any elements contain the -pi or pi value
    mask = input.Infimum < -pi | input.Supremum > pi;
    
    % Create a second layer of the array in the 3rd dimension
    if any(mask,'all')
        [M,N] = size(input);
        input2(M,N) = ciat.RealInterval;
        angInf = ciat.wrapToPi(input.Infimum(mask));
        angSup = ciat.wrapToPi(input.Supremum(mask));
        input2(mask) = ciat.RealInterval(-pi*ones(size(angSup)),angSup);
        input(mask) = ciat.RealInterval( angInf , pi*ones(size(angInf)));
        output = cat(3,input,input2);
    else
        output = ciat.RealInterval(...
                            ciat.wrapToPi(input.Infimum) , ...
                            ciat.wrapToPi(input.Supremum) );
    end

end
        