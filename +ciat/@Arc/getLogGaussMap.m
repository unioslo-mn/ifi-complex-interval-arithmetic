function value = getLogGaussMap(obj)
    % Extract properties
    [M,N] = size(obj);
    isZeroCentered = obj.Center == 0;
    isVertex = obj.Radius == 0;
    isConvex = obj.Radius >= 0;
    isInZero = obj.isin(0);
    isCurveZero = obj.CurveParameter.isin(0);
    isCurvePi = obj.CurveParameter.isin(-pi) | obj.CurveParameter.isin(pi);
   
    % Calculate log-Gauss map values
    gFunc = @(s,R) atan2(R.*sin(s),R.*cos(s))-atan2(R.*sin(s),1+R.*cos(s));
    LGMinf = gFunc(obj.CurveParameter.inf,obj.Radius .* abs(obj.NormFactor));
    LGMsup = gFunc(obj.CurveParameter.sup,obj.Radius .* abs(obj.NormFactor));

    % Assign log-Gauss map value
    value(M,N) = ciat.RealInterval;
        % For zero-centered arcs the LGM is zero
    mask = isZeroCentered;
    if any(mask,'all')
        value(mask) = ciat.RealInterval(zeros(sum(mask,'all'),1));
    end
        % For zero-radius arcs the LGM is calculated from the Gauss map
    mask = isVertex;
    if any(mask,'all')
        LGMinf = ciat.wrapToPi(obj.GaussMap.inf - angle(obj.Center) );
        LGMsup = ciat.wrapToPi(obj.GaussMap.sup - angle(obj.Center) );
        value(mask) = ciat.RealInterval(LGMinf,LGMsup);
    end
        % For non-zero-centered arcs that does not include
        % the origin the LGM is simply the LGM of the endpoints
    mask = ~isZeroCentered & ~isInZero & isConvex;
    if any(mask,'all')
        value(mask) = ciat.RealInterval(LGMinf(mask)-2*pi*isCurvePi(mask), ...
                                        LGMsup(mask));
    end
    mask = ~isZeroCentered & ~isInZero & ~isConvex;
    if any(mask,'all')
        value(mask) = ciat.RealInterval(LGMinf(mask)-2*pi*isCurveZero(mask), ...
                                        LGMsup(mask));
    end
        % For non-zero centered arcs that contain the origin
        % the LGM has an envelope
    mask = ~isZeroCentered & isInZero;
    if any(mask,'all')
          value(mask) = ciat.RealInterval(LGMinf(mask),LGMsup(mask));

        % Check each arc if the LGM includes the function envelope
        list = 1:M*N;
        list = list(mask);
        for idx = 1:length(list)
            % Find envelope
            [sEnv,gEnv] = findLogGaussMapEnvelope(obj(list(idx)));

            % If the parameter interval contains the envelope
            sInt = obj(list(idx)).CurveParameter;
            if any(sInt.isin([sEnv,sEnv+2*pi]))
                value(list(idx)).Supremum = gEnv;
            end
            if any(sInt.isin([-sEnv,-sEnv+2*pi]))
                value(list(idx)).Infimum = -gEnv;
            end 
        end
    end
end

%% Utility functions
function [sEnv,gEnv] = findLogGaussMapEnvelope(obj)
    % Find envelope
    syms s           
    norm = obj.NormFactor;
    rad = obj.Radius * abs(norm);
    g(s) = angle(rad*exp(1i*s)) - angle(1+rad*exp(1i*s));
    dg(s) = (real(exp(s*1i)).^2* ...
            (imag(exp(s*1i)).^2/real(exp(s*1i)).^2 + 1))/...
            (imag(exp(s*1i)).^2 + real(exp(s*1i)).^2) - ...
            (((rad.^2*imag(exp(s*1i)).^2) / ...
            (rad*real(exp(s*1i)) + 1).^2 + ...
            (rad*real(exp(s*1i)))/(rad*real(exp(s*1i)) + 1))* ...
            (rad*real(exp(s*1i)) + 1).^2)/...
            (rad.^2*imag(exp(s*1i)).^2 + ...
            (rad*real(exp(s*1i)) + 1).^2);
    sEnv = abs(ciat.wrapToPi(eval(vpasolve(abs(dg)==0,s))));
    gEnv = eval(g(sEnv));
        
end