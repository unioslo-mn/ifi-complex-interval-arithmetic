function r = timesDouble(xObj,d)

    arx = xObj.Arx{:};
    cen = complex(arx(:,1) , arx(:,2)) * d;
    rad = arx(:,3) * abs(d);
    ang = arx(:,4) + angle(d);
        
    r = ciat.PolyarxInterval([real(cen),imag(cen),rad,ang]);
end