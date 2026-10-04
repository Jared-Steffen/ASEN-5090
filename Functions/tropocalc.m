function tropo = tropocalc(zd, elevation)
% zd: zenith delay (m)
% elevation: satellite elevation (degrees)
% tropo: line-of-sight tropospheric delay (m)

tropo = zd ./ sqrt(1-(cosd(elevation)/1.001).^2);

end