function tropo = tropocalc(zd, elevation)
% zd: zenith delay (m)
% elevation: satellite elevation (degrees)
% tropo: line-of-sight tropospheric delay (m)

tropo = zd ./ sind(elevation);

end