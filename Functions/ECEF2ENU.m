function C_ECEF2ENU = ECEF2ENU(lat,lon)
% Takes a latitude and longitude in degrees and forms the DCM between ECEF
% and ENU frames

lat = deg2rad(lat);
lon = deg2rad(lon);

C_ECEF2ENU = [-sin(lon), cos(lon), 0;
    -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat);
    cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];

end