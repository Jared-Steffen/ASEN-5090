function [AZ, EL, R] = satellitevisparams(userECEF, satECEF)
%userECEF is 3x1
%satECEF is 3xn

%converted to lat lon
userLLA = ecef2lla(userECEF');

userLat = userLLA(1);
userLon = userLLA(2);

%function from 2b
C_ECEF2ENU = ECEF2ENU(userLat, userLon);


%line of sight vector
LOSvecECEF = satECEF - userECEF;

%covertint to ENU
LOSvecENU = C_ECEF2ENU * LOSvecECEF;
east = LOSvecENU(1,:);
north = LOSvecENU(2,:);
up = LOSvecENU(3,:);

%geometrical range (m)
R = sqrt(east.^2 + north.^2 + up.^2);

%Azimuth angle (deg), must be positive
AZ = mod(atan2d(east, north), 360);

%Elevation angle (deg)
EL = asind(up ./ R);

end

function C_ECEF2ENU = ECEF2ENU(ref_lat_deg, ref_lon_deg)

%convert to rad
ref_lat = deg2rad(ref_lat_deg);
ref_lon = deg2rad(ref_lon_deg);

%Transformation Matrix (I took this from lecture 4)
C_ECEF2ENU =[-sin(ref_lon) cos(ref_lon) 0; 
    -sin(ref_lat)*cos(ref_lon), -sin(ref_lat)*sin(ref_lon), cos(ref_lat);
    cos(ref_lat)*cos(ref_lon), cos(ref_lat)*sin(ref_lon), sin(ref_lat)];

end