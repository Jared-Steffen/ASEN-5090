function [az,el,rho] = ecef2az_el_range(ecef_sat,ecef_receiver)
% Takes in an ECEF position of 2 locations and then computes az, el, and
% range

% Get relative position vector
relative_pos_ecef = ecef_sat - ecef_receiver;

% Range calculation -- frame doesn't matter since it is a scalar
rho = vecnorm(relative_pos_ecef);

% Convert receiver ECEF to LLA
lla_receiver = ecef2lla(ecef_receiver(:)');
C_ECEF2ENU = ECEF2ENU(lla_receiver(1),lla_receiver(2));
relative_pos_enu = C_ECEF2ENU*relative_pos_ecef;

% Get Az/El
az = atan2(relative_pos_enu(1,:),relative_pos_enu(2,:));
el = asin(relative_pos_enu(3,:)./rho);

% Return in degrees
az = rad2deg(wrapTo2Pi(az));
el = rad2deg(el);

end