function H = measSensMat(satellite_ecef,receiver_ecef,R)
%{
satellite_ecef : 3xN matrix of satellite positions
receiver_ecef : 3x1 vector of prior estimate for receiver location
R : 1xN matrix of ranges 

N is the number of satellites at a time t
%}

% Extract components
satx = satellite_ecef(1,:);
saty = satellite_ecef(2,:);
satz = satellite_ecef(3,:);

rx = receiver_ecef(1);
ry = receiver_ecef(2);
rz = receiver_ecef(3);

% Compute line-of-sight sensitivity components
H = ones(length(R),4);
for i = 1:length(R)
    H(i,1) = (satx(i) - rx) / R(i);
    H(i,2) = (saty(i) - ry) / R(i);
    H(i,3) = (satz(i) - rz) / R(i);
end

end