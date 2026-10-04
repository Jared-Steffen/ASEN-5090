function corrected_rho = correct_range(t,receiver_ecef,sat_ecef0,broadcast_data,prn)
% Takes in a time vector that is [wn, tow] (outputs from cal2gps), the
% receiver position in ecef (3x1), initial sat position in ecef (3xN), 
% broadcast_data (the output from read_clean_GPSbroadcast) and the PRN
% number of interest

% Set constants
c = 299792458; % Speed of light (meters/s).
wE  = 7.2921151467e-5; % WGS-84 value, rad/s 

% Loop through time
corrected_rho = zeros(1,size(t,1));
for i = 1:length(t)

    % Calculate initial range
    tr = t(i,2);
    wn = t(i,1);
    sat_pos_tr = sat_ecef0(:,i);
    rho = norm(sat_pos_tr-receiver_ecef);
    
    % Correct geometric range for signal travel time
    e = 1;
    % Loop until converged
    while e > 1e-6
    
        % Transmission time
        tt = tr - rho/c;
    
        % Get ECEF position at tt
        [~,sat_pos_tt,~,~,~,~] = eph2pvt2025(broadcast_data,[wn, tt],prn);
        rho2 = norm(sat_pos_tt'-receiver_ecef);
    
        % Rotate to ECEF at tr
        phi = wE*(tr-tt);
        C_ECEFTT2ECEFTR = [cos(phi), sin(phi), 0;
                           -sin(phi), cos(phi) 0;
                           0, 0, 1];
    
        sat_pos_tr = C_ECEFTT2ECEFTR * sat_pos_tt';
    
        % Evalute new range
        new_rho = norm(sat_pos_tr-receiver_ecef);
    
        % Check convergence
        e = abs(new_rho-rho);
        rho = new_rho;
    end
    
    corrected_rho(i) = rho;
end

end