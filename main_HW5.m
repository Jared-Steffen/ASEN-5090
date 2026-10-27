clc; clear; close all

%Jared is working on windows and Brady is working on mac
%This allows us to both use the directory.
Windowsflag = ispc;
Macflag = ismac;

if Windowsflag
    addpath("Functions\");
    addpath("Data\");
elseif Macflag
    addpath("Functions/");
    addpath("Data/");
end

set_constants2025();

% Authors: Brady Sivey and Jared Steffen

%% Problem 1

nistECEF = [-1288398.567; -4721696.932; 4078625.350];
referenceDate = datetime(2026,8,19);

obs = rinexread('NIST00USA_R_20262310000_01D_30S_MO.rnx');
gps = obs.GPS;

epochObs = gps(gps.Time == referenceDate, :);

% Extract C1C and C2W
C1C = epochObs.C1C;
C2W = epochObs.C2W;
obsTime = epochObs.Time;
[PRN,sort_idx] = sort([epochObs.SatelliteID]);

% Sort from highest to lowest PRN
C1C = C1C(sort_idx);
C2W = C2W(sort_idx);

% Time of day
time = hours(obsTime - referenceDate);

%[GPS week number, seconds into the week]
gpsTime = [2432*ones(numel(obsTime),1), ...
    259200 + seconds(obsTime - referenceDate)];

%Broadcast ephemeris
eph = read_clean_GPSbroadcast('brdc2310.26n', true);

% Loop through each PRN
for i = 1:length(PRN)

    % PRIF
    [PRIF(i), ~] = ionocorr(C1C(i), L1, C2W(i), L2);

    % Sat pos, vel, clock correction, and relativity correction
    [~, satPosR,~,BSV(:,i),RELAT(:,i)] = ...
        eph2pvt2025(eph, gpsTime(i,:), PRN(i));

    % Expected range
    [R(i),satPosR_tt(:,i)] = ...
        correct_range(gpsTime(i,:),nistECEF,satPosR,eph,PRN(i));

    [AZ(i),EL(i),~] = ecef2az_el_range(satPosR_tt(:,i),nistECEF);

    % Troposphere
    TROP(i) = tropocalc(2, EL(i));

end

% Make and display table for inspection
prob1_table = table(PRN,PRIF',R',EL',AZ',BSV',RELAT',TROP',...
    'VariableNames', {'PRN','PRIF','R','EL','AZ','BSV','RELAT','TROP'});

format longG
disp(prob1_table)

%% Problem 2

% Calculate prefits
prefits = PRIF - (R - BSV - RELAT + TROP);

prefit_table = table(PRN,prefits','VariableNames',{'PRN','Prefits'});
disp(prefit_table)

% Look at measurement sensitivity matrix
H = measSensMat(satPosR_tt,nistECEF,R);

% Solve for receiver position and clock bias using weighted least squares
stateCorrection = (H' * H) \ (H' * prefits(:));

postfits = prefits' - H*stateCorrection;

% Mean and stdv
mean_postfits = mean(postfits);
stdv_postfits = std(postfits);

postfit_table = table(PRN,postfits,'VariableNames',{'PRN','Postfits'});
disp(postfit_table)

% Investigate better solution
removeSatellites = ismember(PRN, [3, 8,27]);
PRN_filtered = PRN(~removeSatellites);
prefits_filtered = prefits(~removeSatellites);
H_filtered = H(~removeSatellites, :);

stateCorrection_filtered = (H_filtered' * H_filtered) \ ...
    (H_filtered' * prefits_filtered(:));

postfits_filtered = prefits_filtered' - H_filtered *stateCorrection_filtered;

postfit_filtered_table = table(PRN_filtered, postfits_filtered, ...
    'VariableNames', {'PRN','Postfits'});
disp(postfit_filtered_table)

% Print RMSE
receiverPositionEstimate = nistECEF + stateCorrection_filtered(1:3);
positionError = norm(receiverPositionEstimate - nistECEF);
disp(positionError)

%% Problem 3
%Note: I used Chat GPT to make the printed tables and values look pretty
%for this problem

%3a
initialLLA3 = [40, -105, 1600]; % Initial guess: [latitude (deg), longitude (deg), altitude (m)]
initialECEF3 = lla2ecef(initialLLA3)';

% Error convention: true position minus estimated position
initialError3 = nistECEF - initialECEF3;

fprintf('\nQUESTION 3(a)\n');
fprintf('Initial ECEF position (m): [%.2f, %.2f, %.2f]\n', initialECEF3);
fprintf('Initial ECEF error (m): [%.2f, %.2f, %.2f]\n', initialError3);
fprintf('Initial distance from NIST: %.2f m\n', norm(initialError3));

% Use all satellites.
excludedPRNs3 = []; %Can be used to exclude PRN

use3 = ~ismember(PRN(:), excludedPRNs3) & isfinite(PRIF(:));
prn3 = PRN(use3);
prn3 = prn3(:);

pif3 = PRIF(use3);
pif3 = pif3(:);

gpsTime3 = gpsTime(use3,:);

%In case we exclude too many satellites
numSat3 = numel(prn3);
assert(numSat3 >= 4, 'At least four satellites are required.');

% State = [x; y; z; receiver clock bias], all in meters]
state3 = [initialECEF3; 0];

tolerance3 = 0.01; %1 cm
maxIterations3 = 30;
history3 = nan(maxIterations3, 7);
converged3 = false;

%3b, 3c, and 3f
for iter3 = 1:maxIterations3

    range3 = zeros(numSat3,1);
    az3 = zeros(numSat3,1);
    el3 = zeros(numSat3,1);
    bsv3 = zeros(numSat3,1);
    relat3 = zeros(numSat3,1);
    trop3 = zeros(numSat3,1);
    satECEF3 = zeros(3,numSat3);

    % Recompute measurements using the current receiver position
    for j3 = 1:numSat3

        % Satellite position at reception time
        [~, satReceive3] = eph2pvt2025(eph, gpsTime3(j3,:), prn3(j3));

        % Correct range for signal travel time and Earth rotation.
        % Flag = true returns the rotated satellite position.
        [range3(j3), satECEF3(:,j3)] = correct_range(gpsTime3(j3,:), state3(1:3), satReceive3, eph, prn3(j3), true);

        % Satellite clock and relativity at transmission time
        transmitTime3 = gpsTime3(j3,:);
        transmitTime3(2) = transmitTime3(2) - range3(j3)/c;

        [~, ~, ~, bsv3(j3), relat3(j3)] = eph2pvt2025(eph, transmitTime3, prn3(j3));

        % Azimuth, elevation, and tropospheric delay
        [az3(j3), el3(j3)] = ecef2az_el_range(satECEF3(:,j3), state3(1:3));

        trop3(j3) = tropocalc(2, el3(j3));
    end

    % Predicted pseudorange, including the current receiver clock bias
    predicted3 = range3 + state3(4) - bsv3 - relat3 + trop3;

    % Prefit residuals
    prefit3 = pif3 - predicted3;

    % Recompute geometry matrix.
    % Flag = true uses the receiver-position derivative signs.
    G3 = measSensMat(satECEF3, state3(1:3), range3, true);

    % Least squares correction [dx; dy; dz; db]
    correction3 = G3 \ prefit3;

    % Linearized post-fit residuals
    postfit3 = prefit3 - G3*correction3;

    if iter3 == 1
        fprintf('\nQUESTION 3(b): FIRST ITERATION\n');

        firstTable3 = table(prn3, pif3, range3, el3, az3, bsv3, relat3, trop3, ...
            'VariableNames', ...
            {'PRN','PIF','R','EL','AZ','BSV','RELAT','TROP'});
        disp(firstTable3)

        disp(table(prn3, prefit3, postfit3, ...
            'VariableNames', {'PRN','Prefit_m','Postfit_m'}))

        disp('First geometry matrix G:')
        disp(G3)

        disp('First correction [dx; dy; dz; db] (m):')
        disp(correction3)
    end

    % Update position and receiver clock bias
    state3 = state3 + correction3;

    correctionMagnitude3 = norm(correction3);
    distanceError3 = norm(nistECEF - state3(1:3));

    history3(iter3,:) = [iter3, correctionMagnitude3, state3', distanceError3];

    %Stopping point
    if correctionMagnitude3 < tolerance3
        converged3 = true;
        break
    end
end

if ~converged3
    error('Question 3 did not converge within %d iterations.',  maxIterations3);
end

% Display each
history3 = history3(1:iter3,:);

fprintf('\nQUESTION 3(c,f): ITERATION HISTORY\n');
fprintf(['Iter   ||correction|| (m)       X (m)           Y (m)' ...
         '           Z (m)         b (m)    Position error (m)\n']);

for k3 = 1:iter3
    fprintf('%3d       %12.2f   %14.2f  %14.2f  %14.2f  %12.2f  %12.2f\n', ...
        history3(k3,:));
end

%3g
% Final ECEF error: true minus estimated
finalErrorECEF3 = nistECEF - state3(1:3);

fprintf('\nQUESTION 3(g): FINAL SOLUTION\n');
fprintf('Final ECEF position (m): [%.2f, %.2f, %.2f]\n', state3(1:3));
fprintf('Receiver clock bias: %.2f m = %.9e s\n', state3(4), state3(4)/c);
fprintf('Final ECEF error (m): [%.2f, %.2f, %.2f]\n',  finalErrorECEF3);
fprintf('Final position error magnitude: %.2f m\n', norm(finalErrorECEF3));
fprintf('Final correction magnitude: %.6f m\n', correctionMagnitude3);

disp(table(prn3, postfit3, 'VariableNames', {'PRN','FinalPostfit_m'}))

%3h
%ENU to NIST location
nistLLA3 = ecef2lla(nistECEF');
C3 = ECEF2ENU(nistLLA3(1), nistLLA3(2));

finalErrorENU3 = C3*finalErrorECEF3;

fprintf('\nQUESTION 3(h): ENU ERROR\n');
fprintf('[East; North; Up] error (m): [%.2f, %.2f, %.2f]\n', finalErrorENU3);

%3i
% Recompute satellite geometry at the final position for DOP
finalRange3 = zeros(numSat3,1);
finalSatECEF3 = zeros(3,numSat3);

for j3 = 1:numSat3

    [~, satReceive3] = eph2pvt2025(eph, gpsTime3(j3,:), prn3(j3));

    [finalRange3(j3), finalSatECEF3(:,j3)] = correct_range(gpsTime3(j3,:), state3(1:3), satReceive3, eph, prn3(j3), true);
end

finalG3 = measSensMat(finalSatECEF3, state3(1:3), finalRange3, true);

% Geometry covariance for unit measurement variance
Q3 = (finalG3'*finalG3) \ eye(4);

% Rotate the position block from ECEF to ENU
Qenu3 = C3*Q3(1:3,1:3)*C3';

HDOP3 = sqrt(Qenu3(1,1) + Qenu3(2,2));
VDOP3 = sqrt(Qenu3(3,3));

fprintf('\nQUESTION 3(i): DOP\n');
fprintf('HDOP = %.4f\n', HDOP3);
fprintf('VDOP = %.4f\n', VDOP3);

%% Problem 4
%Note: I used Chat GPT to make the printed tables and values look pretty
%for this problem

% Starting guess much farther from NIST
initialLLA4 = [0, -105, 0];
initialECEF4 = lla2ecef(initialLLA4)';

% State: [x; y; z; receiver clock bias], all in meters
state4 = [initialECEF4; 0];

fprintf('\nQUESTION 4\n');
fprintf('Initial distance from NIST: %.2f m\n', norm(nistECEF - initialECEF4));

fprintf('Iter   Correction magnitude (m)   Position error (m)\n');

for iter4 = 1:maxIterations3

    range4 = zeros(numSat3,1);
    bsv4 = zeros(numSat3,1);
    relat4 = zeros(numSat3,1);
    trop4 = zeros(numSat3,1);
    satECEF4 = zeros(3,numSat3);

    for j4 = 1:numSat3

        % Satellite position at reception time
        [~, satReceive4] = eph2pvt2025(eph, gpsTime3(j4,:), prn3(j4));

        % Range and Earth-rotation-adjusted satellite position
        [range4(j4), satECEF4(:,j4)] = correct_range(gpsTime3(j4,:), state4(1:3), satReceive4, eph, prn3(j4), true);

        % Satellite clock and relativity at transmission time
        transmitTime4 = gpsTime3(j4,:);
        transmitTime4(2) = transmitTime4(2) - range4(j4)/c;

        [~, ~, ~, bsv4(j4), relat4(j4)] = eph2pvt2025(eph, transmitTime4, prn3(j4));

        % Elevation and tropospheric delay at the current estimate
        [~, el4] = ecef2az_el_range(satECEF4(:,j4), state4(1:3));

        trop4(j4) = tropocalc(2, el4);
    end

    % Prefit residuals and updated geometry matrix
    predicted4 = range4 + state4(4) - bsv4 - relat4 + trop4;
    prefit4 = pif3 - predicted4;

    G4 = measSensMat(satECEF4, state4(1:3), range4, true);

    % Solve and update
    correction4 = G4 \ prefit4;
    state4 = state4 + correction4;

    fprintf('%3d          %14.2f          %14.2f\n', iter4, norm(correction4), norm(nistECEF - state4(1:3)));

    if norm(correction4) < tolerance3
        break
    end
end

if norm(correction4) >= tolerance3
    error('Question 4 did not converge within the iteration limit.');
end

finalErrorECEF4 = nistECEF - state4(1:3);

fprintf('\nFinal ECEF position (m): [%.2f, %.2f, %.2f]\n', ...
    state4(1:3));
fprintf('Final ECEF error (m): [%.2f, %.2f, %.2f]\n', ...
    finalErrorECEF4);
fprintf('Final position error magnitude: %.2f m\n', ...
    norm(finalErrorECEF4));
fprintf('Final clock bias: %.2f m = %.9e s\n', ...
    state4(4), state4(4)/c);
fprintf('Final correction magnitude: %.6f m\n', norm(correction4));

% Compare with Question 3
fprintf('\nPosition difference from Question 3: %.6f m\n', ...
    norm(state4(1:3) - state3(1:3)));
fprintf('Clock bias difference from Question 3: %.6f m\n', ...
    abs(state4(4) - state3(4)));