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

%% Problem 4