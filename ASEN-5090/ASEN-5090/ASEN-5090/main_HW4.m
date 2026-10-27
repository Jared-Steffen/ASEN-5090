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
prn = 14;
nistECEF = [-1288398.567; -4721696.932; 4078625.350];
referenceDate = datetime(2026,8,19);

obs = rinexread('NIST00USA_R_20262310000_01D_30S_MO.rnx');
gps = obs.GPS;

prnObs = gps(gps.SatelliteID == prn, :);
prnObs = prnObs(2:end, :); %dropping the first epoch

obsTime = prnObs.Properties.RowTimes;
time = hours(obsTime - referenceDate);
C1C = prnObs.C1C;

%[GPS week number, seconds into the week]
gpsTime = [2432*ones(numel(obsTime),1), 259200 + seconds(obsTime - referenceDate)];

%Broadcast ephemeris
eph = read_clean_GPSbroadcast('brdc2310.26n', true);

[~, satPosR] = eph2pvt2025(eph, gpsTime, prn);

[~, ~, R] = satellitevisparams(nistECEF, satPosR.');
R = R.';

tol = 1e-3; %m
delta = Inf;
iteration = 0;

while delta > tol
    iteration = iteration + 1;

    travelTime = R/c;
    txGpsTime = [gpsTime(:,1), gpsTime(:,2) - travelTime];

    %Satellite position at transmission time
    [~, satPosT] = eph2pvt2025(eph, txGpsTime, prn);

    %Rotate satellite position into ECEF
    phi = wE*travelTime;

    xrot =  cos(phi).*satPosT(:,1) + sin(phi).*satPosT(:,2);
    yrot = -sin(phi).*satPosT(:,1) + cos(phi).*satPosT(:,2);
    zrot = satPosT(:,3);

    satPosRot = [xrot, yrot, zrot];

    previousR = R;

    [~, el, R] = satellitevisparams(nistECEF, satPosRot.');
    R = R.';

    delta = max(abs(R - previousR));

    fprintf('Iteration %d: maximum change = %.6f m\n', iteration, delta);
end

%final transmission time
txGpsTime = [gpsTime(:,1), gpsTime(:,2) - R/c];

%Measured minus expected geometric range
dPR0 = C1C - R;

%Plot uncorrected Pseudorange Residual
figure;
plot(time, dPR0, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('dPR_0 = C1C - R (m)');
title('PRN 14 Uncorrected Pseudorange Residual');

% First and last values
fprintf('\nProblem 1:\n');
fprintf('First dPR0: %.3f m at %.6f hours\n',dPR0(1), time(1));
fprintf('Last dPR0: %.3f m at %.6f hours\n', dPR0(end), time(end));

%% Problem 2

%Satellite clock correction at transmission time
[~, ~, ~, bsv] = eph2pvt2025(eph, txGpsTime, prn);

%Plot satellite clock correction
figure;
plot(time, bsv, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('Satellite clock correction b_{sv} (m)');
title('PRN 14 Satellite Clock Correction');

dPR1 = C1C - (R - bsv);

%Plot corrected residual
figure;
plot(time, dPR1, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('dPR_1 = C1C - (R - b_{sv}) (m)');
title('PRN 14 Residual After Satellite Clock Correction');

%first and last values
fprintf('\nProblem 2:\n');
fprintf('First dPR1: %.3f m at %.6f hours\n', dPR1(1), time(1));
fprintf('Last dPR1: %.3f m at %.6f hours\n', dPR1(end), time(end));

%% Problem 3

%edited ephemeris function
[~, ~, ~, bsv, relsv] = eph2pvt2025(eph, txGpsTime, prn);

%Plot relativistic correction
figure;
plot(time, relsv, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('Relativistic correction rel_{sv} (m)');
title('PRN 14 Relativistic Clock Correction');

%Apply both clock and relativistic corrections
dPR2 = C1C - (R - bsv - relsv);

%Plot corrected residual
figure;
plot(time, dPR2, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('dPR_2 = C1C - (R - b_{sv} - rel_{sv}) (m)');
title('PRN 14 Residual After Clock and Relativistic Corrections');

%Print first and last values
fprintf('\nProblem 3:\n');
fprintf('First dPR2: %.3f m at %.6f hours\n', dPR2(1), time(1));
fprintf('Last dPR2: %.3f m at %.6f hours\n', dPR2(end), time(end));

%% Problem 4
%Elevation
[~, el, ~] = satellitevisparams(nistECEF, satPosRot.');
el = el.';

%tropospheric delay model
zd = 2;
tropo = tropocalc(zd, el);

%Plot tropospheric delay
figure;
plot(time, tropo, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('Tropospheric delay (m)');
title('PRN 14 Simple Tropospheric Delay');

dPR3 = C1C - (R - bsv - relsv + tropo);

%Plot corrected residual
figure;
plot(time, dPR3, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('dPR_3 = C1C - (R - b_{sv} - rel_{sv} + tropo) (m)');
title('PRN 14 Residual After Clock, Relativity, and Tropo Corrections');

%first and last values
fprintf('\nProblem 4:\n');
fprintf('First dPR3: %.3f m at %.6f hours\n', dPR3(1), time(1));
fprintf('Last dPR3: %.3f m at %.6f hours\n', dPR3(end), time(end));

%% Problem 5
% Extract
C2L = prnObs.C2L;
L1C = prnObs.L1C;
L2L = prnObs.L2L;

% Ionosphere free model
[pseudorange_IF, iono_corr] = ionocorr(C1C, L1, C2L, L2);

% Plot ionosphere 
figure;
plot(time, iono_corr, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('Ionospheric correction (m)');
title('PRN 14 Ionospheric Correction');

dPR4 = pseudorange_IF - (R - bsv - relsv + tropo);

% Plot corrected residual
figure;
plot(time, dPR4, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('dPR_4 = PRIF - (R - b_{sv} - rel_{sv} + tropo) (m)');
title('PRN 14 Residual After Clock, Relativity, Tropo, and Iono Corrections');

%first and last values
fprintf('\nProblem 4:\n');
fprintf('First dPR4: %.3f m at %.6f hours\n', dPR4(1), time(1));
fprintf('Last dPR4: %.3f m at %.6f hours\n', dPR4(end), time(end));

% Carrier phase
[carrierphase_IF, iono_corr] = ionocorr_CP(L1C, L1, L2L, L2);

% Plot ionosphere 
figure;
plot(time, iono_corr, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('Ionospheric correction (m)');
title('PRN 14 Ionospheric Correction w/ Carrier Phase');

dCP4alt = carrierphase_IF - (R - bsv - relsv + tropo);

% Plot corrected residual
figure;
plot(time, dCP4alt, 'b.');
grid on;
xlabel('Hours into August 19, 2026');
ylabel('dCP4alt = CPIF - (R - b_{sv} - rel_{sv} + tropo) (m)');
title('PRN 14 Residual After Clock, Relativity, Tropo, and Iono Corrections');

%first and last values
fprintf('\nProblem 4:\n');
fprintf('First dCP4alt: %.3f m at %.6f hours\n', dCP4alt(1), time(1));
fprintf('Last dCP4alt: %.3f m at %.6f hours\n', dCP4alt(end-1), time(end));

%% Problem 6

% Plot everything on the same plot
figure();
plot(time, dPR1, 'b.')
hold on
plot(time, dPR2, 'r.')
plot(time, dPR3, 'g.')
plot(time, dPR4, 'k.')
plot(time, dCP4alt, 'm.')
grid on
xlabel('Hours into August 19, 2026')
ylabel('Pseudorange residual (m)')
legend('dPR_1','dPR_2','dPR_3','dPR_4','dCP4alt')
title('PRN 14 Pseudorange Residual Corrections')

%% Problem 7

% Extract needed measurements
C2W = prnObs.C2W;
C5Q = prnObs.C5Q;
L2W = prnObs.L2W;
L5Q = prnObs.L5Q;

% Extract SNRs
S1C = prnObs.S1C;
S2W = prnObs.S2W;
S2L = prnObs.S2L;
S5Q = prnObs.S5Q;

% Calculate multipath
[MP1, CMC1] = mpath(C1C, L1C, L1, L2W, L2);
[MP2a, CMC2a] = mpath(C2L, L2L, L2, L1C, L1);
[MP2b, CMC2b] = mpath(C2W, L2W, L2, L1C, L1);
[MP5, CMC5] = mpath(C5Q, L5Q, L5, L1C, L1);

% Plot
figure;
subplot(2,1,1)
plot(time,CMC1,'b')
grid on
hold on
plot(time,CMC2a,'r')
plot(time,CMC2b,'g')
plot(time,CMC5,'k')
xlabel('Hours into August 19, 2026');
ylabel('CMC (m)')
legend('C1C','C2L','C2W','C5Q')
title('Code Minus Carrier')

subplot(2,1,2)
plot(time,S1C,'.b')
grid on
hold on
plot(time,S2L,'.r')
plot(time,S2W,'.g')
plot(time,S5Q,'.k')
xlabel('Hours into August 19, 2026');
ylabel('SNR (dB-Hz)')
legend('S1C','S2L','S2W','S5Q')
title('Signal-to-Noise Ratio')

figure;
plot(time,MP1,'b')
hold on
grid on
plot(time,MP2a,'r')
plot(time,MP2b,'g')
plot(time,MP5,'k')
xlabel('Hours into August 19, 2026');
ylabel('MP (m)')
legend('C1C','C2L','C2W','C5Q')
title('Multipath')