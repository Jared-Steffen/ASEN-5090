clc; clear; close all

% Jared is working on Windows and Brady is working on Mac.
Windowsflag = ispc;
Macflag = ismac;

if Windowsflag
    addpath("Functions\");
    addpath("Data\");
elseif Macflag
    addpath("Functions/");
    addpath("Data/");
end

set_constants2025;

% Authors: Brady Sivey and Jared Steffen

%% Flags
NISTflag = true; %Switch to false if you would like to use NYA2

if NISTflag
    name = 'NIST';
    obsFile = 'NIST00USA_R_20262310000_01D_30S_MO.rnx';
    receiverECEF = [-1288398.567; -4721696.932; 4078625.350];
    zd = 2;        % zenith delay
else
    name = 'NYA2';
    obsFile = 'NYA200NOR_R_20262310000_01D_30S_MO.rnx';
    receiverECEF = [1202379.3100; 252474.6543; 6237786.5417];
    zd = 2.4;      % zenith delay
end

%% prep
referenceDate = datetime(2026,8,19);
eph = read_clean_GPSbroadcast('brdc2310.26n', true);

obs = rinexread(obsFile);
gps = obs.GPS;

obsTime = gps.Properties.RowTimes;
epochs = unique(obsTime);

time = hours(epochs - referenceDate);
N = numel(epochs);

receiverLLA = ecef2lla(receiverECEF'); %converting to lat-lon
C = ECEF2ENU(receiverLLA(1), receiverLLA(2)); %ECEF to ENU rotation matrix

%preallocating
prefit = cell(N,1);
postfit = cell(N,1);
elevation = cell(N,1);
residualTime = cell(N,1);
ENU = nan(N,3);
clockBias = nan(N,1);
DOP = nan(N,4);
numUsed = zeros(N,1);
discarded = zeros(N,1); 

for k = 1:N

    epochObs = gps(obsTime == epochs(k), :);
    numObserved = height(epochObs);

    PRN = epochObs.SatelliteID;
    n = numel(PRN);

    %ionospheric delay
    [PRIF, ~] = ionocorr(epochObs.C1C, L1, epochObs.C2W, L2);
    PRIF = PRIF(:);

    %August 19, 2026 begins at GPS week 2432, TOW 259200 s.
    gpsTime = [2432, 259200 + seconds(epochs(k)-referenceDate)];

    %Allocation for each epoch
    R = nan(n,1);
    BSV = nan(n,1);
    RELAT = nan(n,1);
    EL = nan(n,1);
    satECEF = nan(3,n);

    for j = 1:n

        % Initial satellite position evaluated at reception time.
        [health, satReceive] = eph2pvt2025(eph, gpsTime, PRN(j));

        %health check
        if health ~= 0 || any(~isfinite(satReceive(:)))
            continue
        end

        % Account for satellite motion during signal travel and Earth's
        [R(j), satECEF(:,j)] = correct_range(gpsTime, receiverECEF, satReceive(:), eph, PRN(j), true);

        % Transmission occurs R/c seconds before reception.
        txTime = gpsTime;
        txTime(2) = txTime(2) - R(j)/c;

        [healthTx, ~, ~, BSV(j), RELAT(j)] = eph2pvt2025(eph, txTime, PRN(j));

        %health check
        if healthTx ~= 0
            continue
        end

        % Elevation measured above the receiver's local horizon.
        [~, EL(j)] = ecef2az_el_range(satECEF(:,j), receiverECEF);
    end

    %retain if elevation >= 10 deg
    use = EL >= 10 & isfinite(R) & isfinite(BSV) & isfinite(RELAT) & all(isfinite(satECEF),1)';

    %tropo delay
    tropo = nan(n,1);
    tropo(use) = tropocalc(zd, EL(use));

    %corrections
    dpr = PRIF - (R - BSV - RELAT + tropo);

    %getting rid of outliers
    if any(use)
        use = use & abs(dpr-median(dpr(use))) <= 100;
    end

    %tracking discards
    discarded(k) = numObserved - sum(use);

    %geometry matrix using the retained satellites.
    G = measSensMat(satECEF(:,use), receiverECEF, R(use), true);

    %correction calc
    correction = G\dpr(use);
    
    %rotate
    ENU(k,:) = (C*correction(1:3))';

    clockBias(k) = correction(4);
    numUsed(k) = sum(use);

    %residuals/pre and post fits
    prefit{k} = dpr(use);
    postfit{k} = prefit{k} - G*correction;
    elevation{k} = EL(use);
    residualTime{k} = repmat(time(k), sum(use), 1);

    %Geometry covariance
    Q = (G'*G)\eye(4);
    Qenu = C*Q(1:3,1:3)*C';

    %4 DOP values
    DOP(k,:) = sqrt([diag(Qenu); Q(4,4)])';
end

%stacking
pre = vertcat(prefit{:});
post = vertcat(postfit{:});
el = vertcat(elevation{:});
tResidual = vertcat(residualTime{:});

%solved epochs
solved = numUsed > 0;

%in case none are solved
assert(any(solved), 'No valid solutions for %s.', name);

fprintf('\n%s \n', name);
fprintf('Solved epochs: %d of %d\n', sum(solved), N);
fprintf('Measurements discarded: %d\n', sum(discarded));
fprintf('Zenith tropospheric delay: %.2f m\n', zd);

%% problem 1

%01:00 GPS time
idx = find(epochs == referenceDate + hours(1), 1);

if ~isempty(idx) && solved(idx)

    fprintf('\nProblem 1: %s at 01:00 GPS time\n', name);

    %Print position corrections
    fprintf('East = %.2f m\n', ENU(idx,1));
    fprintf('North = %.2f m\n', ENU(idx,2));
    fprintf('Up = %.2f m\n', ENU(idx,3));

    %Clock bias
    fprintf('Clock bias = %.2f m = %.9e s\n', clockBias(idx), clockBias(idx)/c);

    fprintf(['EDOP = %.4f, NDOP = %.4f, ' 'VDOP = %.4f, TDOP = %.4f\n'], DOP(idx,:));

else

    warning('%s has no valid solution at 01:00.', name);

end

%% problem 2

%stats
values = {ENU(solved,1), ENU(solved,2), ENU(solved,3), pre, post};

stdValues = zeros(5,1);
rmsValues = zeros(5,1);

for j = 1:5

    %Standard deviation
    stdValues(j) = std(values{j});

    %RMS
    rmsValues(j) = sqrt(mean(values{j}.^2));

end

component = {'East'; 'North'; 'Up'; 'Pre-fit'; 'Post-fit'};

statistics = table(component, stdValues, rmsValues, 'VariableNames', {'Component','Std_m','RMS_m'});

fprintf('\nProblem 2: %s full-day statistics\n', name);
disp(statistics)

%% problem 3
figure;

subplot(2,1,1)
plot(tResidual, pre, 'b.');
grid on;
xlabel('Hours into August 19, 2026 (GPS time)');
ylabel('Pre-fit residual (m)');
title([name, ' Pre-fit Residuals']);

subplot(2,1,2)
plot(el, pre, 'b.');
grid on;
xlabel('Elevation (deg)');
ylabel('Pre-fit residual (m)');

%% problem 4
figure;

subplot(2,1,1)
plot(tResidual, post, 'b.');
grid on;
xlabel('Hours into August 19, 2026 (GPS time)');
ylabel('Post-fit residual (m)');
title([name, ' Post-fit Residuals']);

subplot(2,1,2)
plot(el, post, 'b.');
grid on;
xlabel('Elevation (deg)');
ylabel('Post-fit residual (m)');
