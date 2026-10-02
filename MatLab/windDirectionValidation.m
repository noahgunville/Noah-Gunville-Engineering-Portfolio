% wind direction convention check
%% A. Unit test of the sounding conversion (meteorological "from" -> "toward")
% Same formulas as your sounding code: motionDir = mod(dirFrom+180,360)
convE = @(spd, dirFrom) spd .* sind(mod(dirFrom + 180, 360));
convN = @(spd, dirFrom) spd .* cosd(mod(dirFrom + 180, 360));
status = ["FAIL", "PASS"];

%        speed  dirFrom  expectE  expectN
cases = [ 10     270      10       0;      % from west  -> blows east
    10      90     -10       0;      % from east  -> blows west
    10     180       0      10;      % from south -> blows north
    10       0       0     -10;      % from north -> blows south
    10     225    7.071   7.071];    % from SW    -> blows NE

for k = 1:size(cases,1)
    e = convE(cases(k,1), cases(k,2));
    n = convN(cases(k,1), cases(k,2));
    ok = abs(e - cases(k,3)) < 0.01 && abs(n - cases(k,4)) < 0.01;
    fprintf('%s  wind from %3d deg -> E=%6.2f, N=%6.2f\n', ...
        status(double(ok)+1), cases(k,2), e, n);
end

%% B. HWM component order and sign (climatology sanity check)
% HWM returns [meridional(north), zonal(east)], already "toward" convention.
% At ~41N in December the stratospheric polar-night jet is strongly EASTWARD,
% and meridional winds are comparatively weak. A swapped column order or a
% flipped sign breaks this immediately.
band = atmoAlt >= 20000 & atmoAlt <= 60000;
meanE = mean(windEast(band));
meanN = mean(windNorth(band));
fprintf('\nHWM 20-60 km mean: East = %.1f m/s, North = %.1f m/s\n', meanE, meanN);
ok = meanE > 0 && abs(meanE) > abs(meanN);
fprintf('%s  HWM zonal wind is eastward and dominant in winter stratosphere\n', ...
    status(double(ok)+1));

%% C. Compare HWM against your sounding (the strongest convention test)
% HWM is climatological, so expect similar structure, not identical values.
% A sign or column error shows up as NEGATIVE correlation.
data = readtable("RenoData1208Extracted.csv");   % use your Black Rock sounding if you have it
altS = data.geopotentialHeightM;
spd  = data.windSpeedm_s;
dirS = data.windDirectionDegree;
eS = convE(spd, dirS);
nS = convN(spd, dirS);

valid = isfinite(altS) & isfinite(eS) & isfinite(nS) & altS >= atmoAlt(end) & altS <= 30000;
[altSv, ord] = sort(altS(valid));  eSv = eS(valid); eSv = eSv(ord);
nSv = nS(valid); nSv = nSv(ord);

hwmE = interp1(flipud(atmoAlt), flipud(windEast),  altSv);
hwmN = interp1(flipud(atmoAlt), flipud(windNorth), altSv);

cE = corrcoef(eSv, hwmE); cN = corrcoef(nSv, hwmN);
fprintf('\nSounding vs HWM (up to 30 km): corr East = %.2f, corr North = %.2f\n', cE(1,2), cN(1,2));
fprintf('Mean difference (HWM - sounding): East = %.1f m/s, North = %.1f m/s\n', ...
    mean(hwmE - eSv), mean(hwmN - nSv));
fprintf('%s  East components positively correlated\n', status(double(cE(1,2) > 0)+1));

figure; tiledlayout(1,2);
nexttile; plot(eSv, altSv/1000, hwmE, altSv/1000); grid on
xlabel('East wind (m/s)'); ylabel('Altitude (km)'); legend('Sounding','HWM14');
nexttile; plot(nSv, altSv/1000, hwmN, altSv/1000); grid on
xlabel('North wind (m/s)'); legend('Sounding','HWM14');