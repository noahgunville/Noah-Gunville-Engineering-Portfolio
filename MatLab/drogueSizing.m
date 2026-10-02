
%%parameters for descent sim
dayLabel = '2026-12-05';        
drogueFt = 4:0.5:16;  % drogue diameters to test (ft) (4-16 ft at 0.5 ft increments)

vehicle.mass      = 138; % kg
vehicle.apogee    = 100000; % m
vehicle.drogueCd  = 0.97;
vehicle.mainCd    = 2.2;
vehicle.mainD     = 6.1;  % m (20 ft)
vehicle.deployAlt = 2105; % m MSL 

maxDeploySpeed = NaN;             % m/s: largest drogue speed allowed at main deploy
                                  % (opening-shock limit). NaN = no limit.

%% ---------------- Checks ----------------
if ~exist('weather_HWM', 'var')
    error('weather_HWM is not in the workspace. Run your weather script first.');
end
wx = weather_HWM;
assert(all(isfinite([wx.altitude(:); wx.density(:); wx.windEast(:); wx.windNorth(:)])), ...
    'weather_HWM contains NaN/Inf values');
assert(all(diff(wx.altitude(:)) < 0), 'weather_HWM.altitude must be strictly descending');
assert(vehicle.deployAlt < wx.altitude(1) && vehicle.deployAlt > wx.altitude(end), ...
    'Main deploy altitude is outside the weather profile (is it MSL?)');

%% -Run sim on each drogue size
nD = numel(drogueFt);
landE = nan(1, nD);  landN = nan(1, nD);  drift = nan(1, nD);
tTotal = nan(1, nD); vDeploy = nan(1, nD); vLand = nan(1, nD);

for j = 1:nD
    Dm = drogueFt(j) * 0.3048;                                        % ft -> m
    [lx, ly, d, r] = descentSimulation(vehicle.mass, vehicle.apogee, ...
        vehicle.drogueCd, Dm, vehicle.mainCd, vehicle.mainD, vehicle.deployAlt, wx);
    landE(j) = lx;  landN(j) = ly;  drift(j) = d;
    tTotal(j) = r.time(end);
    vLand(j)  = r.velocity(end);
    iD = find(r.altitude > vehicle.deployAlt, 1, 'last');             % last drogue interval
    vDeploy(j) = r.velocity(iD);
end
bearing = mod(atan2d(landE, landN), 360);                             % deg from north

%% ---------------- Best size ----------------
feasible = true(1, nD);
if ~isnan(maxDeploySpeed), feasible = vDeploy <= maxDeploySpeed; end
driftF = drift;  driftF(~feasible) = Inf;
[bestDrift, ib] = min(driftF);
if isfinite(bestDrift), bestFt = drogueFt(ib); else, bestFt = NaN; end

below = wx.altitude < 35000;
peakWind = max(hypot(wx.windEast(below), wx.windNorth(below)));       % m/s below 35 km

% Linear fit check: in this model each drift component is linear in diameter
pE = polyfit(drogueFt, landE, 1);
pN = polyfit(drogueFt, landN, 1);
linErr = max([max(abs(polyval(pE, drogueFt) - landE)), max(abs(polyval(pN, drogueFt) - landN))]);
Dstar = -(pE(1)*pE(2) + pN(1)*pN(2)) / (pE(1)^2 + pN(1)^2);           % ft, unconstrained minimum
if Dstar < drogueFt(1) || Dstar > drogueFt(end), Dstar = NaN; end
slopeKmPerFt = (drift(end) - drift(1)) / (drogueFt(end) - drogueFt(1)) / 1000;


%% ---------------- Max drogue size under the drift limit ----------------
limitKm    = 16;     % drift limit (km)
marginFrac = 0;      % safety margin on predicted drift (0.30 = 30%); 0 = none

Rm = limitKm*1000 / (1 + marginFrac);                      % allowed predicted drift (m)
a = pE(1)^2 + pN(1)^2;
b = 2*(pE(1)*pE(2) + pN(1)*pN(2));
c = pE(2)^2 + pN(2)^2 - Rm^2;
if c > 0
    DmaxFt = NaN;                                          % limit exceeded even at zero diameter
elseif a == 0
    DmaxFt = Inf;                                          % no drogue-phase drift
else
    DmaxFt = (-b + sqrt(b^2 - 4*a*c)) / (2*a);             % largest D where drift = limit
end
DmaxRound = floor(DmaxFt/0.5) * 0.5;                       % round down to nearest 0.5 ft

fprintf('Max drogue size for drift < %.1f km (margin %.0f%%): %.1f ft (%.1f ft rounded down, %.2f m)\n', ...
    limitKm, 100*marginFrac, DmaxFt, DmaxRound, DmaxRound*0.3048);
if DmaxFt > drogueFt(end)
    fprintf('  Note: this is above your largest tested size (%.0f ft), so it relies on the straight-line trend.\n', drogueFt(end));
elseif isnan(DmaxFt) || DmaxFt < drogueFt(1)
    fprintf('  Note: no tested size meets the limit on this day.\n');
end

%% ---------------- Table ----------------
fprintf('\n=== Drogue sweep: %s ===\n', dayLabel);
fprintf('%8s %8s %11s %10s %10s %9s %10s %12s %9s\n', 'Dia (ft)', 'Dia (m)', 'Drift (m)', ...
    'East (m)', 'North (m)', 'Brg (deg)', 'Time (min)', 'V@main (m/s)', 'OK?');
for j = 1:nD
    if feasible(j), ok = 'yes'; else, ok = 'NO'; end
    fprintf('%8.1f %8.2f %11.0f %10.0f %10.0f %9.0f %10.2f %12.1f %9s\n', drogueFt(j), ...
        drogueFt(j)*0.3048, drift(j), landE(j), landN(j), bearing(j), tTotal(j)/60, vDeploy(j), ok);
end
fprintf('\nPeak wind below 35 km: %.1f m/s\n', peakWind);
fprintf('Best size (min drift%s): %.1f ft -> %.2f km\n', ...
    ternary(isnan(maxDeploySpeed), '', sprintf(', V@main <= %.0f m/s', maxDeploySpeed)), bestFt, bestDrift/1000);
fprintf('Drift sensitivity: %.3f km per ft of drogue diameter (%.2f km at %.0f ft -> %.2f km at %.0f ft)\n', ...
    slopeKmPerFt, drift(1)/1000, drogueFt(1), drift(end)/1000, drogueFt(end));
fprintf('Unconstrained minimum inside range: %s\n', ...
    ternary(isnan(Dstar), 'none (minimum is at the smallest or largest size)', sprintf('%.1f ft', Dstar)));
fprintf('Linearity check (max deviation from a straight line): %.1e m\n', linErr);
fprintf('Landing speed under main: %.2f m/s (set by the main, not the drogue)\n\n', vLand(1));

%% ---------------- Plots ----------------
c1 = [0.00 0.45 0.74];  c2 = [0.85 0.33 0.10];
f = figure('Name', ['Drogue sweep ' dayLabel], 'Position', [60 60 1500 800]);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Drogue diameter sweep: %s  (peak wind %.0f m/s)', dayLabel, peakWind));

% 1: drift vs diameter
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, drogueFt, drift/1000, '-o', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', 'Drift');
%if isfinite(bestFt)
   % plot(ax, bestFt, bestDrift/1000, 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'y', ...
     %   'MarkerEdgeColor', 'k', 'DisplayName', 'Best');
%end
grid(ax, 'on'); xlabel(ax, 'Drogue diameter (ft)'); ylabel(ax, 'Drift distance (km)');
title(ax, 'Drift distance'); legend(ax, 'Location', 'best');

% 2: east / north components
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, drogueFt, landE/1000, '-o', 'Color', c1, 'LineWidth', 1.4, 'DisplayName', 'East');
plot(ax, drogueFt, landN/1000, '-s', 'Color', c2, 'LineWidth', 1.4, 'DisplayName', 'North');
yline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'Drogue diameter (ft)'); ylabel(ax, 'Landing offset (km)');
title(ax, 'Landing components'); legend(ax, 'Location', 'best');

% 3: drogue speed at main deploy
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, drogueFt, vDeploy, '-o', 'Color', c1, 'LineWidth', 1.5, 'DisplayName', 'Drogue speed');
if ~isnan(maxDeploySpeed)
    yline(ax, maxDeploySpeed, '--r', 'limit', 'HandleVisibility', 'off');
end
grid(ax, 'on'); xlabel(ax, 'Drogue diameter (ft)'); ylabel(ax, 'Speed at main deploy (m/s)');
title(ax, 'Drogue descent speed at main deploy');

% 4: total descent time
ax = nexttile(tl);
plot(ax, drogueFt, tTotal/60, '-o', 'Color', c1, 'LineWidth', 1.5);
grid(ax, 'on'); xlabel(ax, 'Drogue diameter (ft)'); ylabel(ax, 'Descent time (min)');
title(ax, 'Total descent time');

% 5: landing points colored by diameter
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, landE/1000, landN/1000, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
scatter(ax, landE/1000, landN/1000, 55, drogueFt, 'filled');
plot(ax, 0, 0, 'o', 'Color', [0.1 0.6 0.1], 'MarkerSize', 9, 'LineWidth', 2);
axis(ax, 'equal'); grid(ax, 'on');
cb = colorbar(ax); cb.Label.String = 'Drogue diameter (ft)';
xlabel(ax, 'East (km)'); ylabel(ax, 'North (km)'); title(ax, 'Landing points (o = apogee point)');

% 6: the day's wind below 35 km, for context
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, wx.windEast(below),  wx.altitude(below)/1000, 'Color', c1, 'LineWidth', 1.3, 'DisplayName', 'East');
plot(ax, wx.windNorth(below), wx.altitude(below)/1000, 'Color', c2, 'LineWidth', 1.3, 'DisplayName', 'North');
xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
yline(ax, vehicle.deployAlt/1000, '--k', 'main deploy', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'Wind (m/s)'); ylabel(ax, 'Altitude (km MSL)');
title(ax, 'Wind below 35 km'); legend(ax, 'Location', 'best');


%% ---------------- Local helper ----------------
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end