function ok = validateSim(weather, params)
%VALIDATESIM  Validation checks for descentSimulation and its weather inputs.
%
%   ok = validateSim(weather)
%   ok = validateSim(weather, params)
%
%   weather : struct with column-vector fields (descending altitude, MSL)
%               .altitude   [m]
%               .density    [kg/m^3]
%               .windEast   [m/s]
%               .windNorth  [m/s]
%   params  : (optional) struct overriding vehicle defaults. Fields:
%               mass, drogueCd, drogueD, mainCd, mainD, deployAlt (MSL, m)
%
%   Returns ok = true if every check passed. Run it on any weather struct
%   (hybrid sounding+HWM, HWM-only, etc.) so it can never read stale
%   workspace variables.
%
%   Note: descentSimulation prints driftDistance if its last assignment has
%   no semicolon. Add one (driftDistance = sqrt(...);) to keep output clean.

if nargin < 2, params = struct(); end
p = setDefaults(params);

nFail = 0;   % shared with the nested check() helper

%% 0. Struct sanity
required = {'altitude','density','windEast','windNorth'};
for k = 1:numel(required)
    assert(isfield(weather, required{k}), 'weather is missing field "%s"', required{k});
end

atmoAlt   = weather.altitude(:);
density   = weather.density(:);
windEast  = weather.windEast(:);
windNorth = weather.windNorth(:);

fprintf('\n=== 1. Input data validity ===\n');
check(numel(density) == numel(atmoAlt) && numel(windEast) == numel(atmoAlt) && ...
      numel(windNorth) == numel(atmoAlt), 'all weather fields have the same length');
check(all(isfinite([atmoAlt; density; windEast; windNorth])), 'no NaN/Inf values');
check(all(density > 0), 'density strictly positive');
check(all(diff(atmoAlt) < 0), 'altitude strictly descending');

% Later checks are meaningless if the inputs are broken
if nFail > 0
    fprintf('\nInput data failed basic checks; stopping before the sim tests.\n');
    ok = false;
    return
end

fprintf('Altitude range: %.0f m (ground) to %.0f m (top), %d points\n', ...
    atmoAlt(end), atmoAlt(1), numel(atmoAlt));
fprintf('Ground density: %.3f kg/m^3 (expect ~1.0-1.15)\n', density(end));
fprintf('Top density:    %.2e kg/m^3 (about 5e-7 if the top is 100 km)\n', density(1));
[wmax, k] = max(hypot(windEast, windNorth));
fprintf('Max wind: %.1f m/s at %.0f m\n', wmax, atmoAlt(k));

% Density should fall smoothly: check for jumps between neighbors
jump = abs(diff(log(density)));
[jmax, jk] = max(jump);
fprintf('Largest log-density step: %.3f between %.0f and %.0f m\n', ...
    jmax, atmoAlt(jk), atmoAlt(jk+1));

%% 2. Density vs US Standard Atmosphere (valid to ~84 km)
fprintf('\n=== 2. Density vs COESA standard atmosphere ===\n');
m = atmoAlt <= 84000;
if any(m)
    [~,~,~,rhoStd] = atmoscoesa(atmoAlt(m));
    ratio = density(m) ./ rhoStd(:);
    fprintf('Density / COESA ratio: min %.2f, max %.2f\n', min(ratio), max(ratio));
    check(all(ratio > 0.5 & ratio < 2), ...
        'density within a factor of 2 of the standard atmosphere');
else
    disp('No points below 84 km; skipped.');
end

%% 3. Unit tests of the simulation with known answers (synthetic atmosphere)
fprintf('\n=== 3. Unit tests of descentSimulation ===\n');
top    = atmoAlt(1);
bottom = atmoAlt(end);
check(p.deployAlt < top && p.deployAlt > bottom, ...
    sprintf('deploy altitude %.0f m (MSL) is inside the profile range', p.deployAlt));

h = (top:-100:bottom)';
if h(end) ~= bottom, h(end+1) = bottom; end
rho0 = 1.0;  g = 9.81;
w.altitude = h;
w.density  = rho0 * ones(size(h));

sim = @(wx) descentSimulation(p.mass, top, p.drogueCd, p.drogueD, ...
    p.mainCd, p.mainD, p.deployAlt, wx);

% A: zero wind -> zero drift
w.windEast = zeros(size(h));  w.windNorth = zeros(size(h));
[~,~,d] = sim(w);
check(d < 1e-9, 'A: zero wind gives zero drift');

% B: uniform 10 m/s east wind -> x = 10 * total descent time
w.windEast = 10 * ones(size(h));
[x,y,~,r] = sim(w);
Vd = sqrt(2*p.mass*g / (rho0*p.drogueCd*pi*(p.drogueD/2)^2));
Vm = sqrt(2*p.mass*g / (rho0*p.mainCd*pi*(p.mainD/2)^2));
T  = (top - p.deployAlt)/Vd + (p.deployAlt - bottom)/Vm;
check(abs(r.time(end) - T)/T < 0.01, ...
    sprintf('B: descent time %.1f s matches analytic %.1f s (within 1%%)', r.time(end), T));
check(abs(x - 10*T)/(10*T) < 0.01 && abs(y) < 1e-9, ...
    sprintf('B: east drift %.1f m matches 10*T = %.1f m, north drift zero', x, 10*T));

% C: north wind only -> drift in y only
w.windEast = zeros(size(h));  w.windNorth = 10 * ones(size(h));
[x,y] = sim(w);
check(abs(x) < 1e-9 && y > 0, 'C: north wind gives drift in +y only');

% D: terminal velocity scales with sqrt(mass)
[~,~,~,r1] = descentSimulation(p.mass,   top, p.drogueCd, p.drogueD, p.mainCd, p.mainD, p.deployAlt, w);
[~,~,~,r2] = descentSimulation(4*p.mass, top, p.drogueCd, p.drogueD, p.mainCd, p.mainD, p.deployAlt, w);
check(abs(r2.velocity(1)/r1.velocity(1) - 2) < 1e-9, 'D: 4x mass doubles terminal velocity');

%% 4. Numerical convergence on the real weather data
fprintf('\n=== 4. Numerical convergence (step size study) ===\n');
steps  = [200 100 50 25];
drifts = zeros(size(steps));
for i = 1:numel(steps)
    hh = (top:-steps(i):bottom)';
    if hh(end) ~= bottom, hh(end+1) = bottom; end
    ww.altitude  = hh;
    ww.density   = exp(interp1(atmoAlt, log(density), hh));
    ww.windEast  = interp1(atmoAlt, windEast,  hh);
    ww.windNorth = interp1(atmoAlt, windNorth, hh);
    [~,~,drifts(i)] = descentSimulation(p.mass, top, p.drogueCd, p.drogueD, ...
        p.mainCd, p.mainD, p.deployAlt, ww);
    fprintf('step %3d m -> drift %.1f m\n', steps(i), drifts(i));
end
relChange = abs(drifts(end) - drifts(end-1)) / drifts(end);
check(relChange < 0.01, ...
    sprintf('drift changes %.2f%% between %d m and %d m steps (want < 1%%)', ...
    100*relChange, steps(end-1), steps(end)));

%% Summary
fprintf('\n');
if nFail == 0
    disp('ALL CHECKS PASSED');
else
    fprintf('%d CHECK(S) FAILED\n', nFail);
end
ok = (nFail == 0);

    function check(cond, msg)
        if cond
            fprintf('PASS  %s\n', msg);
        else
            fprintf('FAIL  %s\n', msg);
            nFail = nFail + 1;
        end
    end

end

function p = setDefaults(params)
    p = struct('mass', 63.5, 'drogueCd', 0.97, 'drogueD', 1.8, ...
               'mainCd', 2.2, 'mainD', 6.1, 'deployAlt', 2105);
    f = fieldnames(params);
    for k = 1:numel(f)
        p.(f{k}) = params.(f{k});
    end
end