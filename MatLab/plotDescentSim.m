function plotDescentSim(resA, wxA, resB, wxB, labelA, labelB, outPrefix)
%PLOTDESCENTCOMPARISON  Plot and tabulate one or two descentSimulation runs.
%
%   plotDescentSim(resA, wxA)                       one run
%   plotDescentSim(resA, wxA, resB, wxB)            two runs (A solid, B dashed)
%   plotDescentSim(resA, wxA, resB, wxB, labelA, labelB)
%   plotDescentSim(..., outPrefix)                  also saves PNGs
%
%   resA/resB : 4th output of descentSimulation (the "results" struct)
%   wxA/wxB   : the weather struct that was passed to that run
%   labelA/B  : legend names, e.g. 'Sounding + HWM', 'HWM only'
%   outPrefix : optional file prefix, e.g. 'dec08' -> dec08_atmosphere.png, dec08_drift.png
%
%   Figure 1: east/north wind, wind speed, wind difference, density,
%             descent speed, elapsed time (all vs altitude)
%   Figure 2: east/north drift, ground track, distance from launch,
%             drift difference, drift per altitude band
%   Also prints a summary table and a drift-by-band table to the command window.
%
%   Positions are cumulative from the start of the descent, so each value is
%   attributed to the altitude at the BOTTOM of the interval it was computed for.

if nargin < 3, resB = []; wxB = []; end
if nargin < 5 || isempty(labelA), labelA = 'Model A'; end
if nargin < 6 || isempty(labelB), labelB = 'Model B'; end
if nargin < 7, outPrefix = ''; end

runs = {prepRun(resA, wxA, labelA)};
if ~isempty(resB)
    runs{2} = prepRun(resB, wxB, labelB);
end
n = numel(runs);
two = (n == 2);

colors = [0.00 0.45 0.74; 0.85 0.33 0.10];
styles = {'-', '--'};
eColor = [0.20 0.60 0.20];     % east component in difference plots
nColor = [0.50 0.20 0.60];     % north component in difference plots

if two
    hdr = sprintf('A (solid): %s   |   B (dashed): %s', runs{1}.label, runs{2}.label);
else
    hdr = runs{1}.label;
end

% ---- common altitude grid (only needed to difference two runs) ----
if two
    gTop = min(runs{1}.h(1),   runs{2}.h(1));
    gBot = max(runs{1}.h(end), runs{2}.h(end));
    g = (gTop:-100:gBot)';
    dWindE  = onGrid(runs{1}.h, runs{1}.wE, g) - onGrid(runs{2}.h, runs{2}.wE, g);
    dWindN  = onGrid(runs{1}.h, runs{1}.wN, g) - onGrid(runs{2}.h, runs{2}.wN, g);
    dDriftE = onGrid(runs{1}.h, runs{1}.x,  g) - onGrid(runs{2}.h, runs{2}.x,  g);
    dDriftN = onGrid(runs{1}.h, runs{1}.y,  g) - onGrid(runs{2}.h, runs{2}.y,  g);
end

% ---- altitude bands (MSL, m) for drift accumulation ----
lo = max(cellfun(@(r) r.h(end), runs));
hi = min(cellfun(@(r) r.h(1),   runs));
cand = [3000 6000 9000 12000 20000 35000 60000];
edges = [lo, cand(cand > lo & cand < hi), hi];           % ascending
nb = numel(edges) - 1;
bandLabels = cell(nb, 1);
for i = 1:nb
    bandLabels{i} = sprintf('%.1f-%.1f', edges(i)/1000, edges(i+1)/1000);
end
band = struct('dE', cell(1,n), 'dN', cell(1,n), 't', cell(1,n));
for k = 1:n
    r  = runs{k};
    xe = onGrid(r.h, r.x, edges');
    ye = onGrid(r.h, r.y, edges');
    te = onGrid(r.h, r.t, edges');
    % cumulative values grow toward lower altitude, so band = lower edge - upper edge
    band(k).dE = xe(1:end-1) - xe(2:end);
    band(k).dN = ye(1:end-1) - ye(2:end);
    band(k).t  = te(1:end-1) - te(2:end);
end

%% ================= Figure 1: atmosphere and descent speed =================
f1 = figure('Name', 'Atmosphere and descent', 'Position', [40 60 1500 700]);
tl = tiledlayout(f1, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, hdr, 'Interpreter', 'none');

ax = nexttile(tl);
addLines(ax, runs, 'wE', 'h', 1, 1000, colors, styles);
xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
styleAxes(ax, 'East wind (m/s)', 'Altitude (km MSL)', 'East wind', two);

ax = nexttile(tl);
addLines(ax, runs, 'wN', 'h', 1, 1000, colors, styles);
xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
styleAxes(ax, 'North wind (m/s)', 'Altitude (km MSL)', 'North wind', two);

ax = nexttile(tl);
addLines(ax, runs, 'wS', 'h', 1, 1000, colors, styles);
styleAxes(ax, 'Wind speed (m/s)', 'Altitude (km MSL)', 'Wind speed', two);

if two
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, dWindE, g/1000, '-',  'Color', eColor, 'LineWidth', 1.4, 'DisplayName', 'East');
    plot(ax, dWindN, g/1000, '--', 'Color', nColor, 'LineWidth', 1.4, 'DisplayName', 'North');
    xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    styleAxes(ax, 'Wind difference, A - B (m/s)', 'Altitude (km MSL)', 'Wind difference', true);
end

ax = nexttile(tl);
addLines(ax, runs, 'rho', 'h', 1, 1000, colors, styles);
set(ax, 'XScale', 'log');
styleAxes(ax, 'Density (kg/m^3)', 'Altitude (km MSL)', 'Air density', two);

ax = nexttile(tl);
addLines(ax, runs, 'v', 'hv', 1, 1000, colors, styles);
set(ax, 'XScale', 'log');
styleAxes(ax, 'Descent speed (m/s)', 'Altitude (km MSL)', 'Terminal velocity', two);

ax = nexttile(tl);
addLines(ax, runs, 't', 'h', 60, 1000, colors, styles);
styleAxes(ax, 'Elapsed time (min)', 'Altitude (km MSL)', 'Descent time', two);

%% ================= Figure 2: drift =================
f2 = figure('Name', 'Drift', 'Position', [60 40 1500 800]);
tl2 = tiledlayout(f2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, hdr, 'Interpreter', 'none');

ax = nexttile(tl2);
addLines(ax, runs, 'x', 'h', 1000, 1000, colors, styles);
styleAxes(ax, 'East drift (km)', 'Altitude (km MSL)', 'Cumulative east drift', two);

ax = nexttile(tl2);
addLines(ax, runs, 'y', 'h', 1000, 1000, colors, styles);
xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
styleAxes(ax, 'North drift (km)', 'Altitude (km MSL)', 'Cumulative north drift', two);

ax = nexttile(tl2); hold(ax, 'on');
for k = 1:n
    r = runs{k};
    plot(ax, r.x/1000, r.y/1000, styles{k}, 'Color', colors(k,:), 'LineWidth', 1.6, ...
        'DisplayName', r.label);
    plot(ax, r.x(end)/1000, r.y(end)/1000, 'x', 'Color', colors(k,:), ...
        'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
end
plot(ax, 0, 0, 'o', 'Color', [0.1 0.6 0.1], 'MarkerSize', 8, 'LineWidth', 2, ...
    'DisplayName', 'Apogee point');
axis(ax, 'equal');
styleAxes(ax, 'East (km)', 'North (km)', 'Ground track (x = landing)', true);

ax = nexttile(tl2);
addLines(ax, runs, 'd', 'h', 1000, 1000, colors, styles);
styleAxes(ax, 'Distance from start (km)', 'Altitude (km MSL)', 'Cumulative horizontal drift', two);

if two
    ax = nexttile(tl2); hold(ax, 'on');
    plot(ax, dDriftE/1000, g/1000, '-',  'Color', eColor, 'LineWidth', 1.4, 'DisplayName', 'East');
    plot(ax, dDriftN/1000, g/1000, '--', 'Color', nColor, 'LineWidth', 1.4, 'DisplayName', 'North');
    xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    styleAxes(ax, 'Drift difference, A - B (km)', 'Altitude (km MSL)', ...
        'Drift difference (where the models diverge)', true);
end

ax = nexttile(tl2);
mag = zeros(nb, n);
for k = 1:n
    mag(:, k) = hypot(band(k).dE, band(k).dN) / 1000;
end
b = bar(ax, mag);
for k = 1:n
    b(k).FaceColor = colors(k,:);
    b(k).DisplayName = runs{k}.label;
end
set(ax, 'XTick', 1:nb, 'XTickLabel', bandLabels);
xtickangle(ax, 30);
styleAxes(ax, 'Altitude band (km MSL)', 'Horizontal drift in band (km)', ...
    'Drift accumulated per altitude band', two);

%% ================= Console summary =================
fprintf('\n=================== Descent summary ===================\n');
fprintf('%-34s', 'Metric');
for k = 1:n, fprintf('%-22s', runs{k}.label); end
fprintf('\n');

metrics = {'Landing east (m)', 'Landing north (m)', 'Drift distance (m)', ...
    'Bearing (deg from N)', 'Total descent time (min)', 'Landing speed (m/s)', ...
    'Max wind speed (m/s)', 'Altitude of max wind (m)', ...
    'Net drift below 20 km (m)'};
vals = nan(numel(metrics), n);
for k = 1:n
    r = runs{k};
    [wm, iw] = max(r.wS);
    vals(1,k) = r.x(end);
    vals(2,k) = r.y(end);
    vals(3,k) = hypot(r.x(end), r.y(end));
    vals(4,k) = mod(atan2d(r.x(end), r.y(end)), 360);
    vals(5,k) = r.t(end) / 60;
    vals(6,k) = r.v(end);
    vals(7,k) = wm;
    vals(8,k) = r.h(iw);
    if 20000 > r.h(end) && 20000 < r.h(1)
        x20 = onGrid(r.h, r.x, 20000);
        y20 = onGrid(r.h, r.y, 20000);
        vals(9,k) = hypot(r.x(end) - x20, r.y(end) - y20);
    end
end
for m = 1:numel(metrics)
    fprintf('%-34s', metrics{m});
    for k = 1:n, fprintf('%-22.1f', vals(m,k)); end
    fprintf('\n');
end
if two
    sep = hypot(runs{1}.x(end) - runs{2}.x(end), runs{1}.y(end) - runs{2}.y(end));
    fprintf('\nLanding separation between models: %.0f m (%.1f%% of A''s drift)\n', ...
        sep, 100 * sep / max(vals(3,1), eps));
end

fprintf('\n=================== Drift by altitude band ===================\n');
for k = 1:n
    fprintf('%s\n', runs{k}.label);
    fprintf('  %-12s %10s %11s %11s %11s\n', 'Band (km)', 'Time (s)', 'East (m)', 'North (m)', 'Horiz (m)');
    for i = 1:nb
        fprintf('  %-12s %10.1f %11.0f %11.0f %11.0f\n', bandLabels{i}, band(k).t(i), ...
            band(k).dE(i), band(k).dN(i), hypot(band(k).dE(i), band(k).dN(i)));
    end
end
fprintf('\n');

%% ================= Optional save =================
if ~isempty(outPrefix)
    try
        exportgraphics(f1, [outPrefix '_atmosphere.png'], 'Resolution', 200);
        exportgraphics(f2, [outPrefix '_drift.png'],      'Resolution', 200);
    catch err
        warning('Could not save figures: %s', err.message);
    end
end

end


%% ------------------------------------------------------------------------
function r = prepRun(res, wx, label)
    r.label = label;
    r.h   = wx.altitude(:);
    r.rho = wx.density(:);
    r.wE  = wx.windEast(:);
    r.wN  = wx.windNorth(:);
    r.wS  = hypot(r.wE, r.wN);
    % results hold one value per interval; prepend the start point so that each
    % cumulative value lines up with the altitude at the bottom of its interval
    r.x = [0; res.x(:)];
    r.y = [0; res.y(:)];
    r.t = [0; res.time(:)];
    r.d = hypot(r.x, r.y);
    r.hv = res.altitude(:);       % altitude at the top of each interval
    r.v  = res.velocity(:);
    if numel(r.x) ~= numel(r.h)
        error('plotDescentComparison:sizeMismatch', ...
            ['results and weather do not match for "%s" (%d intervals vs %d altitudes). ' ...
             'Pass the weather struct that produced these results.'], label, numel(res.x), numel(r.h));
    end
end

function v = onGrid(h, val, g)
    % interpolate val(h) onto g; h may be descending and may contain duplicates
    [hu, iu] = unique(h);
    v = interp1(hu, val(iu), g, 'linear');
end

function addLines(ax, runs, xf, yf, xScale, yScale, colors, styles)
    hold(ax, 'on');
    for k = 1:numel(runs)
        plot(ax, runs{k}.(xf) / xScale, runs{k}.(yf) / yScale, styles{k}, ...
            'Color', colors(k,:), 'LineWidth', 1.4, 'DisplayName', runs{k}.label);
    end
end

function styleAxes(ax, xl, yl, ttl, showLegend)
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl);
    if showLegend
        legend(ax, 'Location', 'best', 'Interpreter', 'none');
    end
end