function [landingX, landingY, driftDistance, results] = descentSimulation( ...
    mass_kg, ...
    apogee_m, ...
    drogueCd, ...
    drogueDiameter_m, ...
    mainCd, ...
    mainDiameter_m, ...
    mainDeploy_m, ...
    weather)

% define variables
g = 9.81;
x = 0;
y = 0;
totalTime = 0;
h = weather.altitude;

% Store results
results.altitude = [];
results.velocity = [];
results.time = [];
results.x = [];
results.y = [];




% Loop through every altitude interval
for i = 1:length(h)-1

    % Determine parachute
    if h(i) > mainDeploy_m
        Cd = drogueCd;
        diameter = drogueDiameter_m;
    else
        Cd = mainCd;
        diameter = mainDiameter_m;
    end

    % Parachute area
    A = pi*(diameter/2)^2;

    % Atmospheric density
    rho = weather.density(i);

    % Terminal velocity
    Vt = sqrt((2*mass_kg*g)/(rho*Cd*A));

    % Altitude change
    dh = h(i) - h(i+1);

    % Time spent in this altitude interval
    dt = dh/Vt;

    % Wind
    windEast = weather.windEast(i);
    windNorth = weather.windNorth(i);

    % Horizontal movement
    dx = windEast*dt;
    dy = windNorth*dt;

    % Update position
    x = x + dx;
    y = y + dy;

    % Update time
    totalTime = totalTime + dt;

    % Store results
    results.altitude(end+1) = h(i);
    results.velocity(end+1) = Vt;
    results.time(end+1) = totalTime;
    results.x(end+1) = x;
    results.y(end+1) = y;

end

% Final landing location
landingX = x
landingY = y

driftDistance = sqrt(landingX^2 + landingY^2)

end
    