data = readtable("RenoData12_06_00.csv");

%data from CSV
altitude = data.geopotentialHeightM;
pressure = data.pressurehPa;
temperature = data.temperatureC;
windDir = data.windDirectionDegree;
windSpeed = data.windSpeedm_s;

%calculate density
pressurePa = pressure .* 100;
temperatureK = temperature + 273.15;
R = 287.05;
density = pressurePa ./ (R .* temperatureK);

%convert wind direction to E/W
motionDir = mod(windDir + 180,360);
theta = deg2rad(motionDir);
windEast = windSpeed .* sind(motionDir);
windNorth = windSpeed .* cosd(motionDir);

%clean sounding: drop bad rows, drop the part overlapping the upper layer,
%sort descending, and remove duplicate altitudes
valid = isfinite(altitude) & isfinite(density) & isfinite(windEast) & ...
        isfinite(windNorth) & altitude < 35000;
altitude  = altitude(valid);
density   = density(valid);
windEast  = windEast(valid);
windNorth = windNorth(valid);

[altitude, order] = sort(altitude, 'descend');
density   = density(order);
windEast  = windEast(order);
windNorth = windNorth(order);

[altitude, iu] = unique(altitude, 'stable');   % keeps the first of each duplicate
density   = density(iu);
windEast  = windEast(iu);
windNorth = windNorth(iu);

%vars for descentSimulation
weather_data.altitude = altitude;
weather_data.density = density;
weather_data.windEast = windEast;
weather_data.windNorth = windNorth;

%create upper atmosphere layers
upperAlt = (100000:-100:35000)';



%HWM parameters
lat = 40.91;
lon = -119.056;
latVec = lat*ones(size(upperAlt));
lonVec = lon*ones(size(upperAlt));


day = 340; %dec 8: 342 (predicted launch day). used other days for variability checks
seconds = 0; % 8AM UTC in black rock: 57600. 12UTC : 43200
dayVec     = day   * ones(size(upperAlt));
secondsVec = seconds * ones(size(upperAlt));

%run HWM
windHWM = atmoshwm(latVec, lonVec, upperAlt, 'day',dayVec,'seconds',secondsVec, 'version','14');
upperWindEast_HWM = windHWM(:,2);
upperWindNorth_HWM = windHWM(:,1);

%density model parameters
lstVec = mod(secondsVec/3600 + lonVec/15, 24);   % local solar time in hours
f107A = 130 * ones(size(upperAlt));      % 81-day average F10.7
f107  = 130 * ones(size(upperAlt));      % daily F10.7
aph   = 4   * ones(numel(upperAlt), 7);  % M-by-7 Ap matrix

%density model
year = 2026;                                  
yearVec = year * ones(size(upperAlt));
[~, rhoNRL] = atmosnrlmsise00(upperAlt, latVec, lonVec, yearVec, dayVec, secondsVec);
upperDensity = rhoNRL(:,6);    % total mass density, kg/m^3



%combinded data for sim
weather_HWM.altitude = [upperAlt; weather_data.altitude];
weather_HWM.windEast = [upperWindEast_HWM; weather_data.windEast];
weather_HWM.windNorth = [upperWindNorth_HWM; weather_data.windNorth];
weather_HWM.density = [upperDensity; weather_data.density];

%sim test
[landingX, landingY, driftDistance, results] = descentSimulation( ...
    138, ... % mass kg
    100000, ... % apogee m
    0.97, ... % drogue Cd (example)
    3.05, ...% drogue diameter m (example)
    2.2, ... % main Cd
    7.3, ... % main diameter m (24 ft)
    2105, ... % deploy altitude m
    weather_HWM);

resComb = results;
wxComb  = weather_HWM;