% LLA2ECEF - convert latitude, longitude, and altitude to
%            earth-centered, earth-fixed (ECEF) cartesian
% 
% USAGE:
% [x,y,z] = lla2ecef(lat,lon,alt)
% 
% x = ECEF X-coordinate (m)
% y = ECEF Y-coordinate (m)
% z = ECEF Z-coordinate (m)
% lat = geodetic latitude (radians)
% lon = longitude (radians)
% alt = height above WGS84 ellipsoid (m)
% 
% Notes: This function assumes the WGS84 model.
%        Latitude is customary geodetic (not geocentric).
% 
% Source: "Department of Defense World Geodetic System 1984"
%         Page 4-4
%         National Imagery and Mapping Agency
%         Last updated June, 2004
%         NIMA TR8350.2
% 
% Michael Kleder, July 2005

function ecef =lla2ecef(lla)

% WGS84 ellipsoid constants:
a = 6378137;
e = 8.1819190842622e-2;

% unpack input
if size(lla,2) == 3
    lat = lla(:,1) * pi / 180; % geodetic latitude (radians)
    lon = lla(:,2) * pi / 180; % longitude (radians)
    alt = lla(:,3); % height above WGS84 ellipsoid (m)
elseif size(lla,1) == 3
    lat = lla(1)* pi / 180; % geodetic latitude (radians)
    lon = lla(2)* pi / 180; % longitude (radians)
    alt = lla(3); % height above WGS84 ellipsoid (m)
else
    error('Input must be a 3-element vector or Nx3 matrix');
end
% intermediate calculation
% (prime vertical radius of curvature)
N = a ./ sqrt(1 - e^2 .* sin(lat).^2);

% results:
x = (N+alt) .* cos(lat) .* cos(lon);
y = (N+alt) .* cos(lat) .* sin(lon);
z = ((1-e^2) .* N + alt) .* sin(lat);
ecef = [x(:), y(:), z(:)];
return