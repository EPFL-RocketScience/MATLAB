function [T, a, p, rho] = correcStdAtmos(alt, hc, Tc, pc)
% stdAtmos
%
% INPUT:    - alt   : altitude above ground (hc) [m]
%           - hc    : altitude above see where correction values are
%           measured
%           - Tc    : local temperature
%           - pc    : local pressure
%
% OUPTUT:   - T     : local standard temperature [?K]
%           - a     : local speed of sound [m/s]
%           - p     : local standard pressure [Pa]
%           - rho   : local standard density [kg/m^3]
%
% ASSUMPTIONS:
% - hydrostatic approximation of atmosphere
% - linear temperature variation with altitude
% - homogenous composition
%
% LIMITATIONS:
% - troposphere: 10km
%
% AUTOR:    ERIC BRUNNER
% LAST UPDATE: 05/03/2017

% CHECK ALTITUDE RANGE
if alt > 1e4
    error('stdAtmos:outOfRange', 'The altitude is out of range: max 10km.')
end
    
% CONSTANTS
R = 287.04; %[M^2/?K/sec^2] real gas constant of air
gamma = 1.4; %[-] specific heat coefficient of air

% START CONDITIONS
p0 = pc; %[Pa]
rho0 = p0/R/T0; %[kg/m^3]
T0 = Tc; %[?K]
a0 = 340.294; %[m/sec]
g0 = 9.80665; %[m/sec^2]    

% ALTITUDE DIFFERENCE
Dh = alt - hc
if Dh < 0
    error('correcStdAtmos:valueOutOfRange', 'altitude is lower then ground level.')
end

% DATA
% stations
dTdh = -6.5; %[?K/km] temperature variation in troposphere

% TEMPERATURE MODEL
T = T0 + dTdh*Dh/1000;

% PRESSURE MODEL
p = p0*(1+dTdh/1000*Dh/T0).^(-g0/R/dTdh*1000);

% DENSITY MODEL
rho = p./R./T;

% SPEED OF SOUND
a = sqrt(gamma*R*T);
end