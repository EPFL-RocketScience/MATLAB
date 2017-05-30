function R = RocketJuju()

% parameters
rho_polystyr    = 1050; % kg/m^3
rho_kraft       = 950;
rho_birch       = 630;

% Motor
data = load('AeroTech_H148R.mat');
thrustCurve = cleanThrustCurve(data.data);
bt = 1.5; % burn time

% Create a rocket
R = Rocket();
% nose(L, D, e, rho)
R.nose(0.42, 0.102, 0.0035, rho_polystyr);
% stage(id, z, L, Dout, e, rho)
R.stage('tube', 0.42, 1.24, 0.102, 0.0035, rho_kraft);
% tail(D1, D2, L, e, z, rho)
R.tail(0.102, 0.101, 0.01, 0.0035, 1.66, 0);
% motor(id, z, D, L, e, m, mp, rho, thrustCurve, bt)
R.motor('motor', 1.47, 0.038, 0.17, 1.5e-3, 0.309, 0.148, 1772, thrustCurve, bt);
% payload(id, z, m, L, D)
R.cylinder('main payload', 0.52, 0.4, 0.3, 0.102);
% parachute(id, z, L, m, D, Cd)
R.parachute('main', 1.12, 0.14, 0.050, 1, 1.5);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(1.35, 3, 0.086, 0.3, 0.166, 0.096, 0.051, 0.008, rho_birch);
% point(id, z, L, m)
R.point('trim_nose', 0.3, 0.1, 0.3);
% coupler(obj, id, z, L, Dout, e, rho)
R.coupler('coupler', 0.81, 0.2, 0.098, 0.003, rho_kraft);
R.coupler('nose shoulder', 0.42, 0.1, 0.098, 0.003, rho_polystyr)
R.coupler('motor tube', 1.26, 0.4, 0.042, 0.003, rho_kraft);

% Additional parameters:
% lifting bodys correction factor
R.k = 1.1; % default 
% Drag coefficient
R.CD = 0.55;

% Check definition
R.check();

% Print specs
R.printSpecs();

end

