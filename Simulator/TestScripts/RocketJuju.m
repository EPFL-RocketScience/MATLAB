function R = TestRocket()

% parameters
rho_polystyr    = 1050; % kg/m^3
rho_kraft       = 950;
rho_birch       = 630;

data = load('AeroTech_H123W.mat');
thrustCurve = cleanThrustCurve(data.data);
bt = 2.37; % burn time


% Create a rocket
R = Rocket();
% nose(L, D, e, rho)
R.nose(0.42, 0.102, 0.0022, rho_polystyr);
% stage(id, z, L, Dout, e, rho)
R.stage('tube', 0.42, 1.24, 0.102, 0.003, rho_kraft);
% tail(D1, D2, L, e, z, rho)
R.tail(0.102, 0.101, 0.001, 0.002, 1.66, rho_kraft);
% motor(id, z, D, L, e, m, mp, rho, thrustCurve, bt)
R.motor('motor', 1.47, 0.038, 0.17, 1.5e-3, 0.309, 0.148, 1772, thrustCurve, bt);
% payload(id, z, m, L, D)
R.cylinder('main payload', 0.52, 0.4, 0.3, 0.102);
% parachute(id, z, L, m, D, Cd)
R.parachute('main', 1.12, 0.14, 0.005, 1, 1.5);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(1.35, 3, 0.086, 0.3, 0.166, 0.096, 0.051, 0.0065, rho_birch);
% point(id, z, L, m)
R.point('coupler', 0.9, 0.4, 0.172);
R.point('trim_nose', 0.3, 0.3, 0.5);
% Check definition
R.check();

% Print specs
R.printSpecs();

end

