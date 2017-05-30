function R = RocketIREC()

% parameters
rho_carbon      = 1780; 
rho_aluminum    = 2700;
rho_SRF         = 1772;

% Motor
data = load('Aerotech_N2220.mat');
thrustCurve = cleanThrustCurve(data.data);
bt = 5.4; % burn time

% Create a rocket
R = Rocket();
% nose(L, D, e, rho)
R.nose(0.5, 0.17, 0.0013, rho_carbon);
% stage(id, z, L, Dout, e, rho)
R.stage('upper_tube', 0.5, 1.2, 0.17, 0.001, rho_carbon);
R.stage('motor_tube', 1.7, 1.3, 0.17, 0.001, rho_carbon);
% tail(D1, D2, L, e, z, rho)
R.tail(0.17,0.12,0.08,0.0015,3,rho_aluminum);
% motor(id, z, D, L, e, m, mp, rho, thrustCurve, bt)
R.motor('motor', 2.08, 0.098, 1, 0.005, 11.997, 7.183, rho_SRF, thrustCurve, bt);
% cylinder(id, z, m, L, D)
R.cylinder('Payload_main', 1.02, 5, 0.3, 0.15);
R.cylinder('Payload_secondary', 0.3, 1, 0.1,0.1);
R.cylinder('Recovery_system', 1.32, 1, 0.15, 0.13);
R.cylinder('Avionics_main', 1.47, 0.5, 0.2, 0.13);
R.cylinder('Avionics_control', 1.72, 0.5, 0.15, 0.13);
% parachute(id, z, L, m, D, Cd)
R.parachute('Parachute_main', 0.63, 0.4, 1 ,2.7, 1.5);
R.parachute('Parachute_drogue', 0.584, 0.04, 0.062, 1, 1.5);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(2.65, 3, 0.15, 0.35, 0.15, 0.17, 0.085, 0.005, rho_carbon);
% point(id, z, L, m)
R.point('GNSS', 0.4, 0.05, 0.1);
R.point('Altimeter_Bkp', 0.45, 0.05, 0.1);
R.point('Aerobreaks', 3, 0.2, 0.2);
% coupler(id, z, L, Dout, e, rho)
R.coupler('Coupler-upper', 0.5, 0.115, 0.168, 0.0055, rho_aluminum);
R.coupler('Coupler-lower', 1.54, 0.115, 0.168, 0.0055, rho_aluminum);
R.coupler('Motor-mount', 2.03, 1.05, 0.102, 0.001, rho_carbon);

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

