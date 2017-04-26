function R = TestRocket()

% parameters
rho_carbon = 1550; % kg/m^3
data = load('HyperTEK_Mgrain.mat');
thrustCurve = data.Tdat;
bt = 10; % burn time


% Create a rocket
R = Rocket();

% nose(L, D, e, rho)
R.nose(0.5, 0.15, 0.002, rho_carbon);
% stage(id, z, L, Dout, e, rho)
R.stage('tube', 0.5, 2.5, 0.15, 0.002, rho_carbon);
% tail(D1, D2, L, e, z, rho)
R.tail(0.15, 0.12, 0.08, 0.002, 3, rho_carbon);
% motor(id, z, D, L, e, m, mp, rho, thrustCurve, bt)
R.motor('motor', 2.5, 0.098, 0.75, 1e-3, 9.037, 4.237, 1772, thrustCurve, bt);
% payload(id, z, m, L, D)
R.cylinder('main payload', 0.5, 8, 0.3, 0.15);
% parachute(id, z, L, m, D, Cd)
R.parachute('main', 0.3, 0.1, 1, 2, 1.5);
R.parachute('secondary', 0.5, 0.1, 0.5, 0.5, 0.75);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(2.65, 3, 0.2, 0.35, 0.15, 0.15, 0.075, 0.005, rho_carbon);
% point(id, z, L, m)
R.point('tuyere', 2.9, 0.2, 1);

% Check definition
R.check();

% Print specs
R.printSpecs();

end

