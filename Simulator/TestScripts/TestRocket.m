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
% motor(id, z, m, D, L, thrustCurve, bt)
m = @(t) (10 - 0.4*t).*(t<=bt)+6.*(t>bt);
R.motor('grain', 1.65, m, 0.098, 1.43, thrustCurve, bt);
% payload(id, z, m, L, D)
R.payload('main', 0.5, 5, 0.3, 0.15);
% parachute(id, z, m, D, Cd)
R.parachute('main', 0.3, 1, 2, 1.5);
R.parachute('secondary', 0.5, 0.5, 0.5, 0.75);
% fins(z, N, Ct, Cr, xt, S, r, e, rho);
R.fins(2.65, 3, 0.2, 0.35, 0.15, 0.15, 0.075, 0.005, rho_carbon);
% point(id, z, m)
R.point('tuyere', 2.9, 1);

% Check definition
R.check();

% Print specs
R.printSpecs();

end

