%% Cas de charge
clear all;
close all;
clc;

%% Parametres

% General
L_tot;           % total length of the rocket
D_ext;           % diameter of the rocket
e;               % thickness of the tube wall
E_t;             % tube young's modulus
E_s;             % shock cord young's modulus
M_tot;           % Lift-off weight
M_bout;          % mass at burnout
S_s;             % shock cord section
S_l;             % shock cord length
rho;             % air density
Dt;              % half periode oscilation time after a pertubation

% Sizing
M1;              % Nose cone mass
L1;              % Nose cone length
D1;              % Nose cone top position (0)
M2;              % Chute + payload length
L2;              % Chute + payload length
D2;              % Chute top position
M3;              % Avionic mass
L3;              % Avionic length
D3;              % Avionic top position
M4;              % Tank mass (function of time)
L4;              % Tank length
D4;              % Tank top position
M5;              % Motor mass (function of time)
L5;              % Motor length
D5;              % Motor top position
Mt;              % Masse tube
Lt;              % Length of the tube
Dt;              % Tube top position


% Thrust
Thrust;          % Thrust, function of time (+20% uncertainty)
Mass;            % Total mass, function of time
A_max;           % Max acceleration due to thrust
v;               % Vitesse de la fusee

% Recovery
Cx;              % Chute drag coefficient (2.2 fruity chutes)
V_f;             % Desired speed at touchdown (5m/s)
V_op;            % Maximal speed with first parachute (IREC RULES)
V_im;            % Maximal speed at touchdown (IREC RULES)

%% Flambage d'Euler

% Nc = pi^2*E*I/l^2;                  % Charge critique fondamentale
% I = pi/64*(D_ext^4-(D_ext-e)^4);
% Ia = pi/16*(D_ext^3*e);             % Premiere approximation D>>e
Nc = (M1+M2+M3+M4+M5+Mt)*A_max;
e_lim = 16*l^2/(pi^3*D_ext^3*E)*Nc    % Epaisseur limite

%% Ouverture parachute

%approche energetique (toute l'energie cinetique perdue instantanement se dissipe 
%dans l'elongation de la corde : Ec = F*elongation)

Ec = 1/2*M_bout(V_op^2-V_f^2);
%el = F/S*E
Fm1 = (M_bout(V_op^2-V_f^2)*S_s*E_s/(2*S_l))^(1/2)


%approche drag (force de drag du parachute à Vmax)
S = 2*M_bout*9.81/(roh*Cx*V_f^2);  %Surface nocessaire pour obtenir 5m/s à l'atterissage
Fm2 = roh*S/2*Cx*V_op^2

%% Impact atterrissage
%Disspiation de toute l'energie cinetique lors du choc
Emax = 1/2*M_bout*V_im^2

% capacite d'absorbtion de la fibre de carbone par surface ?

%% Acceleration laterale
% Attention toutes les longueures doivent être en metre
% hypothese d'une trajectoire circulaire pour la sequence de stabilisation
% apres une perurbation qui induit un angle d'attaque alpha

Rayon = Dt*v/(2*alpha);
Ac = v^2/Rayon;                 %acceletarion laterale

m = [M1 M2 M3 M4 M5 Mt];
L = [L1 L2 L3 L4 L5 Lt];
d = [D1 D2 D3 D4 D5 Dt];
if Ac>3
    P = m*Ac./L;                %forces reparties
else
   g=9.81; 
   Acn = 3*g;                   %si acceleration laterale trop faible on prend 3g pour avoir une valeur de base pour le transport
   P = m*Acn./L;                    
end
el = 0.01;                      %increment de position (taille du maillage en m)

% analyse separee en 3 parties pour chaque force repartie ensuite appliquer
% principe de superposition.
% Calcul des reactions Ay, By, calcul de T(x) et M(x) le long de la fusee en
% fonction de la position de la force repartie
xx = 0:el:L_tot;
lm = length(m);                 %permet de rendre le code independant du nombre d'element
lx = length(xx);
M_tot = zeros(lm,lx);
T_tot = zeros(lm,lx);

for i=1:lm
    Ay = L(i)*P(i)*(1-(d(i)+L(i)/2)/L_tot);
    By = L(i)*P(i)*(d(i)+L(i)/2)/L_tot;
    T1 = []; T2 = []; T3 = [];
    Mo1 =[]; Mo2 = []; Mo3 = [];

    for x=xx
        % x<d
        if x<=d(i)
        T1 = [T1 Ay];
        Mo1 = [Mo1 Ay*x];
        
        elseif x>d(i) && x<=d(i)+L(i)
        %pour d<x<d+L
        T2 = [T2 Ay-P(i)*(x-d(i))];
        Mo2 = [Mo2 Ay*x-P(i)*(x-d(i))*(x-d(i))/2];
        
        else
        %pour x>d+L
        T3 = [T3 Ay-P(i)*L(i)];
        Mo3 = [Mo3 Ay*x-P(i)*L(i)*(x-d(i)-L(i)/2)];
        end
    end

    T_tot(i,:) = [T1 T2 T3];
    M_tot(i,:) = [Mo1 Mo2 Mo3];      %valeur des moments pour chaque force

end

Mfin = sum(M_tot);
Tfin = sum(T_tot);

figure(1)
plot(xx,Tfin,'b','DisplayName','T(x)')
hold on
plot(xx,Mfin,'r','DisplayName','M(x)')
hold off
legend show

%% Calcul des contraintes de traction maximales
I = pi/64*(D_ext^4-(D_ext-e).^4);
S = pi*D_ext*e;                     % surface de coupe
ymax = D_ext/2;
ld=length(Mfin);
le=length(e);                       % possible de passer un veteur d'epaisseur

sigma=zeros(le,ld);
tau = zeros(le,ld);
for i=1:length(e)
    sigma(i,:) = ymax*Mfin/I(i);    % sigma le long de la fibre la plus extérieure
    tau(i,:) = 2*Tfin/S(i);         % valeur de tau max pour un tube fin
end
figure(2)
plot(xx,sigma(1,:),'b','LineWidth',1,'DisplayName','\sigma_{ymax} (x) e=1')
hold on
grid on
plot(xx,sigma(2,:),'k','LineWidth',1,'DisplayName','\sigma_{ymax} (x) e=2')
plot(xx,sigma(3,:),'g','LineWidth',1,'DisplayName','\sigma_{ymax} (x) e=3')
plot(xx,tau(1,:),'r','LineWidth',1,'DisplayName','\tau_{max} (x)')
hold off
legend show

sigma_max = max(max(sigma))         %contrainte maximale en traction