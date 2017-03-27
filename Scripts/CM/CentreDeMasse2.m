clear all;
close all;

%** CATIA FileName
CADFILEPATH='C:\Users\Jonathan\Documents\GrabCAD\ROCKETS\MODEL_CAD_Emilien\DESIGN_TABLE_JULIEN\DESIGN_TABLE_ROCKET_V_2.txt';

CADFILE = importdata(CADFILEPATH);

%% ** parameter match
CONE = 0    ;  %cone(0) VS ogive(1) (!!NOT YET A AUTO NATCH)

LN  = CADFILE.data(13)  ;    %length of nose  
d  = CADFILE.data(8)   ; %diameter at base of nose  
dF  = d   ;   %diameter at front of transition  
dR  =  CADFILE.data(12)   ;  %diameter at rear of transition  
LT  = CADFILE.data(11)    ;%length of transition  
XP  =  CADFILE.data(9)+LN  ;   %distance from tip of nose to front of transition (=tube lenght + nose as transition at the end) 
CR  =  CADFILE.data(1)    ; %fins root chord  
CT  =  CADFILE.data(2) ;    %fins tip chord  
S  =  (CR-CT)/(1/tand(CADFILE.data(4))+1/tand(CADFILE.data(3)))   ;%fins semispan  
LF  =  (S^2+(CR/2-CT/2-S*tand(CADFILE.data(3)))^2)^.5  ;   %length of fin mid-chord line  
R  = d/2      ;%radius of body at aft end  
XR  = S/tand(CADFILE.data(3))     ; %distance between fin root leading edge and fin tip leading edge parallel to body  
XB  = XP-CR-CADFILE.data(7)     ;%distance from nose tip to fin root chord leading edge  
N  = CADFILE.data(6)      ;%number of fins  
E = CADFILE.data(10); %tube_thickness (mm)
E_fins = CADFILE.data(5); % Fins thickness (mm)
TL = CADFILE.data(9); %longueur du tube
gamma = CADFILE.data(3)/360*2*pi;%fins_front_angle [rad]
phi = CADFILE.data(4)/360*2*pi;%fins_back_angle [rad]


%% ** center of mass
% INPUT from datasheet
b = CR;
a = CT;

%INPUT
rho = 2.7*10^(-6); % density of the rocket in kg/mm^3 Attention if fins different material !!!
rho_fins = rho;
rho_tank = 4*10^(-6); %density of tank [kg/mm^3]
h_pos_tank = 800; % position par rapport au haut du réservoir
h_length_tank = 400; % taille du réservoir
typeTank = 1; % 1 => cylinder, 2 => cylinder+sphère
d_tank = 80; %diameter of the tank
e_tank = 2; %thickness of the tank

%* nose cone
z_cm_cone = 3*LN/4; % poisition of the center of mass in [mm]
m_cone = rho*pi/3*(((d/2)^2*LN)-(((d/2)-E)^2*(LN-E))); %mass of nose cone

%* fins
h = (b - a)/(1/tan(gamma)+1/tan(phi)); %high fins
z_cm_fins = XB + (h/3)*(2*a+b)/(a+b);
m_fins = rho_fins * N * (1/2)*(a+b)*h*E_fins; % mass for N fins

%* empty engine
z_cm_tank = h_pos_tank + h_length_tank/2;
%compute mass of the empty tank
if typeTank == 1
    m_tank = rho_tank*pi*(((d_tank/2)^2*h_length_tank)-((d_tank/2-e_tank)^2*(h_length_tank-2*e_tank)));
else
    m_tank = rho_tank*pi*(((d_tank/2)^2*(4/3*(d_tank/2)+h_length_tank-d_tank))-(((d_tank/2-e_tank)^2*(4/3*(d_tank/2-e_tank)+h_length_tank-d_tank-2*e_tank))));
end

%* empty tube
z_cm_tube = LN + TL/2;
m_tube = pi*rho*((d/2)^2-(d/2-E)^2)*TL;

%* payload
m_payload = 4; %in [kg]
z_cm_payload = 300; % Same like position of the payload

%* chute
m_chute = 2.5; %in [kg]
z_cm_chute = 350; % Same like position of the payload

% Compute of the mass and the center mass
m_tot_empty = m_payload + m_tube + m_chute + m_tank + m_cone + m_fins;
z_cm_empty = (z_cm_fins*m_fins + z_cm_cone*m_cone + z_cm_payload * m_payload + z_cm_chute * m_chute + z_cm_tank * m_tank + z_cm_tube * m_tube)/(m_tot_empty);

%variation of the mass in flight
t_stop = 10;%temps de vol
m_i = 10; % mass inital en [kg]
m_f = 0; % mass final [kg]

for t = 0:t_stop
    m_moteur(t+1) = mass( t, t_stop, m_i, m_f );
    m_vol(t+1) = m_tot_empty + m_moteur(t+1);
    z_cm(t+1) = (m_moteur(t+1)* z_cm_tank + z_cm_fins*m_fins + z_cm_cone*m_cone + z_cm_payload * m_payload + z_cm_chute * m_chute + z_cm_tank * m_tank + z_cm_tube * m_tube)/(m_vol(t+1));
   
end
t = [0:1:t_stop];
plot(t,z_cm);


%% ** center of pressure (using Barrowman simplified equ.)
%variables
 Xn=0  ; %nose cone pressure point hight(zero @ top of rocket)
 Xt=0   ;%conical end pressure point hight
 Xf=0   ;%fins pressure point hight
 Xp=0   ;%rocket pressure point hight
 CNn=0  ;%nose cone
 CNt=0  ;%conical end
 CNf=0  ;%fins
 CNr=0 ; %rocket 

%* body conical
CNt=2*((dR/d)^2-(dF/d)^2);
Xt=Xp+(1+(1-dF/dR)/(1-(dF/dR)^2))*LT/3;

Xt=0 ;% ¨¨"NO CONICAL @ the moment

%* nose cone
CNn = 2;
if CONE == 1
Xn = 0.666*LN;
else %ogive case
Xn = 0.466*LN ;
end

%* fins
CNf=(1+R/(S+R))*((4*N*(S/d)^2)/(1+sqrt(1+(2*LF/(CR+CT))^2)));
Xf=XB+XR*(CR+2*CT)/3/(CR+CT)+1/6*((CR+CT)-(CR*CT)/(CR*CT));

%*total
 CNr = CNn + CNt + CNf ;
 Xp = (Xn*CNn+Xt*CNt+Xf*CNf)/CNr;
 
 
 %% Marge statique à vide
 marge = (Xp-z_cm_empty)/(d)
 marge_time = (Xp-z_cm)/(d);
 plot(t,marge_time);