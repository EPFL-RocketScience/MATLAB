%** CATIA FileName
CADFILEPATH='./DESIGN_TABLE_ROCKET_V_2.txt';

CADFILE = importdata(CADFILEPATH);

%** parameter match
CONE = 1    ;  %cone(0) VS ogive(1) (!!NOT YET A AUTO NATCH)

LN  = CADFILE.data(14)  ;    %length of nose  
d  = CADFILE.data(9)   ; %diameter at base of nose  
LB = CADFILE.data(10)    ; %length of Body tube
dF  = d   ;   %diameter at front of transition  
dR  =  CADFILE.data(13)   ;  %diameter at rear of transition  
LT  = CADFILE.data(12)    ;%length of transition  
XP  =  CADFILE.data(10)+LN  ;   %distance from tip of nose to front of transition (=tube lenght + nose as transition at the end) 
CR  =  CADFILE.data(1)    ; %fins root chord  
CT  =  CADFILE.data(2) ;    %fins tip chord  
S  =  (CR-CT)/(1/tand(CADFILE.data(4))+1/tand(CADFILE.data(3)))   ;%fins semispan  
LF  =  (S^2+(CR/2-CT/2-S*tand(CADFILE.data(3)))^2)^.5  ;   %length of fin mid-chord line  
R  = d/2      ;%radius of body at aft end  
XR  = S/tand(CADFILE.data(3))     ; %distance between fin root leading edge and fin tip leading edge parallel to body  
XB  = XP-CR-CADFILE.data(8)     ;%distance from nose tip to fin root chord leading edge  
N  = CADFILE.data(7)      ;%number of fins  



%** center of mass

%* body
XCG_B = LB/2+LN;
%* nose cone
%* fins
%* empty engine
%* payload
%* chute




%** center of pressure (using Barrowman simplified equ.)
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
Xt=XP+(1+(1-dF/dR)/(1-(dF/dR)^2))*LT/3;

%Xt=0 ;% ??"NO CONICAL @ the moment

%* nose cone
CNn = 2;
if CONE == 1
Xn = 2/3*LN;
else %ogive case
Xn = 0.466*LN ;
end

%* fins
CNf=(1+R/(S+R))*((4*N*(S/d)^2)/(1+sqrt(1+(2*LF/(CR+CT))^2)));
Xf=XB+XR*(CR+2*CT)/3/(CR+CT)+1/6*((CR+CT)-(CR*CT)/(CR+CT));

%*total
 CNr = CNn + CNt + CNf ;
 Xp = (Xn*CNn+Xt*CNt+Xf*CNf)/CNr
 