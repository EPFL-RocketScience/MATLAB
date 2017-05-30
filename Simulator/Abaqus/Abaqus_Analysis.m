function [ output_args ] = Abaqus_Analysis(epaisseur_bague,taille_maillage,epaisseur_CFRP,acceleration,R)
% Fonction qui rempli un formulaire pour cr?ation de fus?e dans abaqus avec
% mod?lisation des efforts en acc?l?ration verticale. 
% Abaqus_Analysis(epaisseur_bague,taille_maillage,epaisseur_CFRP,acceleration,R(Fus?e charg?e sous R));


%% Computation du tube de la fus?e
fid = fopen('./Abaqus/Param.py','wt'); % Le document est plac? dans le dossier courant
fprintf(fid,'# -*- coding: mbcs -*-');fprintf(fid, '\n'); % Abaqus demande cette commande
Lcone = R.Nose.L*1000; % N?cessite la longueur du cone pour le placement des objets
g = acceleration+1;
% Envoi de la taille du maillage
fprintf(fid, 'h = ');
fprintf(fid, num2str(taille_maillage));fprintf(fid, '\n');

% Envoi de l'?paisseur des plies de CFRP
fprintf(fid, 'plythickness = ');
fprintf(fid, num2str(epaisseur_CFRP));fprintf(fid, '\n');

% Envoi de la longueur du tube
Ltube = sum([R.Stage.L])*1000;
fprintf(fid, 'Ltube = ');
fprintf(fid, num2str(Ltube));fprintf(fid, '\n');

% Envoi du diam?tre du tube
Dtube = R.d*1000;
fprintf(fid, 'Dtube = ');
fprintf(fid, num2str(Dtube/2));fprintf(fid, '\n');

% Envoi de ebague
fprintf(fid, 'ebague = ');
fprintf(fid, num2str(-epaisseur_bague));fprintf(fid, '\n');

% Envoi toutes les pos de ref des objets par rapport au d?but du tube
posup = [[R.Cylinder.z]*1000-Lcone zeros(1,6-length([R.Cylinder.z]))]; % La longueur, [mm]
fprintf(fid,'posup = [') 
for i = 1:length(posup)
   fprintf(fid, num2str(-posup(i)));
   if(i<length(posup))
       fprintf(fid,',');
   end
end
fprintf(fid,']\n');

% Envoi toutes les pos d'arret des objets par rapport au d?but du tube
l = [R.Cylinder.L]*1000; % La longueur de tous les ?l?ment, [mm]
posdown = posup+[[R.Cylinder.L]*1000 zeros(1,6-length([R.Cylinder.L]))];
fprintf(fid,'posdown = [') 
for i = 1:length(posdown)
   fprintf(fid,num2str(-posdown(i)));
   if(i<length(posdown))
       fprintf(fid,',');
   end
end
fprintf(fid,']\n');

% Envoi de posuptank
posuptank = R.Motor.z*1000-Lcone;
fprintf(fid, 'posuptank = ');
fprintf(fid, num2str(-posuptank));fprintf(fid, '\n');

% Envoi de posdowntank
posdowntank = posuptank + R.Motor.L*1000;
fprintf(fid, 'posdowntank = ');
fprintf(fid, num2str(-posdowntank));fprintf(fid, '\n');

% Envoi des charges
mass_cyl = [[R.Cylinder.m] zeros(1,6-length([R.Cylinder.m]))]; % On veut les efforts au temps t=0, soit le d?collage, [kg]
A = 2*epaisseur_bague*pi*Dtube; % Surface ou s'applique chaque traction, [mm^2]
m_tank = R.Motor.m(0);
l_tank = R.Motor.L*1000;
load = [m_tank*9.81/(l_tank*pi*Dtube) 9.81/A.*mass_cyl]*g;
fprintf(fid,'load = [');
for i = 1:length(load)
   fprintf(fid,num2str(load(i)));
   if(i<length(load))
       fprintf(fid,',');
   end
end
fprintf(fid,']\n');


%% Computation des ailerons
%% R.fins(z, N, Ct, Cr, xt, S, r, e, rho);
%L_in = R.fins.Cr

fclose(fid);
%% Abaqus 
system('C:\SIMULIA\Abaqus\Commands\Abaqus cae script=Abaqus_Rocket.py')% Lance abaqus depuis la fen?tre de commande, le fichier doit ?tre dans le dossier courant.


end

