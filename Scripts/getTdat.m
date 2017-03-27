clear all
close all

% Code to get data points on thrust curve image
filename = 'Thrust_Curves/Contrail_O.png';
im = imread(filename); 
figure;

haxim = axes;    % haxim est le 'handle' qui permet d'acceder aux parametres des axes crees dans la figure
hvid = image(im); % hvid est le 'handle' qui permet d'acceder aux donnes de l'image affiche

% L'axe vertical de l'image affiche par matlab croit vers le bas. Ici on le change
yplot=get(hvid,'YData');   % on recupere les coordones 'y' des pixels affiches
set(hvid,'YData',yplot(end:-1:1)); % on inverse ces coordones 
set(haxim,'YDir','normal'); % on affiche l'axe verticale croissant vers le haut
hold on  

% Ajustement des axes pour avoir les vraies unit?s
disp('Recuperez le coin bas gauche puis le coin haut droit');
bords=ginput(2);  % on recupere les bords des axes (X0,Y0), (Xmax,Ymax) 
X0 = bords(1,1); Y0 = bords(1,2);
Xmax=bords(2,1); Ymax=bords(2,2);
totY=Ymax-Y0;     % l'echelle de hauteur 
totX=Xmax-X0;     % l'echelle de temps (totale)
yplot = get(hvid,'YData'); % vecteur des coordones verticales (pixels)
xplot = get(hvid,'XData'); % vecteur des coordones horizontales (pixels)

% Dans quels uunites travaille-t-on?
unit = input('Quels sont les unites de la courbe de poussee (N=0/lb=1)?')
if unit<0 || unit>1
   error('getTdat:outOfRange', 'entered value for units is out of range. Allowed values are N = 0 or lb = 1'); 
end

% avant l'ajustement, on demande ? l'utilisateur l'?chelle des axes
disp('Entrez les valeurs extr?mes des axes pour connaitre les unit?s')

tmin = input('tmin:');
tmax = input('tmax:');
Nmin = input('Tmin:');
Nmax = input('Tmax:'); 

Zax = (yplot-Y0)/totY*(Nmax-Nmin);     
Tax = (xplot-X0)/totX*(tmax-tmin);

set(hvid,'XData',Tax,'YData',Zax); % on met les coordonees physiques
set(haxim,'Xlim',[tmin tmax],'Ylim',[Nmin Nmax]); % on etablit les limites d'affichage correspondants

% L'utilisateur r?cup?re des points de la courbe de pousee
disp('Recuperez maintenant les points de mesure puis appuyez sur entree');
Tdat=ginput; % on recupere les donnees

% correction is on travaille en lb
if unit == 1
    Tdat =[Tdat(:,1), Tdat(:,2)*4.45];
end

% Les donnees sont sauv?es dans un fichier:
save(filename(1:end-4), 'Tdat');