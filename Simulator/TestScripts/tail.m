function Tail = tail(D1, D2, L, e, rho)
    % INPUTS
    % D1 :  diametre du haut
    % D2 :  diametre du bas
    % L :   longueur du tail
    % e :   epaisseur de paroie
    % z :   postion par rapport au haut de la fus?e du haut du tail
    % rho : densite

    %Calcule intermediaire
    hF = L*D2/(D1-D2);
    hTot = L*(D1/(D1-D2));
    vTot = e*pi*sqrt(1+D1^2/(4*L^2))*hTot*D1/2;
    vF = e*pi*sqrt(1+D2^2/(4*L^2))*hF*D2/2;
    Ixprime = rho/4*(vTot*(D1^2/4+2*hTot^2)-vF*(D2^2/4+2*hF^2));
    hPrime = hF+L*((D2+2*D1)/(3*(D2+D1)));

    % Calcule des proprietes de masse
    %on utilise un cone - un plus petit cone
    Tail.m = rho*(vTot-vF);
    Tail.cm = L*(2*D2+D1)/(3*(D2+D1));
    Tail.Iz = rho/8*(vTot*D1^2-vF*D2^2);
    Tail.Ir = Ixprime - (Tail.m)*(hPrime^2-(Tail.cm)^2);            

    % Calcule des proprietes aerodynamiques
    Tail.CN = 0; % coefficient aerodynamique normal
    Tail.zCP = 0; % position du centre de pression relatif

end
