function drawRocket(R)
% dessine la fusee

    figure; hold on;

    % cone
    triplot([1 2 3], [0 R.Nose.L R.Nose.L], [0 R.Nose.D/2 -R.Nose.D/2], 'k');

    % stage
    for stage = R.Stage
        rectangle('Position',[stage.z,-stage.Dout/2,stage.L,stage.Dout], 'EdgeColor', 'k');
    end

    % tail
    plot([R.Tail.z R.Tail.z R.Tail.z+R.Tail.L R.Tail.z+R.Tail.L R.Tail.z],...
         [-R.Tail.D1/2 R.Tail.D1/2 R.Tail.D2/2 -R.Tail.D2/2 -R.Tail.D1/2], 'k');

    % fins
    plot([R.Fins.z R.Fins.z+R.Fins.xt R.Fins.z+R.Fins.xt+R.Fins.Ct R.Fins.z+R.Fins.Cr R.Fins.z],...
         [R.Fins.r R.Fins.r+R.Fins.S  R.Fins.r+R.Fins.S            R.Fins.r           R.Fins.r], 'k');
    plot([R.Fins.z R.Fins.z+R.Fins.xt R.Fins.z+R.Fins.xt+R.Fins.Ct R.Fins.z+R.Fins.Cr R.Fins.z],...
         -1*[R.Fins.r R.Fins.r+R.Fins.S  R.Fins.r+R.Fins.S            R.Fins.r           R.Fins.r], 'k');
    
    [CN, CP] = R.aeroCoeff(0,0,0);
    % cp
    plot(CP, 0, 'ro');
    % cm 
    plot(R.cm(0), 0, '*b');
     
    daspect([1 1 1]);
    set(gca,'FontSize',30);
    set(gcf, 'Position', [0, 0, 3840, 2160]);
    saveas(gcf,'rocket','eps')

end

