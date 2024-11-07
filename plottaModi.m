% file di esempio input
clear
close all

%% modi di vibrare
% N.B.
% - in questo set di modi i torsionali sono il 11 (1° torsionale)
%   e il 14 (2° torsionale) 

load('struct_1.mat')
for nmodo = 1:size(fisez,2)
    
    phi_y = squeeze(fisez(1,nmodo,:)); % componente y del modo nmodo
    phi_z = squeeze(fisez(2,nmodo,:)); % componente z del modo nmodo
    phi_t = squeeze(fisez(3,nmodo,:)); % componente theta del modo nmodo --> in RADIANTI
    
    figure
    hold on
    plot(xyzdeck(:,1),phi_y)
    plot(xyzdeck(:,1),phi_z)
    plot(xyzdeck(:,1),phi_t)
    title(sprintf('mode %02i, f = %6.4f Hz, modal mass  = %7.1f, damping ratio = %5.4f ',nmodo,parmod.freq(nmodo),parmod.m(nmodo),parmod.h))
    legend('\phi_y','\phi_z','\phi_\theta')
    xlabel('Deck axis [m]')
    ylabel('\phi')
    grid on
    set(gca,'box','on')
end




