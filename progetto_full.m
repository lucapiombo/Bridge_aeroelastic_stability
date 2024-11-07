clear
close all
clc

%% RISULTATO
% va in Flutter a 61 [m/s], modale full

tic

load('struct_1.mat')
L_sez = L/243;
B1y=0*B; 
B1z=0*B;
B1t=0.3*B;

%inizializzo matrici diagonali
M_diag = diag(parmod.m);
r = parmod.h*2.*parmod.m.*parmod.freq*2*pi; %r = h*2*m*w, dove w = 2*pi*f
R_diag = diag(r);
K_diag = diag(parmod.k);

%ricavo matrice autovettori PHI da fisez ((3*243)(gdl))*14(modi)
PHI = [];
for i = 1:1:243
    PHI = [ PHI; fisez(:,:,i)];
end

%% STATICA
%Aerodynamic coefficients, li ricavo da file esterno.mat
load('STA_1');   %alpha is in [deg]!

Cd=DLM(:,1);
Cl=DLM(:,2);
Cm=DLM(:,3);
Cd_alpha=KDLM(:,1);
Cl_alpha=KDLM(:,2);
Cm_alpha=KDLM(:,3);
rho = 1.225;

vel_min = 0;  vel_max = 80;  nvel = 81;
vel = linspace(vel_min,vel_max,nvel);  %wind speed vector 

% nonlinear solution
sol_old=zeros(14,1);%guess iniziale, non è a caso 0, infatti il vettore velocità  parte da zero
q0_statico=zeros(14,length(vel));%inizializzo

hh = waitbar(0,'Computing static displacements, full bridge...'); 
steps = length(vel);

for ii=1:length(vel)
    U = vel(ii);
    
	% fsolve find the nonlinear solution of the problem defined in the function 'statica', that is in another m-file
    [ampiezze_modali] = fsolve(@(q0) statica_modale_full (K_diag,q0,rho,U,B,L_sez,Cd,Cl,Cm,alpha, PHI),sol_old);%
    q0_statico(: , ii)= ampiezze_modali;
    sol_old=ampiezze_modali;       
    waitbar(ii/steps)
    
end
close(hh)
%%

%una volta ricavate le ampiezze dei modi ÃƒÂ¨ un attimo ricavare gli
%spostamenti fisici
y_statico = zeros(ii, 243);
z_statico = zeros(ii, 243);
theta_statico = zeros(ii, 243);
sezioni = [1:1:243];
%%

for ii=1:length(vel)
    x = PHI*q0_statico(:, ii);
        for k=1:243
            y_statico(ii,k) = [x(1+3*(k-1))];
            z_statico(ii,k) = [x(2+3*(k-1))];
            theta_statico(ii,k) = [x(3+3*(k-1))];
        end 
%     plot(sezioni, y_statico(ii,:));
%     hold on
%     plot(sezioni, z_statico(ii,:));
%     plot(sezioni, rad2deg(theta_statico(ii,:)));
%     grid on
%     grid minor
%     legend('y statico','z statico','theta statico')
%     xlabel('sezioni')
%     ylabel('spostamenti: [m]')
%     pause(0.05)
end

%% plot singole variabili [y, z, theta]

figure(101)
    plot(sezioni, y_statico(10,:), 'Linewidth', 2);
    hold on
    plot(sezioni, y_statico(35,:), 'Linewidth', 2);
    plot(sezioni, y_statico(60,:), 'Linewidth', 2);
    grid on
    grid minor
    legend('10 [m/s]','35 [m/s]','60 [m/s]','Location', 'best')
    xlabel('sections')
    ylabel('y dispacement [m]')
    pause(0.05)
    saveas(gcf, 'Immagini_risultati/full_bridge/y_statico.png')

figure(102)
    plot(sezioni, z_statico(10,:), 'Linewidth', 2);
    hold 
    plot(sezioni, z_statico(35,:), 'Linewidth', 2);
    plot(sezioni, z_statico(60,:), 'Linewidth', 2);
    grid on
    grid minor
    legend('10 [m/s]','35 [m/s]','60 [m/s]','Location', 'best')
    xlabel('sections')
    ylabel('z displacement [m]')
    pause(0.05)
    saveas(gcf, 'Immagini_risultati/full_bridge/z_statico.png')

figure(103)
    plot(sezioni, rad2deg(theta_statico(10,:)), 'Linewidth', 2);
    hold on
    plot(sezioni, rad2deg(theta_statico(35,:)), 'Linewidth', 2);
    plot(sezioni, rad2deg(theta_statico(60,:)), 'Linewidth', 2);
    grid on
    grid minor
    legend('10 [m/s]','35 [m/s]','60 [m/s]','Location', 'best')
    xlabel('sections')
    ylabel('theta rotation [deg]')
    saveas(gcf, 'Immagini_risultati/full_bridge/theta_statico.png')
    pause(0.05)


%% Dinamica

hh = waitbar(0,'Calcolo velocità critica di flutter FULL BRIDGE...');
steps = length(vel);

for ii = 1:1:length(vel)

    U = vel(ii);
    t_eq_rad = theta_statico(ii,:);
    t_eq_deg = rad2deg(theta_statico(ii,:)); 
    
    for k = 1:1:243
        Cd_eq(k) = interp1(alpha, Cd, t_eq_deg(k));
        Cl_eq(k) = interp1(alpha, Cl, t_eq_deg(k));
        Cm_eq(k) = interp1(alpha, Cm, t_eq_deg(k));
        Cd_alpha_eq(k) = interp1(alpha, Cd_alpha, t_eq_deg(k));
        Cl_alpha_eq(k) = interp1(alpha, Cl_alpha, t_eq_deg(k));
        Cm_alpha_eq(k) = interp1(alpha, Cm_alpha, t_eq_deg(k));
    end
    
    %costruisco matrici aerodinamiche
    qBL_sez = 0.5*rho*U^2*B*L_sez;
    
    %RIGIDEZZA
    z0 = zeros(1, 243);
    K_aero_temp_scazzata = qBL_sez*[ z0, z0, -Cd_alpha_eq;
                                     z0, z0, -Cl_alpha_eq;
                                     z0, z0, -B*Cm_alpha_eq];
                                 
        K_aero_temp = [];
    for k = 1:1:243
        vett_K = K_aero_temp_scazzata (:,243*2+k);
        matrice_K = [zeros(3,1), zeros(3,1), vett_K]; 
        K_aero_temp = [K_aero_temp, matrice_K];
    end
    
    %SMORZAMENTO
    
    R_aero_temp = [];
    
    for k = 1:1:243
        R_aero_temp_sezione = L_sez * [  0.5*rho*U*B*2*Cd_eq(k), 0.5*rho*U*B*(Cd_alpha_eq(k)-Cl_eq(k)), 0.5*rho*U*B*(Cd_alpha_eq(k)-Cl_eq(k))*B1y;
                                         0.5*rho*U*B*2*Cl_eq(k), 0.5*rho*U*B*(Cl_alpha_eq(k)+Cd_eq(k)), 0.5*rho*U*B*(Cl_alpha_eq(k)+Cd_eq(k))*B1z;
                                         0.5*rho*U*B*B*2*Cm_eq(k),   0.5*rho*U*B*B*Cm_alpha_eq(k),    0.5*rho*U*B*B*Cm_alpha_eq(k)*B1t              ];
        R_aero_temp = [R_aero_temp, R_aero_temp_sezione];
    end
    
    %abbiamo ricavato K_aero e R_aero per tutte le sezioni alla velocità ii [sono le  _temp]
    
    
    
    K_aero_modale_sum_matrice = [];
    R_aero_modale_sum_matrice = [];
    
    for n = 1:1:14
        
        for j = 1:1:14
            
                K_aero_modale_sum = 0;
                R_aero_modale_sum = 0;
                
            for k = 1:3:729
                
                phi_j = PHI(k:k+2, j);
                phi_n_trasposto = PHI(k:k+2, n)';
                
                K_aero_k = K_aero_temp(:, k:k+2);
                R_aero_k = R_aero_temp(:, k:k+2);
                
                K_aero_modale_k = phi_n_trasposto*K_aero_k*phi_j;
                R_aero_modale_k = phi_n_trasposto*R_aero_k*phi_j;
                
                K_aero_modale_sum = K_aero_modale_sum + K_aero_modale_k;
                R_aero_modale_sum = R_aero_modale_sum + R_aero_modale_k;
                
            end
            
            K_aero_modale_sum_vettore(j) = K_aero_modale_sum;
            R_aero_modale_sum_vettore(j) = R_aero_modale_sum;
            
        end
        
        K_aero_modale_sum_matrice = [ K_aero_modale_sum_matrice; K_aero_modale_sum_vettore];
        R_aero_modale_sum_matrice = [ R_aero_modale_sum_matrice; R_aero_modale_sum_vettore];
        
    end

M_din = M_diag;
R_din = R_diag + R_aero_modale_sum_matrice;
K_din = K_diag + K_aero_modale_sum_matrice;


[autovet, autoval] = polyeig(K_din, R_din, M_din);

%il riordino magico

                if ii==1
                    f_old = parmod.freq; 
                    smo_old = r;
                    lambda_old(1:14) = -smo_old.* (2*pi*f_old) + 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
                    lambda_old(15:28) = -smo_old.* (2*pi*f_old) - 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
                end
                
                indici = 1:28;
                for jj = 1:28
                    [val,dove]=min(abs(autoval(indici)-lambda_old(jj)));
                    ind(jj) = indici(dove);
                    indici = setdiff(indici,ind(jj));
                end
                autoval = autoval(ind);
                autovet = autovet(:,ind);
     
        smo_1(ii)= -real(autoval(1))/norm(autoval(1));
        smo_2(ii)= -real(autoval(3))/norm(autoval(3));
        smo_3(ii)= -real(autoval(5))/norm(autoval(5));
        smo_4(ii)= -real(autoval(7))/norm(autoval(7));
        smo_5(ii)= -real(autoval(9))/norm(autoval(9));
        smo_6(ii)= -real(autoval(11))/norm(autoval(11));
        smo_7(ii)= -real(autoval(13))/norm(autoval(13));
        smo_8(ii)= -real(autoval(15))/norm(autoval(15));
        smo_9(ii)= -real(autoval(17))/norm(autoval(17));
        smo_10(ii)= -real(autoval(19))/norm(autoval(19));
        smo_11(ii)= -real(autoval(21))/norm(autoval(21));
        smo_12(ii)= -real(autoval(23))/norm(autoval(23));
        smo_13(ii)= -real(autoval(25))/norm(autoval(25));
        smo_14(ii)= -real(autoval(27))/norm(autoval(27));
        
        
        freq_1(ii)= abs(imag(autoval(1)))/(2*pi);
        freq_2(ii)= abs(imag(autoval(3)))/(2*pi);
        freq_3(ii)= abs(imag(autoval(5)))/(2*pi);
        freq_4(ii)= abs(imag(autoval(7)))/(2*pi);
        freq_5(ii)= abs(imag(autoval(9)))/(2*pi);
        freq_6(ii)= abs(imag(autoval(11)))/(2*pi);
        freq_7(ii)= abs(imag(autoval(13)))/(2*pi);
        freq_8(ii)= abs(imag(autoval(15)))/(2*pi);
        freq_9(ii)= abs(imag(autoval(17)))/(2*pi);
        freq_10(ii)= abs(imag(autoval(19)))/(2*pi);
        freq_11(ii)= abs(imag(autoval(21)))/(2*pi);
        freq_12(ii)= abs(imag(autoval(23)))/(2*pi);
        freq_13(ii)= abs(imag(autoval(25)))/(2*pi);
        freq_14(ii)= abs(imag(autoval(27)))/(2*pi);
        
        autov(:,ii) = autoval;
        lambda_old  = autoval;
        
waitbar(ii/steps)

end

close(hh)

%% Plot finali

figure()
hold on; grid on
plot(vel,freq_1,'--','LineWidth',2)
plot(vel,freq_2,'--','LineWidth',2)
plot(vel,freq_3,'--','LineWidth',2)
plot(vel,freq_4,'--','LineWidth',2)
plot(vel,freq_5,'--','LineWidth',2)
plot(vel,freq_6,'--','LineWidth',2)
plot(vel,freq_7,'--','LineWidth',2)
plot(vel,freq_8,'--','LineWidth',2)
plot(vel,freq_9,'--','LineWidth',2)
plot(vel,freq_10,'--','LineWidth',2)
plot(vel,freq_11,'--','LineWidth',2)
plot(vel,freq_12,'--','LineWidth',2)
plot(vel,freq_13,'--','LineWidth',2)
plot(vel,freq_14,'--','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Frequency [Hz]')
xlim([0, 110])
legend('freq 1','freq 2','freq 3','freq 4','freq 5','freq 6','freq 7','freq 8','freq 9','freq 10','freq 11','freq 12','freq 13','freq 14', 'Location','east')
title('Frequency')
saveas(gcf, 'Immagini_risultati/full_bridge/frequency.png')

%%

figure()
hold on; grid on
plot(vel,smo_1,'--','LineWidth',2)
plot(vel,smo_2,'--','LineWidth',2)
plot(vel,smo_3,'--','LineWidth',2)
plot(vel,smo_4,'--','LineWidth',2)
plot(vel,smo_5,'--','LineWidth',2)
plot(vel,smo_6,'--','LineWidth',2)
plot(vel,smo_7,'--','LineWidth',2)
plot(vel,smo_8,'--','LineWidth',2)
plot(vel,smo_9,'--','LineWidth',2)
plot(vel,smo_10,'--','LineWidth',2)
plot(vel,smo_11,'--','LineWidth',2)
plot(vel,smo_12,'--','LineWidth',2)
plot(vel,smo_13,'--','LineWidth',2)
plot(vel,smo_14,'--','LineWidth',2)
plot(vel,zeros(length(vel),1),'k','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Damping \zeta [-]')
legend('smo 1','smo 2','smo 3','smo 4','smo 5','smo 6','smo 7','smo 8','smo 9','smo 10','smo 11','smo 12','smo 13','smo 14','Location', 'northwest')
title('Damping')
saveas(gcf, 'Immagini_risultati/full_bridge/smo.png')
%%
figure()
hold on; grid on
plot(vel,smo_1,'--','LineWidth',2)
plot(vel,smo_6,'--','LineWidth',2)
xline(66)

plot(vel,zeros(length(vel),1),'k','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Damping \zeta [-]')
legend('smo 1','smo 6','Location', 'northwest')
title('Damping')
saveas(gcf, 'Immagini_risultati/full_bridge/smo_focus.png')
%%
figure()
hold on
for jj = 1:28
plot3(real(autov(jj,:)),imag(autov(jj,:)),vel,'o','LineWidth',1)
end
view(2)
grid on
xline(0)
xlim([-0.275 0.2])
xlabel('Real')
ylabel('Imag')
title('Eigenvalues')
legend('modo 1','modo 2','modo 3','modo 4','modo 5','modo 6','modo 7','modo 8','modo 9','modo 10','modo 11','modo 12','modo 13','modo 14')
saveas(gcf, 'Immagini_risultati/full_bridge/modi.png')
%%
figure()
hold on
for jj = [1,2,7,8,11,12]
plot3(real(autov(jj,:)),imag(autov(jj,:)),vel,'o','LineWidth',1)
end
view(2)
grid on
xline(0)
xlim([-0.275 0.2])
xlabel('Real')
ylabel('Imag')
title('Eigenvalues')
legend('modo 1','modo 2','modo 7','modo 8','modo 11','modo 12')
saveas(gcf, 'Immagini_risultati/full_bridge/modi_focus.png')

toc