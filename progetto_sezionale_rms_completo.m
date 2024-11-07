%%Wind Eng Lab - Bridge Aeroelasticity
clc
clear 
close all
%%
tic

%Constants
rho = 1.225; % air density [kg/m3]   
%B, deck chord [m]
%L, length [m]
L = 1;
load('struct_1.mat')

%Structural characteristics 
my = 23160;% [kg/m]
mz = my;
m = my;
J  = 2.77e6;% [kg m^2 / m] 

%frequenza
fy = 0.05;      
fz = 0.0884;
ft = 0.259;

zeta  = 0.004;     % damping ratio:  zeta = r / (2 m omega)

% N.B. [w = 2*pi*f]
wy = 2*pi*fy;
wz = 2*pi*fz;
wt = 2*pi*ft;

ry =2*m*wy*zeta;
rz =2*m*wz*zeta;
rt =2*J*wt*zeta;

%da definizione: omega = sqrt(k/m)
ky = (wy)^2*my;
kz = (wz)^2*mz;
kt = (wt)^2*J;

%Structural Matrices
M_stru = [m 0 0;
          0 m 0;
          0 0 J];
R_stru = [ry 0 0;
          0 rz 0;
          0 0 rt];
K_stru = [ky 0 0;
          0 kz 0;
          0 0 kt];
    
%Aerodynamic coefficients, li ricavo da file esterno.mat
load('STA_1.mat'); %alpha is in [deg]!

Cd=DLM(:,1);
Cl=DLM(:,2);
Cm=DLM(:,3);
Cd_alpha=KDLM(:,1);
Cl_alpha=KDLM(:,2);
Cm_alpha=KDLM(:,3);

figure(1)
plot(alpha,Cd,'-or','LineWidth',1)
hold on
plot(alpha,Cl,'-og','LineWidth',1)
plot(alpha,Cm,'-ob','LineWidth',1)
legend('CD','CL','CM', 'Location','southeast');
grid on;title('D,L,M')   
xlabel('Angle of attack \alpha [deg]','FontSize',10)
ylabel('[-]','FontSize',10)
saveas(gcf, 'Immagini_risultati/section/Aero_coeff.png')



% Plot derivatives of the aerodynamic coefficients w.r.t. the A.o.A
figure(2)
plot(alpha,Cd_alpha,'-or','LineWidth',1)
hold on
plot(alpha,Cl_alpha,'-og','LineWidth',1)
plot(alpha,Cm_alpha,'-ob','LineWidth',1)
legend('CD_a','CL_a','CM_a', 'Location','southeast');
grid on;title('D,L,M derivatives')   
xlabel('Angle of attack \alpha [deg]','FontSize',10)
ylabel('[-]','FontSize',10)
saveas(gcf, 'Immagini_risultati/section/Aero_coeff_derivatives.png')

%% STATIC SOLUTION

% NON CAMBIARE MAI NELLA VITA! (oppure rileggi il codice da capo e riscrivilo)
% se proprio devi: tutti numeri interi e con
% nvel = vel_max+1. ovviamente vel_min deve rimanere uguale a zero!

vel_min = 0;  vel_max = 70;  nvel = 71;
vel = linspace(vel_min,vel_max,nvel);  %wind speed vector 

% simplified static solution CM
CM0 = interp1(alpha,Cm,0); %calcola il Cm in alpha = zero facendo una interpolazione lineare
Cm_alpha0 = interp1(alpha,Cm_alpha,0); %calcola la derivata del Cm in alpha = zero  facendo una interpolazione lineare
% Linear (ideal) behaviour of Cm
CM_lin = CM0 + Cm_alpha0 * [-10:0.1:10]*pi/180; %costruisco vettore del CM_linearizzato che va da -10 a 10 gradi

% Plot ideal VS real behaviour of C_M
figure(3)
plot([-10:0.1:10], CM_lin,'-r','LineWidth',2,'DisplayName','C_{M0} + K_{M0} \theta_{0}')
hold on
plot(alpha,Cm,'-ob','LineWidth',1)
legend('CM0+Cm_alpha0*theta','CM(theta)')
xlabel('Angle of attack \theta_{0} [deg]','FontSize',10)
title('CM linearizzato')
grid on

% --------- Linearized solution: ---------
% Determine the solution assuming Cm linear (not function of theta)

%implement the linearized solution
th_simple = (0.5*rho*vel.^2*B^2*CM0)./(kt-0.5*rho*vel.^2*B^2*Cm_alpha0); %[rad]
th_simple_deg = rad2deg(th_simple);

% Plot the linearized solution
figure(4)
plot(vel,th_simple_deg)
xlabel('Velocity U [m/s]')
ylabel('Angle of attack \theta_{0} [deg]')
title('Linearized behaviour of the twist angle')
grid on

% nonlinear solution
sol_old=0;%guess iniziale, non è a caso 0, infatti il vettore velocità parte da zero
theta_NON_linear_deg=zeros(1,length(vel));%inizializzo
hh = waitbar(0,'Computing static displacemnts (section model)...'); %serve per far comparire la finestra per visualizzare lo stato d'avanzamento del ciclo

steps = length(vel);
for ii=1:length(vel)
    % ii-th velocity
    U = vel(ii);
    % NB : THETA and ALPHA are IN DEG, torsional stiffness is in Nm/rad
	% fsolve find the nonlinear solution of the problem defined in the function 'statica', that is in another m-file
    [equilibrio] = fsolve(@(theta) statica(K_stru,theta,rho,U,B,alpha,Cm), sol_old);%
    theta_NON_linear_deg(ii)= equilibrio; % theta static, in deg
    sol_old=equilibrio;       
    waitbar(ii/steps)
end
close(hh)

% Now that you have the theta modified based on the velocity you can
% compute the new quasi-linear coefficients

% The other displacements and coefficients come as consequence since you got all the theta 
Cd_theta = interp1(alpha,Cd,theta_NON_linear_deg); % Cd(theta)
Cl_theta = interp1(alpha,Cl,theta_NON_linear_deg); % Cl(theta)
y_statico = 0.5*rho*vel.^2*B.*(Cd_theta)./ky;
z_statico = 0.5*rho*vel.^2*B.*(Cl_theta)./kz;

% Plot of the coefficients
figure()
subplot(311)

plot(vel,y_statico,'r','LineWidth',1)
hold on
grid on
xlabel('Mean velocity U [m/s]')
ylabel('y [m]')
title('Static displacements under mean velocity');

subplot(312)
plot(vel,z_statico,'g','LineWidth',1)
hold on
grid on
xlabel('Mean velocity U [m/s]')
ylabel('z [m]')

subplot(313)
plot(vel,theta_NON_linear_deg,'b','LineWidth',1)
hold on
grid on
xlabel('Mean velocity U [m/s]')
xlabel('Mean velocity U [m/s]')
ylabel('\theta [deg]')
hold on
plot(vel,th_simple_deg,'--b')
legend('non linear', 'linear', 'Location','southwest')

saveas(gcf, 'Immagini_risultati/section/Stati_displacement.png')


%% FLUTTER
% B1 for QuasiSteadyTheory
B1y=0*B; 
B1z=0*B;
B1t=0.3*B;

% Damping vectors
damp_y=zeros(1,length(vel)); 
damp_z=zeros(1,length(vel));
damp_t=zeros(1,length(vel));

% Frequency vectors
freq_y=zeros(1,length(vel)); 
freq_z=zeros(1,length(vel));
freq_t=zeros(1,length(vel));

% Eigenvalues
autoval=zeros(6,length(vel));

flutter_speed = [];

hh = waitbar(0,'Calcolo velocità critica di flutter (SEZIONALE)...');

steps = length(vel);
%qua inizia un ciclo per ogni velocità
for ii=1:length(vel)
        % ii-th velocity
        U    = vel(ii); 
        % Dynamic pressure
        qBL = 0.5*rho*U.^2*B*L ;   
        
        % Non-linear twist (twist associated to the ii-th speed)
        theta_eq = theta_NON_linear_deg(ii);
        
        %coeff. aerodinamici di equilibrio
        Cd_eq = interp1(alpha, Cd, theta_eq);
        Cl_eq = interp1(alpha, Cl, theta_eq);
        Cm_eq = interp1(alpha, Cm, theta_eq);
        Cd_alpha_eq = interp1(alpha, Cd_alpha, theta_eq);
        Cl_alpha_eq = interp1(alpha, Cl_alpha, theta_eq);
        Cm_alpha_eq = interp1(alpha, Cm_alpha, theta_eq);

        % Define the Aero Matrices
        M_aero = zeros(3,3);
        R_aero = 0.5*rho*U*B  * [2*Cd_eq, (Cd_alpha_eq-Cl_eq), (Cd_alpha_eq-Cl_eq)*B1y;
                                 2*Cl_eq, (Cl_alpha_eq+Cd_eq), (Cl_alpha_eq+Cd_eq)*B1z;
                                 B*2*Cm_eq,   B*Cm_alpha_eq,    B*Cm_alpha_eq*B1t];
                             
        K_aero = 0.5*rho*U^2*B* [0, 0, -Cd_alpha_eq;
                                 0, 0, -Cl_alpha_eq;
                                 0, 0, -B*Cm_alpha_eq];
                             
        % Build up the total matrices of the resolutive equation
        M = M_stru + M_aero ; % Mass
        R = R_stru + R_aero ; % Damping
        K = K_stru + K_aero ; % Stiffness
       
        [autovet,autoval]=polyeig(K,R,M);
                
                if ii==1
                    f_old = [fy fz ft]; 
                    smo_old = [zeta zeta zeta];
                    lambda_old(1:3) = -smo_old.* (2*pi*f_old) + 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
                    lambda_old(4:6) = -smo_old.* (2*pi*f_old) - 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
                end
                
                indici = 1:6;
                for jj = 1:6
                    [val,dove]=min(abs(autoval(indici)-lambda_old(jj)));
                    ind(jj) = indici(dove);
                    indici = setdiff(indici,ind(jj));
                end
                autoval = autoval(ind);
                autovet = autovet(:,ind);
     
        damp_y(ii)= -real(autoval(1))/norm(autoval(1));
        damp_z(ii)= -real(autoval(3))/norm(autoval(3));
        damp_t(ii)= -real(autoval(5))/norm(autoval(5));
        freq_y(ii)= abs(imag(autoval(1)))/(2*pi);
        freq_z(ii)= abs(imag(autoval(3)))/(2*pi);
        freq_t(ii)= abs(imag(autoval(5)))/(2*pi);


        autov(:,ii) = autoval;
        lambda_old  = autoval;
        
        if (damp_t(ii)<=0)
            flutter_speed = [flutter_speed, vel(ii)];
        end

        waitbar(ii/steps)
end
close(hh)
        

% Plot the results of FLUTTER 

figure() % Frequency
hold on; grid on
plot(vel,freq_y,'--','LineWidth',2)
plot(vel,freq_z,'--','LineWidth',2)
plot(vel,freq_t,'--','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Frequency [Hz]')
legend('y','z','t')
title('Frequency')

saveas(gcf, 'Immagini_risultati/section/Frequency.png')


figure() % Damping
hold on; grid on
plot(vel,damp_y,'--','LineWidth',2)
plot(vel,damp_z,'--','LineWidth',2)
plot(vel,damp_t,'--','LineWidth',2)
plot(vel,zeros(length(vel),1),'k','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Damping \zeta [-]')
legend('y','z','t', 'Location','west')
title('Damping')
saveas(gcf, 'Immagini_risultati/section/Damping.png')
%%

figure() % Eigenvalues
hold on
for jj = 1:6
plot3(real(autov(jj,:)),imag(autov(jj,:)),vel,'-o','LineWidth',2)
end
view(2)
grid on
xlabel('Real')
ylabel('Imag')
xline(0)
title('Eigenvalues')

saveas(gcf, 'Immagini_risultati/section/Eigenvalues.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- Buffeting Response of Bridge section ----------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%% per ogni velocità di buffeting vanno prese le relative u e ww

qpsd_totale = [];
rms_totale = [];

for kk = 1:1:vel_max
        
    U = kk;
    
    theta_eq = theta_NON_linear_deg(U+1);
    
    Cd_eq = interp1(alpha, Cd, theta_eq);
    Cl_eq = interp1(alpha, Cl, theta_eq);
    Cm_eq = interp1(alpha, Cm, theta_eq);
    Cd_alpha_eq = interp1(alpha, Cd_alpha, theta_eq);
    Cl_alpha_eq = interp1(alpha, Cl_alpha, theta_eq);
    Cm_alpha_eq = interp1(alpha, Cm_alpha, theta_eq);
            
    Matrix_buff = 0.5*rho*U^2*B*[2*Cd_eq (Cd_alpha_eq - Cl_eq)
                                 2*Cl_eq (Cl_alpha_eq + Cd_eq)
                                 2*Cm_eq*B  B*Cm_alpha_eq];
    % Compute the wind speeds u and w with the Von Karman Spectrum
    [f1,u_totale] = u_piccolo_totale(vel_max);% [m/s]
    [f2,w_totale] = w_piccolo_totale(vel_max);% [m/s]
    
    Vect_buff = [u_totale(kk, :)/U; w_totale(kk, :)/U];
                         
    F_buff = Matrix_buff * Vect_buff;
    
    R_aero = 0.5*rho*U*B  * [2*Cd_eq, (Cd_alpha_eq-Cl_eq), (Cd_alpha_eq-Cl_eq)*B1y;
                             2*Cl_eq, (Cl_alpha_eq+Cd_eq), (Cl_alpha_eq+Cd_eq)*B1z;
                             B*2*Cm_eq,   B*Cm_alpha_eq,    B*Cm_alpha_eq*B1t];
                         
    K_aero = 0.5*rho*U^2*B* [0, 0, -Cd_alpha_eq;
                             0, 0, -Cl_alpha_eq;
                             0, 0, -B*Cm_alpha_eq];
    M_10 = M_stru;
    R_85 = R_stru + R_aero;
    K_85 = K_stru + K_aero;
    
    q = [];
    
    n_plot = 300;
    om_plot = 2*pi*f1(1:1:n_plot);
    
    for a=1:n_plot
        w = 2*pi*f1(a);
        q = [q, abs((-w^2*M_10+1i*w*R_85+K_85)\F_buff(:, a))];
    end
    
    f_plot = om_plot/2/pi; 
    
    if (kk == 10) 
        figure(101)
        plot(f_plot, q(1,:),'Linewidth', 2)
        hold on
        grid on
        grid minor
        plot(f_plot, q(2,:),'Linewidth', 2)
        plot(f_plot, q(3,:)*B/2,'Linewidth', 2)
        xline(fy)
        xline(fz)
        xline(ft)
        legend('y','z','t', 'fy','fz','ft')
        title('fft 10 m/s')
        xlabel('f [Hz]')
        ylabel('[m]')

        saveas(gcf, strcat('Immagini_risultati/section/fft_', num2str(kk), 'ms.png'))
    

        figure(102)
        plot(f_plot, qpsd(1,:),'Linewidth', 2)
        hold on
        plot(f_plot, qpsd(2,:),'Linewidth', 2)
        hold on
        plot(f_plot, qpsd(3,:)*B/2,'Linewidth', 2)
        grid on
        grid minor
        xline(fy)
        xline(fz)
        xline(ft)
        legend('y','z','t', 'fy','fz', 'ft')
        tit = title('psd 10 m/s');
        xlabel('f [Hz]')
        ylabel('[m^2/Hz]')

        saveas(gcf, strcat('Immagini_risultati/section/psd_', num2str(kk), 'ms.png'))
    end
    
    if (kk == 30)
        figure(301)
        plot(f_plot, q(1,:),'Linewidth', 2)
        hold on
        grid on
        grid minor
        plot(f_plot, q(2,:),'Linewidth', 2)
        plot(f_plot, q(3,:)*B/2,'Linewidth', 2)
        xline(fy)
        xline(fz)
        xline(ft)
        legend('y','z','t', 'fy','fz', 'ft')
        title('fft 30 m/s')
        xlabel('f [Hz]')
        ylabel('[m]')

        saveas(gcf, strcat('Immagini_risultati/section/fft_', num2str(kk), 'ms.png'))
    

        figure(302)
        plot(f_plot, qpsd(1,:),'Linewidth', 2)
        hold on
        plot(f_plot, qpsd(2,:),'Linewidth', 2)
        hold on
        plot(f_plot, qpsd(3,:)*B/2,'Linewidth', 2)
        grid on
        grid minor
        xline(fy)
        xline(fz)
        xline(ft)
        legend('y','z','t', 'fy','fz', 'ft')
        title('psd 30 m/s')
        xlabel('f [Hz]')
        ylabel('[m^2/Hz]')

        saveas(gcf, strcat('Immagini_risultati/section/psd_', num2str(kk), 'ms.png'))
    end
    
    if (kk == 61)
        figure(601)
        plot(f_plot, q(1,:),'Linewidth', 2)
        hold on
        grid on
        grid minor
        plot(f_plot, q(2,:),'Linewidth', 2)
        plot(f_plot, q(3,:)*B/2,'Linewidth', 2)
        xline(fy)
        xline(fz)
        xline(ft)
        legend('y','z','t', 'fy','fz', 'ft')
        title('fft 61 m/s')
        xlabel('f [Hz]')
        ylabel('[m]')

        saveas(gcf, strcat('Immagini_risultati/section/fft_', num2str(kk), 'ms.png'))

        
        figure(602)
        plot(f_plot, qpsd(1,:),'Linewidth', 2)
        hold on
        plot(f_plot, qpsd(2,:),'Linewidth', 2)
        hold on
        plot(f_plot, qpsd(3,:)*B/2,'Linewidth', 2)
        grid on
        grid minor
        xline(fy)
        xline(fz)
        xline(ft)
        legend('y','z','t', 'fy','fz', 'ft')
        title('psd 61 m/s')
        xlabel('f [Hz]')
        ylabel('[m^2/Hz]')

        saveas(gcf, strcat('Immagini_risultati/section/psd_', num2str(kk), 'ms.png'))
    end
    
    df = f1(1);
    qpsd = [q(1,:).^2 / (2*df); q(2,:).^2 / (2*df); q(3,:).^2 / (2*df); ];
    rms = sqrt(trapz(f_plot,qpsd')); % ogni termine è riferito a ogni riga di qpsd
    
    qpsd_totale = [qpsd_totale; qpsd];
    rms_totale = [rms_totale; rms];
    
    clear q
end




%% plot RMS
n = vel_max;

figure(111)
plot(1:1:n, rms_totale(:,1),'Linewidth',2)
hold on
plot(1:1:n, rms_totale(:,2),'Linewidth', 2)
hold on
plot(1:1:n, rms_totale(:,3),'Linewidth', 2)
xline(61)
grid on
grid minor
legend('y','z','t','61 [m/s]')
title('rms')
xlabel('U [m/s]')
ylabel('[m]')
toc

saveas(gcf, 'Immagini_risultati/section/RMS.png')