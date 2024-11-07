function [f1,u_totale] = u_piccolo_totale(vel_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate the fluctuation with Von Karman spectrum at z=68m
% INPUTS: 
% U --> wind-speed
% 
% OUTPUTS:
% w_totale --> Von Karman spectrum for w
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
data = readmatrix('bridge_1\wind&turbulence profiles.xlsx');
z0 = 0.01; % Roughness length (category II) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
z = 68; % [m] eight bridge

% Turbulence intensity vertical profile
idx = find(data(:,1)==70);
Iu = data(idx,3);%[-]

% Integral length scale vertical profile

Lu = data(idx,5);%[m]

for k = 1:1:vel_max

U_mean = k;%[m/s]
U = U_mean*1;
% Lu= 32.9; %[m]
% Iu= 0.062; %[-]

df=1/600; %frequency resolution
fsamp=10;

f=[df:df:fsamp/2]; %frequency array [Hz]
fn=f*Lu/U; %non dimensional frequency
SuuN=4*fn./(1+70.8*fn.^2).^(5/6); %Von Karman spectrum

%Plot Von Karman spectrum
% figure
% loglog(fn,SuuN,'b'),grid on
% xlabel('f*L/U')
% ylabel('f*Suu/\sigma_u^2')

%computation of dimensional spectrum
varu=(Iu*U)^2; %variance of u - da definizione di turbulence intensity [slide 57 di 02b]
Suu=SuuN./f*varu; %dimensional spectrum

%plot dimensional spectrum
% figure
% loglog(f,Suu,'b'),grid on
% xlabel('f [Hz]')
% ylabel('Suu [(m/s)^2/Hz]')

%computation of the harmonics
A=sqrt(Suu*df*2); %[m/s] harmonics of the spectrum

%plot harmonics
% figure
% stem(f,A,'b'),grid on
% xlabel('f [Hz]')
% ylabel('A [m/s]')

u_totale(k, :) = A;
f1 = f;

end


%%
%save sim1_totale.mat f1 u_totale