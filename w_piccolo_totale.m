
function [f2,w_totale] = w_piccolo_totale(vel_max)
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
Iw = data(idx,4);%[-]

% Integral length scale vertical profile

Lw = data(idx,6);%[m]

for k = 1:1:vel_max

U_mean = k;%[m/s]
U = U_mean*1;
% Lu= 32.9; %[m]
% Iu= 0.062; %[-]

df=1/600; %frequency resolution
fsamp=10;

f=[df:df:fsamp/2]; %frequency array [Hz]
fn=f*Lw/U; %non dimensional frequency
SwwN=4*fn./(1+70.8*fn.^2).^(5/6); %Von Karman spectrum

%Plot Von Karman spectrum
% figure
% loglog(fn,SwwN,'b'),grid on
% xlabel('f*L/U')
% ylabel('f*Sww/\sigma_u^2')

%computation of dimensional spectrum
varw=(Iw*U)^2; %variance of u - da definizione di turbulence intensity [slide 57 di 02b]
Sww=SwwN./f*varw; %dimensional spectrum

%plot dimensional spectrum
% figure
% loglog(f,Sww,'b'),grid on
% xlabel('f [Hz]')
% ylabel('Sww [(m/s)^2/Hz]')

%computation of the harmonics
A=sqrt(Sww*df*2); %[m/s] harmonics of the spectrum

%plot harmonics
% figure
% stem(f,A,'b'),grid on
% xlabel('f [Hz]')
% ylabel('A [m/s]')

w_totale(k, :) = A;
f2 = f;

end


%%
% save sim2_totale.mat f2 w_totale