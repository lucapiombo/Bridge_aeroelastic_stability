function g = statica_modale_full (K_diag,q0,rho,U,B,L_sez,Cd,Cl,Cm,alpha,PHI)

x = PHI*q0;

theta = [];

for k = 1:243
theta = [theta; x(3*k)];
end
theta = rad2deg(theta);

F_vett_sezionale = [];
for k = 1:243
Cd_eq(k) = (interp1(alpha,Cd,theta(k)));
Cl_eq(k) = (interp1(alpha,Cl,theta(k)));
Cm_eq(k) = (interp1(alpha,Cm,theta(k)));
F_vett_sezionale = [F_vett_sezionale; 0.5*rho*U^2*L_sez*B*[Cd_eq(k) ; Cl_eq(k) ; B*Cm_eq(k)]];
end

F_vett_totale = F_vett_sezionale;

g = K_diag*q0-PHI'*F_vett_totale;
return