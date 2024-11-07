function f = statica(K_stru,theta,rho,U,B,alpha,Cm)
% theta is in [deg], but should be converted in [rad] when you compute K_theta*theta
% to find Cm(theta) use the 'interp1' function
% fsolve finds theta that makes f=0, so write the static equilibrium accordingly

theta_rad = deg2rad(theta);
Cm_theta = interp1(alpha,Cm,theta);
f = K_stru(end,end)*theta_rad - 0.5*rho*U^2*B^2*Cm_theta;

return