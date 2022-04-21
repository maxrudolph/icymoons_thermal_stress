function dUdz = ocean_dUdz(Ri,Rc,rho_i,g,Pex,X,X0,V0)
% Ri: inner radius in m
% Rc: core radius in m
% P: total pressure in Pa
% X: mass fraction ammonia
% calculates energy released per unit thickening
%
Cp_h20 = 4182; %J/kg/K
Cp_Nh3 = 4744; %J/kg/K

Cp_eff = Cp_h20*(1-X) + Cp_Nh3*X;% J/kg/K
rho = ammonia_density(X,P,T);

V = 4/3*pi*(Ri^3-Rc^3); % current volume.

% calculate change in ocean temperature with respect to unit change in
% thickness. dT/dz = dT/dP * dP/dz where dP/dz = rho_i*g
% dX/dz = concentration of ammonia in solidifying ocean.

% dV/dz = -4*pi*Ri^2
dVdz = -4*pi*Ri^2;
% dX/dV = -X0*V0/V^2
dXdV = -X0*V0/V^2;

dUdz = Cp*rho
