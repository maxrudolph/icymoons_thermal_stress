function [P,T,z,rho] = ammonia_adiabatic_profile(X,Pocean,g)
% (dT/dP)_S = T*alpha/(rho*Cv)
% ocean pressure should be in Pascal
Tm = ammonia_melting(X);
% use the ammonia_density function to calculate an approximate thermal
% expansivity
n = 1000; % number of sample points
% make a pressure linspace
P = linspace(Pocean,.006,1000);
T = zeros(size(P));
T(1) = Tm;
rho = zeros(size(P));
rho(1) = ammonia_density(X,P(1),T(1));
z = zeros(size(P));

for i=2:n
    dTdP = adiabat(X,P(i-1),T(i-1));
    T(i) = T(i-1) + dTdP*(P(i)-P(i-1));
    rho(i) = ammonia_density(X,P(i),T(i));
    z(i) = z(i-1) + -1/(rho(i)*g)*(P(i)-P(i-1));
end


end

function dTdP = adiabat(X,P,T)
% T must be in Kelvin
% P must be in Pascal
% X is mass fraction ammonia.
V_c = 1/ammonia_density(X,P,T);
delta = 1e-3;
V_Tp = 1/ammonia_density(X,P,T+delta);
V_Tm = 1/ammonia_density(X,P,T-delta);
alpha = 1/V_c * (V_Tp-V_Tm)/(2*delta);

dTdP = T*alpha/(1/V_c * ammonia_cp(X));% note 1/Vc = rho
end