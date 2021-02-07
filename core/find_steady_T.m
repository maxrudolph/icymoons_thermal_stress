function [qb,T] = find_steady_T(Ri,Ro,Tb,Ts,r)
% Solve for the steady-state basal heat flux q0.
% Assumes k = k1/T0 (Petrenko and Whitworth 1999)
% k = 651./T
% Analytic solution is worked out in OneNote Enceladus document.
% ln(T) = c1+c0/r.
k1 = 651.0;

% h = 30000;
% Ro = 1560*1000;
% Ri = Ro-h;
% Ts = 100;
% Tb = 273;

c0 = Ri*Ro/(Ri-Ro)*(log(Ts)-log(Tb));
c1 = log(Ts)-c0/Ro;

lnT = c1 + c0./r;
T = exp(lnT);
figure();
plot(r,T);

qb = k1*c0/Ri^2;
qs = k1*c0/Ro^2;
