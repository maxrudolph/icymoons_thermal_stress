function [h] = find_steady_h(Ro,Tb,Ts,qs)
% Solve for the steady-state ice shell thickness.
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

% To solve using qb (basal heat flux - doesn't work for thick ice shells)
% c = -k1*Ro./qb*log(Ts/Tb);
% h = (Ro - sqrt(Ro^2-4*c))/2;

h = (-qs - Ro*k1/Ro^2*log(Ts/Tb)).^-1 * k1/Ro*log(Ts/Tb);

% Uncomment for validation:
% Ri = Ro-h;
% r = linspace(Ri,Ro);
% 
% c0 = Ri*Ro/(Ri-Ro)*(log(Ts)-log(Tb));
% c1 = log(Ts)-c0/Ro;

% lnT = c1 + c0./r;
% T = exp(lnT);
% figure();
% plot(r,T);

% qb = k1*c0/Ri^2
% qs = k1*c0/Ro^2;
