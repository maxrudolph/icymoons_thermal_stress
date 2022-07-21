function g = gravity(r,parameters)
% compute g at depth z
Ro = parameters.Ro;
g0 = parameters.g;
G = 6.6743015e-11; % N m^2/kg^2
rho_i = 917;

% compute m
M = g0*Ro^2/G;
% compute mass exterior to r

Moutside = rho_i * 4/3*pi*(Ro^3-r^3);
Mice = rho_i * 4/3*pi*(Ro^3-parameters.Rc^3);
Mrock = M-Mice;
rho_rock = Mrock/(4/3*pi*parameters.Rc^3);
if r < parameters.Rc
    Moutside = rho_i * 4/3*pi*(Ro^3 - parameters.Rc^3);
    Moutside = Moutside + rho_rock*4/3*pi*(parameters.Rc^3 - r^3);
end

g = G*(M-Moutside)/(r^2);
