function [ifail] = tensile_failure_criterion(z,sigma_t,rho,g,tensile_strength)
% Determine whether failure occurs for a given choice of parameters
% z is a (vector) of depth values
% sigma_t is a (vector) of stress values in the ice shell
% rho, g are density and gravitational acceleration
% tensile_strength is the tensile strength.
assert(all(size(z) == size(sigma_t)),'Error: sigma_t and z must have the same dimensions')
lithostatic_pressure = rho*g*z;
ifail = (sigma_t - lithostatic_pressure) > tensile_strength;