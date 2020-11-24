function T = solve_stefan_analytic(depth,k,rho,Cp,L,Tm,T0)
% function solve_stefan
% input kappa
kappa = k/rho/Cp;
% lambda1 is calculated following Turcotte and Schubert equation 4.141.
f = @(lambda1) L*sqrt(pi)/Cp/(Tm-T0) - exp(-lambda1^2)/(lambda1*erf(lambda1));
% find lambda1 such that f(lambda1) == 0
lambda1 = fzero(f,1.0); % initial guess is 1.0;
fprintf('lambda1 = %.02f\n',lambda1);

tm = (max(depth)/2/lambda1/sqrt(kappa))^2;
eta = depth/2/sqrt(kappa*tm);
theta = erf(eta)/erf(lambda1);
T = theta*(Tm-T0) + T0;
