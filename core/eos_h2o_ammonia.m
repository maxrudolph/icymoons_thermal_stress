
p=1.;

n=30;
for i=1:n+1
    x=0.3*(i-1)/n;
    [tm]=melting(x);
    [rho]=density(x,p,tm);
    tmplot(i)=tm;
    rhoplot(i)=rho;
    xplot(i)=x;
end

plot(xplot,rhoplot)
hold on


function [tm] = melting(X)
% melting point of ammonian-water mixtures
% X is the mass fraction of ammonia
% equation (3a) from Croft et al. (1988)
% valid for X<0.329
%
tm=273.15-53.07*X-1651.4*X^2+11842.0*X^3-46269.0*X^4+56695.0*X^5;
end

function [rho] = density(X,P,T)
% density of ammonia-water mixtures from EOS in Croft et al. (1988)
% X mass fraction ammonia
% P pressure [bar]
% T temperature [Kelvin]
% ko bulk modulus [kbar I think]
% ko d(ko)/dP unitless

[tm]=melting(X);

% calculate bulk moduli and derivatives
konh3=48.503*exp(-1.0134e-3*T^1.3067); % equation 4a
kopnh3=(4.0858e16*T^-6.6083)+4.1831; % equation 4b
koh2o=-73.184+0.591*T-9.139e-4*T^2; % equation 4c
koph2o=39.999-0.2058*T+3.111e-4*T^2; % equation 4d
C=-638.89+1.9519*T; % equation 5b
D=316.84-0.99616*T; % equation 5c
ko=koh2o*(1.0-X)+konh3*X+C*(X^2-X)+D*(X^3-X); % equation 5a
kop=(kopnh3^X)*(koph2o^(1.0-X)); % equation 5d

% calculate zero pressure density rho0
a1=-1.0e-6*(-92.88+1371.05*X+185.91*X^2); % equation 2c
a2=-1.0e-6*(14.51-47.5*X+42.35*X^2); % equation 2d
a3=-1.0e-6*(-0.0764+0.3118*X-0.2465*X^2); % equation 2e
rhor=0.9991-0.4336*X+0.3303*X^2+0.2833*X^3-1.9716*X^4+2.1396*X^5-0.7294*X^6; % equation 2b
rho0=rhor*exp(a1*(T-288.15))+0.5*a2*((T-tm)^2-(288.15-tm)^2)+(a3/3.0)*((T-tm)^3-(288.15-tm)^3); % euqation 2a

rho=rho0*((kop*P/ko)+1.0)^(1.0/kop); % eqiation (1)

end