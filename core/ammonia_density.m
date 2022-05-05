function rho = ammonia_density(X,P,T)
% density of ammonia-water mixtures from EOS in Croft et al. (1988)
% X mass fraction ammonia
% P pressure [kbar]
P = P/1e8;% convert pressure from kbar to Pa (1 kbar = 10^8 Pa)
% T temperature [Kelvin]
% ko bulk modulus [kbar]
% ko d(ko)/dP unitless


[tm]=ammonia_melting(X);

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

rho = rho*1000.0;% convert g/cc -> kg/m^3
end

