% test implementation of the Stefan solution
clear
close all

nz = 100;
H=0;
Ts=100;
Tb=270;
rho=1000;
Cp = 2100; %heat capacity of ice J/kg/K
Lf = 334*1000; % latent heat of fusion (J/kg)
H0=2400; % initial thickness in m
seconds_in_year = 3.15e7;
zmax = 1e6;
zmin = zmax-H0;
dt = 1000*seconds_in_year;
tmax = 10e6*seconds_in_year;
kappa=1e-6;
k=kappa*rho*Cp;

% find lambda1 (dimensionless parameter in stefan solution) by bisection:
% first, solve for lambda1
lam1_max = 1.6;
lam1_min = 0.1;
residual = @(lam1) exp(-lam1^2)/lam1/erf(lam1) - Lf*sqrt(pi)/(Cp*(Tb-Ts)); % turcotte and schubert 4.141 expressed as a residual
while lam1_max-lam1_min > 1e-6
    lam1_guess = 1/2*(lam1_max+lam1_min)
    ub = residual(lam1_max);
    lb = residual(lam1_min);
    if sign(residual(lam1_guess)) == sign(ub)
        lam1_max = lam1_guess;
    else
        lam1_min = lam1_guess;
    end
end
lam1 = lam1_guess;

cooling_age = @(zm) (zm/2/lam1)^2/kappa;
initial_cooling_age = cooling_age(H0);
zm = @(t) 2*lam1*sqrt(kappa*t);
dTdt = @(z,t) -(Tb-Ts)/erf(lam1)*exp(-z.^2/4/kappa/t).*z/(2*sqrt(pi*kappa*t^3));

z = linspace(zmin,zmax,nz);
z_last = z;
% stefan initial condition
T_last = zeros(nz,1);
eta = (zmax-z)/2/sqrt(kappa*initial_cooling_age);
theta = erf(eta)/erf(lam1);
T_last(:) = theta*(Tb-Ts)+Ts;

f1=figure(1);
plot(z,T_last); hold on

time=0;
itime=1;
while( time < tmax )
    % advance mesh
    % 1. calculate dTdz
    dTdz = zeros(nz,1);
    for i=2:nz-1
        dz = z_last(i+1)-z_last(i-1);
        dTdz(i) = (T_last(i+1)-T_last(i-1))/dz;
    end
    Tg = Tb - (T_last(2)-Tb);
    dTdz(1) = (T_last(2)-Tg)/2/(z_last(2)-z_last(1));
    
    q = -k*dTdz;
    % calculate thickening rate
       delta_z = q(1)/rho/Lf*dt
%     delta_z = 4;
    if( delta_z < 0 )
        itime
        time
        keyboard
    end
    
    
    %     % interpolate onto new grid
    %     z = linspace(z_last(1)-delta_z,zmax,nz);
    ztmp = [z_last(1)-delta_z z_last];
    %      ztmp(1) = ztmp(1)-delta_z;
    z=linspace(ztmp(1),zmax,nz);
    T_last = interp1(ztmp,[T_last(1) T_last'],z);
    z_last = z;
    
    %     Ttmp = T_last;
    % %     Ttmp(1) = Tb;
    %     T_last = interp1([ztmp],[Ttmp],z);
    %     z = z_last;
    %     z(1) = z(1)-delta_z;
    %     T_last = T_last;
    %
    L = zeros(nz,nz);
    R = zeros(nz,1);
    for i=1:nz
        if( i==1 )
            dzp = z(2)-z(1);
            dzm = z(2)-z(1);
            dz  = z(2)-z(1);
        elseif i==nz
            dzp = z(i)-z(i-1);
            dzm = z(i)-z(i-1);
            dz  = z(i)-z(i-1);
        else
            dzp = z(i+1)-z(i);
            dzm = z(i)-z(i-1);
            dz = 1/2*(dzp+dzm);
        end
        
        coef_minus  =            -k/dzm/dz;
        coef_center = rho*Cp/dt +k/dzm/dz + k/dzp/dz;
        coef_plus   =            -k/dzp/dz;
        rterm       = rho*Cp/dt*T_last(i);
        if i==1
%                         L(i,i) = coef_center;
%                         L(i,i+1) = coef_plus - coef_minus;
%                         R(i) = rterm + H -2*coef_minus*Tb;
            L(i,i) = coef_center;
            R(i)   = coef_center*Tb;
        elseif i==nz
%                         L(i,i-1) = coef_minus - coef_plus;
%                         L(i,i) = coef_center;
%                         R(i) = rterm + H -2*coef_plus*Ts;
            L(i,i) = coef_center;
            R(i)   = coef_center*Ts;
        else
            L(i,i-1) = coef_minus;
            L(i,i)   = coef_center;
            L(i,i+1) = coef_plus;
            R(i) = rterm + H;
        end
    end
    T = L\R;
    
    % advance timestep;
    time_store(itime) = time+dt;
    z_store(itime) = z(1);
    T_last = T;
    z_last=z;
    time = time+dt;
    itime = itime+1;
end

% plot analytic and numerical thickness
figure()
plot(time_store,(zmax-H0)-z_store); hold on;
plot(time_store, zm(initial_cooling_age+time_store)-H0,'r--');
err = norm((zmax-H0)-z_store - (zm(initial_cooling_age+time_store)-H0))/norm( zm(initial_cooling_age+time_store)-H0 )

figure, plot(z,T);

figure, plot(z,dTdz)

% 2nd method - note that even though T is not correct at the boundary, the
% approximation dt/dx is consistent.
dTdz = zeros(nz,1);
for i=2:nz-1
    dTdz(i) = (T(i+1)-T(i-1))/dz/2;
end

dTdz(nz) = (T(nz)-T(nz-1))/dz;
figure, plot(z,dTdz)
