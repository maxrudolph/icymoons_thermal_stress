% Script to solve coupled ice shell thermal and stress evolution
% This version of the code is meant to reproduce the benchmark in section
% 163 of Timoshenko and Goodier. The benchmark problem involves the steady
% conduction of heat across a spherical shell of finite thickness.
%
% Max Rudolph, March 19, 2020
clear;
close all;

% Numerical parameters
nr = 101; % number of grid points
relaxation_parameter=.01; % used in nonlinear loop.
maxiter=300;
% Define physical constants and parameters
% Physical constants
seconds_in_year = 3.1558e7;
R=8.314e-3;     % in kJ/mol/K
% Boundary conditions and internal heating
H=0; % internal heating rate.
Tb=270;
Ts=100;
Ro = 1000;             % outer radius of ice shell (m)
initial_thickness = 200;
Ri = Ro-initial_thickness;         % inner radius of ice shell (m)
Rc = 1.60e5;             % core radius (m)
% Elastic and Viscous properties
E = 5e9;        % shear modulus of ice (Pa)
nu = 0.3;      % Poisson ratio of ice (-)
beta_w = 4e-10; % Compressibility of water (1/Pa)
alpha_l = 1e-4; % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
rho_i=900;      % density of ice (kg/m^3)
rho_w=1000;     % density of water (kg/m^3)
Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=1e15;       % basal viscosity (Pa-s)
% mu = @(T) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
mu = @(T) 1e99;
% Thermal properties
Cp = 2100; %heat capacity of ice J/kg/K
Lf = 334*1000; % latent heat of fusion (J/kg)
kappa = 1e-6;% m/s/s
k=kappa*rho_i*Cp;

% Benchmark-specific things
% temperature solution
lam1 = 0.65;
cooling_age = @(zm) (zm/2/lam1)^2/kappa;
initial_cooling_age = cooling_age(initial_thickness);
zm = @(t) 2*lam1*sqrt(kappa*t);
dTdt = @(z,t) -(Tb-Ts)/erf(lam1)*exp(-z.^2/4/kappa/t).*z/(2*sqrt(pi*kappa*t^3));

% calculate maxwell time at 100, 270
fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(100)/E,mu(Tb)/E);
fprintf('Thermal diffusion timescale %.2e\n',(4e4)^2/kappa);
% set end time and grid resolution
% t_end = 90e6*seconds_in_year;
t_end = 1;
% dt = 1e3*seconds_in_year; % time step in seconds
% dt1 = 3600;
% times = logspace(log10(dt1),log10(t_end+dt1),1e4)-dt1;
times = [0 1];
plot_interval = t_end;
save_interval = 1e5*seconds_in_year;
save_depths = [0 1 2 3 5]*1000;
nsave = ceil(t_end/save_interval) + 1; nsave_depths = length(save_depths);
sigma_t_store = zeros(nsave_depths,nsave);

% set up the grid
grid_r = linspace(Ri,Ro,nr); % set up the grid

% initialize solution vectors (IC)
sigma_r_last = zeros(nr,1); % initial stresses
sigma_t_last = zeros(nr,1); % initial stresses
% stefan initial condition
eta = (Ro-grid_r)/2/sqrt(kappa*initial_cooling_age);
theta = erf(eta)/erf(lam1);
T_last = ones(nr,1)*Ts;
% initialize with a stefan solution corresponding to initial thickness
er_last = zeros(nr,1); % strains
et_last = zeros(nr,1);
ur_last = zeros(nr,1);      % displacement
z_last = 0; % total amount of thickening
Pex_last =0; %initial overpressure

% Set up plot
hf2=figure();

plot_times = [0.8 1.6 3.2 6.4 12 24 46 90]*1e6*seconds_in_year; iplot=1;
hf=figure(2);
subplot(1,4,1); % sigma_r and sigma_t
h=plot(sigma_r_last,Ro-grid_r); hold on;
plot(sigma_r_last,Ro-grid_r,'--','Color',h.Color);
% h=legend('\sigma_r','\sigma_t','Interpreter','tex'); h.AutoUpdate=false;
title('Stress (Pa)','Interpreter','tex');
ylabel('r (m)');
set(gca,'YDir','reverse');
subplot(1,4,2); % e_r and e_t
h=plot( sigma_r_last,Ro-grid_r); hold on;
plot( sigma_r_last,Ro-grid_r,'--','Color',h.Color); hold on;
% h=legend('r','t'); h.AutoUpdate=false;
title('Strain (-)','Interpreter','tex');
set(gca,'YDir','reverse');
subplot(1,4,3); % temperature
plot(T_last,Ro-grid_r); hold on; title('T (K)','Interpreter','tex'); set(gca,'YDir','reverse');
subplot(1,4,4); % radial displacement (u)
plot(ur_last,Ro-grid_r); hold on; title('u_r');
set(gca,'YDir','reverse');
last_plot_time = 0;

fig1a.h = figure(); % Nimmo's Figure 1a
[ax,h1,h2]=plotyy((Ro-grid_r)/1e3,sigma_t_last/1e6,(Ro-grid_r)/1e3,T_last);
fig1a.ax = ax;
h2.Color = h1.Color;
h2.LineStyle = '--';
hold(ax(1)); hold(ax(2));
set(ax,'Xlim',[0 20]);
set(ax(1),'YLim',[-10 40]);
set(ax(1),'YTick',[-10:5:40]);
set(ax(2),'YTick',[100:20:180]);
ylabel(ax(1),'Tangential Stress (MPa)');
xlabel(ax(1),'Depth (km)');
ylabel(ax(2),'Temperature (K)');
set(ax(2),'YLim',[100 180]);

time=0; itime=1;
% save initial state
isave = 1;
sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
time_store(isave) = time;
last_store = time; isave = isave+1;


while time < t_end
    dt = times(itime+1)-times(itime); itime = itime+1;
    % In each timestep, we do the following
    % 1. Calculate the amount of basal freeze-on and advance the mesh
    % 2. Solve the heat equation using an implicit method
    % 3. Solve for sigma_r
    % 4. Calculate sigma_t
    % 5. Calculate the radial displacements u(r)
    
    % 1. form discrete operators and solve the heat equation
    %     L = zeros(nr,nr);
    %     R = zeros(nr,1);
    %     for i=1:nr
    %         r = grid_r(i);
    %         if i==1
    %             drm = grid_r(i+1)-grid_r(i);
    %         else
    %             drm = grid_r(i)-grid_r(i-1);
    %         end
    %         if i==nr
    %             drp = drm;
    %         else
    %             drp = grid_r(i+1)-grid_r(i);
    %         end
    %         rA = r + drp/2;
    %         rB = r - drm/2;
    %         kA = k;% thermal conductivities
    %         kB = k;
    %         dr = rA-rB;
    %         coef_plus = -kA*rA^2/r^2/drp/dr;
    %         coef_center = rho_i*Cp/dt + kA*rA^2/r^2/drp/dr + kB*rB^2/r^2/drm/dr;
    %         coef_minus = -kB*rB^2/r^2/drm/dr;
    %         R(i) = rho_i*Cp/dt*T_last(i) + H;
    %         L(i,i) =  coef_center;
    %         if( i==1 )
    %             L(i,i+1) = coef_plus-coef_minus;
    %             R(i) = R(i) - 2*Tb*coef_minus;
    %         elseif i==nr
    %             L(i,i-1) = coef_minus-coef_plus;
    %             R(i) = R(i) - 2*Ts*coef_plus;
    %         else
    %             L(i,i-1) = coef_minus;
    %             L(i,i+1) = coef_plus;
    %         end
    %     end
    %     T = L\R;
    
    % 2. Calculate thickening
    % basal heat flux
    %     qb = -k*(T(2)-T(1))/(grid_r(2)-grid_r(1)); % W/m^2
    % thickening would be dx/dt = qb/(L*rho_i)
    %     delta_rb = dt*qb/Lf/rho_i;
    %     z = z_last + delta_rb;
    %     new_thickness = zm(initial_cooling_age+time+dt);
    new_thickness = initial_thickness;
    z=new_thickness-initial_thickness;
    delta_rb = new_thickness-initial_thickness;
    f = (1-2*initial_thickness/Ro);
    delta_rho = (rho_w-rho_i);
    surface_uplift = (delta_rb)*(delta_rho)/rho_w*f/(1+(f*(delta_rho)/rho_w));
       
    % calculate d/dr(Tdot)
    %     Tdot = (T-T_last)/dt;
    Tdot = zeros(size(T_last));
    Tdot(:) = dTdt(Ro-grid_r,initial_cooling_age+time+dt);
    T = T_last + Tdot*dt;
    
    % Timoshenko and Goodier problem 136
    T(:) = Ts + (Tb-Ts)*Ri/(Ro-Ri)*(Ro./grid_r-1);
    Tdot = (T-T_last)/dt;
    
    
    dTdotdr = zeros(nr,1);
    for i=2:nr-1
        dTdotdr(i) = (Tdot(i+1)-Tdot(i))/(grid_r(i+1)-grid_r(i));
    end
    dTdotdr(1) = (Tdot(2)-Tdot(1))/(grid_r(i+1)-grid_r(i));
    dTdotdr(nr) = (Tdot(nr)-Tdot(nr-1))/(grid_r(nr)-grid_r(nr-1));
    
    
    
    % 3. Nonlinear loop over pressure.
    % because the ocean pressure depends on the uplift, we make a guess
    % (above). Using this guess, we calculate stresses, strains, and
    % displacements. Then we re-calculate the pressure using the new value
    % of radial displacement. We continue until the pressure used in the
    % calculations has converged to the pressure consistent with the
    % calculated displacement;
    for iter=1:1%maxiter
        Pex = Pex_last;
        % 3a. Assemble and solve equations for radial stress

        % calculate viscosity at each node (needed by solver)
        mu_node = zeros(nr,1);
        for i=1:nr
            mu_node(i) = mu(T(i));
        end
        [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,Pex,E,nu,dt); 
                                               
        % 5. Calculate the strains        
        dT = T-T_last;
        %         dT(1) = delta_rb*dTdr_b;
        
        dsigma_t = sigma_t - sigma_t_last;
        dsigma_r = sigma_r - sigma_r_last;
        
        de_t = 1/E*(dsigma_t-nu*(dsigma_t+dsigma_r)) -alpha_l*dT + dt/2*(sigma_tD./mu_node); % change in tangential strain
        de_r = 1/E*(dsigma_r-2*nu*dsigma_t)          -alpha_l*dT + dt/2*(sigma_rD./mu_node); % HS91 equations A5-6
        er = er_last + de_r;
        et = et_last + de_t;
        ur = grid_r'.*et; %radial displacement
    end%end nonlinear loop
    
    % 6. advance to next time step and plot (if needed)
    sigma_r_last = sigma_r;
    sigma_t_last = sigma_t;
    T_last = T;
    er_last = er;
    et_last = et;
    z_last = z;
    Pex_last = Pex;
    
    time = time + dt;
    %     if( time-last_plot_time >= plot_interval)
    if (time >= plot_times(iplot) || time >= t_end)
        iplot = iplot+1;
        
        figure(hf);
        subplot(1,4,1);
        h=plot(sigma_r,Ro-grid_r);
        plot(sigma_t,Ro-grid_r,'--','Color',h.Color);
        subplot(1,4,2);
        h=plot(er,Ro-grid_r);
        plot(et,Ro-grid_r,'--','Color',h.Color);
        subplot(1,4,3);
        plot(T,Ro-grid_r);
        subplot(1,4,4);
        plot(ur,Ro-grid_r);
        
        last_plot_time = time;
        drawnow();
    end
end
%% plot numerical and analytic solutions
b=Ro;
a=Ri;
Ti = Tb-Ts;
sig_r_analytic = alpha_l*E*Ti/(1-nu)*a*b/(b^3-a^3)*(a+b-1./grid_r .* (b^2+a*b+a^2) + a^2*b^2./(grid_r.^3) );
sig_t_analytic = alpha_l*E*Ti/(1-nu)*a*b/(b^3-a^3)*(a+b-1./grid_r/2 .* (b^2+a*b+a^2) - a^2*b^2./(2*grid_r.^3) );
% sig_r_analytic = - a^3/(b^3-a^3)*(b^3./grid_r.^3-1)*(-1e6);
% sig_t_analytic =   a^3/(b^3-a^3)*(b^3./(2*grid_r.^3)+1)*(-1e6);
figure, subplot(1,3,1);
plot(sig_r_analytic,Ro-grid_r,'DisplayName','analytic');
hold on
plot(sigma_r,Ro-grid_r,'--','DisplayName','numerical');
legend()
ylabel('Depth (km)');
xlabel('Radial Stress (Pa)')
set(gca,'YDir','reverse');
subplot(1,3,2)
plot(sig_t_analytic,(Ro-grid_r)/1e3,'DisplayName','analytic');
hold on
plot(sigma_t,(Ro-grid_r)/1e3,'--','DisplayName','numerical');
legend()
set(gca,'YDir','reverse');
xlabel('Depth (km)');
ylabel('Tangential Stress (Pa)')
subplot(1,3,3);
plot(ur,(Ro-grid_r)/1e3);
set(gca,'YDir','reverse');
xlabel('Radial displacement');