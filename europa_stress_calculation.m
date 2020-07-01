% Script to solve coupled ice shell thermal and stress evolution
% Max Rudolph, March 19, 2020
clear;
close all;
addpath core; % this is where the helper functions live.

% Numerical parameters
nr = 500; % number of grid points
relaxation_parameter=1e-4; % used in nonlinear loop.
maxiter=300;
% Define physical constants and parameters
% Physical constants
seconds_in_year = 3.1558e7;
R=8.314e-3;     % in kJ/mol/K
% Boundary conditions and internal heating
H=0; % internal heating rate.
Tb=270;
Ts=100;
Ro = 1.561e6;             % outer radius of ice shell (m)
Ri = Ro-2.4e3;         % inner radius of ice shell (m)
Rc = Ro-1.3e5;             % core radius (m)
% Elastic and Viscous properties
E = 5e9;        % shear modulus of ice (Pa)
nu = 0.3;      % Poisson ratio of ice (-)
beta_w = 4e-10; % Compressibility of water (1/Pa)
alpha_l = 1e-4; % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
rho_i=900;      % density of ice (kg/m^3)
rho_w=1000;     % density of water (kg/m^3)
Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=1e15;       % basal viscosity (Pa-s)
mu = @(T) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
% Failure criterion:
g = 1.3;        % used to plot a failure curve
tensile_strength = 3e6; % tensile strength, Pa
% Thermal properties
Cp = 2100; %heat capacity of ice J/kg/K
Lf = 334*1000; % latent heat of fusion (J/kg)
kappa = 1e-6;% m/s/s
k=kappa*rho_i*Cp;

% calculate maxwell time at 100, 270
fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(100)/E,mu(Tb)/E);
fprintf('Thermal diffusion timescale %.2e\n',(4e4)^2/kappa);
% set end time and grid resolution
t_end = 0.1e6*seconds_in_year;
% dt = 1e4*seconds_in_year; % time step in seconds
dtmax = 1e2*seconds_in_year;
dtmin = seconds_in_year;
% dt1 = 3600; % size of first timestep
% times = logspace(log10(dt1),log10(t_end+dt1),1e4)-dt1;
plot_interval = t_end;
save_interval = 1e5*seconds_in_year;
save_depths = [0 0.5 1 1.5 2 2.5]*1000;
nsave = ceil(t_end/save_interval) + 1; nsave_depths = length(save_depths);
sigma_t_store = zeros(nsave_depths,nsave);

% set up the grid
grid_r = linspace(Ri,Ro,nr); % set up the grid

% initialize solution vectors (IC)
sigma_r_last = zeros(nr,1); % initial stresses
sigma_t_last = zeros(nr,1); % initial stresses
T_last = linspace(Tb,Ts,nr)';
er_last = zeros(nr,1); % strains
et_last = zeros(nr,1);
ur_last = zeros(nr,1);      % displacement
z_last = 0; % total amount of thickening
Pex_last = 0; %initial overpressure

% Set up plot
hf2=figure();

% plot_times = [0.0 0.1 0.2 0.3 0.4 0.5]*1e6*seconds_in_year; iplot=2;
plot_times = linspace(0,t_end,5); iplot=2;
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
subplot(2,1,1);
[ax,h1,h2]=plotyy((Ro-grid_r)/1e3,sigma_t_last/1e6,(Ro-grid_r)/1e3,T_last);
fig1a.ax = ax;
h2.Color = h1.Color;
h2.LineStyle = '--';
hold(ax(1)); hold(ax(2));
set(ax,'Xlim',[0 10]);
set(ax(1),'YLim',[-10 40]);
set(ax(1),'YTick',[-10:5:40]);
set(ax(2),'YTick',[100:20:180]);
set(ax(1),'YTickLabelMode','auto');
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

failure_occurred = 0;
while time < t_end && ~failure_occurred
    % In each timestep, we do the following
    % 1. Calculate the amount of basal freeze-on and advance the mesh
    % 2. Solve the heat equation using an implicit method
    % 3. Solve for sigma_r
    % 4. Calculate sigma_t
    % 5. Calculate the radial displacements u(r)
    
    % 1. Calculate basal freeze-on and interpolate old solution onto new mesh
    % calculate heat flux
    Tg = Tb-(T_last(2)-Tb);
    dTdr_b_last = (T_last(2)-Tg)/2/(grid_r(2)-grid_r(1));
    qb = -k*dTdr_b_last;
    
    % determine the timestep
    dt = dtmax;
    if abs(qb/Lf/rho_i*dt) > (grid_r(2)-grid_r(1))/2
       dt =  (grid_r(2)-grid_r(1))/2/(qb/Lf/rho_i);
    end
    if dt < dtmin
        dt = dtmin;
        disp('Setting dt = dtmin');
    end
    % thickening would be dx/dt = qb/(L*rho_i)    
    delta_rb = dt*qb/Lf/rho_i;
    z = z_last + delta_rb;
    
    % calculate new ocean pressure (Manga and Wang 2007, equation 5)
    Pex_pred = Pex_last + 3*Ri^2/beta_w/(Ri^3-Rc^3)*(delta_rb*(rho_w-rho_i)/rho_w-ur_last(1)); % ur_last because we don't yet know the uplift
    
    % re-mesh onto new grid
    new_grid_r = linspace(Ri-z,Ro,nr);
    [T_last,sigma_r_last,sigma_t_last,er_last,et_last] = interpolate_solution(new_grid_r,grid_r,T_last,sigma_r_last,sigma_t_last,er_last,et_last,Tb);
    grid_r = new_grid_r; % end interpolation step
    
    % 2. form discrete operators and solve the heat equation
    T = solve_temperature_shell(grid_r,T_last,Tb,Ts,k,rho_i,Cp,H,dt);
    
    %Compute d(Tdot)/dr
    Tdot = (T-T_last)/dt;
    dTdr_b = (T(2)-T(1))/(grid_r(2)-grid_r(1));
    Tdot(1) = delta_rb*dTdr_b/dt; % this is an approximation to the Eulerian cooling rate at the ocean-ice interface
    
    dTdotdr = zeros(nr,1);
    for i=2:nr-1
        dTdotdr(i) = (Tdot(i+1)-Tdot(i-1))/(grid_r(i+1)-grid_r(i-1));
    end
    dTdotdr(1) = (Tdot(2)-Tdot(1))/(grid_r(2)-grid_r(1));
    dTdotdr(nr) = (Tdot(nr)-Tdot(nr-1))/(grid_r(nr)-grid_r(nr-1));
    
    % 3. Nonlinear loop over pressure.
    % because the ocean pressure depends on the uplift, we make a guess
    % (above). Using this guess, we calculate stresses, strains, and
    % displacements. Then we re-calculate the pressure using the new value
    % of radial displacement. We continue until the pressure used in the
    % calculations has converged to the pressure consistent with the
    % calculated displacement;
    for iter=1:maxiter
        if iter>1
            Pex = Pex + relaxation_parameter*(Pex_post-Pex);
        else
            Pex = Pex_last;
        end
        
        % calculate viscosity at each node
        mu_node = zeros(nr,1);
        for i=1:nr
            mu_node(i) = mu(T(i));
        end
        
        [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-Pex,E,nu,dt);
               
        % 5. Calculate the strains        
        dT = T-T_last;
        dT(1) = delta_rb*dTdr_b;
        
        dsigma_t = sigma_t - sigma_t_last;
        dsigma_r = sigma_r - sigma_r_last;
        
        de_t = 1/E*(dsigma_t-nu*(dsigma_t+dsigma_r))+alpha_l*dT + dt/2*(sigma_tD./mu_node); % change in tangential strain
        de_r = 1/E*(dsigma_r-2*nu*dsigma_t)         +alpha_l*dT + dt/2*(sigma_rD./mu_node); % HS91 equations A5-6
        er = er_last + de_r;
        et = et_last + de_t;
        ur = grid_r'.*et; %radial displacement
        % re-calculate excess pressure using new uplift
        Pex_post = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur(1));% ur_last because we don't yet know the uplift
        %fprintf('iter %d. Pex_post %.2e Pex %.2e\n',iter,Pex_post,Pex);
        
        % check for convergence
        if abs( Pex_post-Pex )/abs(Pex) < 1e-6
            fprintf('dt=%.2e yr, time=%.3e Myr, Pex_post %.6e Pex %.6e, converged in %d iterations\n',dt/seconds_in_year,(time+dt)/seconds_in_year/1e6,Pex_post,Pex,iter);            
            break;
        elseif iter==maxiter
            error('Nonlinear loop failed to converge');
        end
    end%end nonlinear loop
    % 5. Determine whether tensile failure has occurred
    failure = tensile_failure_criterion(Ro-grid_r',sigma_t,rho_i,g,tensile_strength);
    if( any(failure))
        disp(['Failure criterion has been reached']);
        idx_shallow = find(failure,1,'last');
        idx_deep = find(failure,1,'first');
        fprintf('Shallowest, deepest failure: %f, %f\n\n',Ro-grid_r(idx_shallow),Ro-grid_r(idx_deep));
        fprintf('Failure time: %f Myr\n',time / seconds_in_year / 1e6);
        fprintf('Surface stress at failure: %f MPa\n',sigma_t(end)/1e6);
        failure_occurred = failure_occurred + 1;
    end
    
    % 6. advance to next time step and plot (if needed)
    sigma_r_last = sigma_r;
    sigma_t_last = sigma_t;
    T_last = T;
    er_last = er;
    et_last = et;
    z_last = z;
    ur_last = ur;
    Pex_last = Pex;
    
    time = time + dt;
    
    if (time >= plot_times(iplot) || time >= t_end || failure_occurred)
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
        
        
        figure(hf2);
        plot(ur(end),sigma_t(end),'.'); hold on;
        
        figure(fig1a.h); % Nimmo's Figure 1a
        h=plot(fig1a.ax(1),(Ro-grid_r)/1e3,sigma_t_last/1e6);
        plot(fig1a.ax(2),(Ro-grid_r)/1e3,T_last,'--','Color',h.Color);
        
        
        last_plot_time = time;
        drawnow();
    end
    if (time-last_store >= save_interval || time >= t_end || failure_occurred)
        sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
        time_store(isave) = time;
        last_store = time; isave = isave+1;       
    end
end
%% add legends
for i=1:length(plot_times)
    labels{i} = sprintf('%.2f Myr',plot_times(i)/seconds_in_year/1e6);
end
figure( fig1a.h );
axis(fig1a.ax(1));
legend(labels,'Location','southeast','AutoUpdate','off');
plot(fig1a.ax(1),(Ro-grid_r)/1e3,(tensile_strength + rho_i*g*(Ro-grid_r))/1e6,'k');

%% Nimmo's figure 1b
figure(fig1a.h);
labels = {};
for i=1:length(save_depths)
    labels{i} = sprintf('%.1f km',save_depths(i)/1000);
end
subplot(2,1,2);
plot(time_store(1:isave-1)/seconds_in_year/1e6,sigma_t_store(:,1:isave-1)/1e6);
legend(labels,'Location','southeast');
xlabel('Time (Myr)');
ylabel('Tangential Stress (MPa)');
set(gcf,'Color','w');
h=gcf;
h.Position(3:4) = [390 580];
saveas(gcf,'Europa_initial2.4km.eps','psc2')

%% Plastic yielding?
figure();
Pnhs = -1/3*(sigma_r + 2*sigma_t); % compression = positive for pressure
Phydrostatic = rho_i*g*(Ro-grid_r');
Ptotal = Pnhs + Phydrostatic;
Sr = sigma_r - 1/3*(sigma_r + 2*sigma_t);
St = sigma_t - 1/3*(sigma_r + 2*sigma_t);
J2 = 1/2*(sigma_rD.^2 + 2*sigma_tD.^2);
sig_eff = sqrt(3*J2);

tau_yield = 3e6 + 0.0*Ptotal;
subplot(2,1,1);
plot(Ro-grid_r,Phydrostatic); hold on
plot(Ro-grid_r,Pnhs,'r--');
plot(Ro-grid_r,Pnhs+Phydrostatic,'k');
legend('Hydrostatic','Nonhydrostatic','Total');
subplot(2,1,2);
plot(Ro-grid_r,tau_yield,'g'); hold on;
plot(Ro-grid_r,sig_eff,'k');
legend('Yield Stress','Effective Stress');


