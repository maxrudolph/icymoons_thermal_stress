% Script to solve coupled ice shell thermal and stress evolution
% Max Rudolph, March 19, 2020
clear;
close all;

% Numerical parameters
nr = 200; % number of grid points
relaxation_parameter=.01; % used in nonlinear loop.
maxiter=100;
% Define physical constants and parameters
% Physical constants
seconds_in_year = 3.15e7;
R=8.314e-3;     % in kJ/mol/K
% Boundary conditions and internal heating
H=0; % internal heating rate.
Tb=270;
Ts=100;
Ro = 2.52e5;             % outer radius of ice shell (m)
Ri = 2.52e5-3e4;         % inner radius of ice shell (m)
Rc = 1.60e5;             % core radius (m)
% Elastic and Viscous properties
E = 5e9;        % shear modulus of ice (Pa)
nu = 0.35;      % Poisson ratio of ice (-)
beta_w = 4e-10; % Compressibility of water (1/Pa)
alpha_l = 1e-4; % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
rho_i=917;      % density of ice (kg/m^3)
rho_w=1000;     % density of water (kg/m^3)
Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=1e14;       % basal viscosity (Pa-s)
mu = @(T) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
% Thermal properties
Cp = 2100; %heat capacity of ice J/kg/K
Lf = 334*1000; % latent heat of fusion (J/kg)
kappa = 1e-6;% m/s/s
k=kappa*rho_i*Cp;

% calculate maxwell time at 100, 270
fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(100)/E,mu(Tb)/E);
fprintf('Thermal diffusion timescale %.2e\n',(4e4)^2/kappa);
% set end time and grid resolution
t_end = 1e15;
dt = 1e11; % time step in seconds
plot_interval = 10*dt;

% set up the grid
grid_r = linspace(Ri,Ro,nr); % set up the grid

% initialize solution vectors (IC)
sigma_r_last = zeros(nr,1); % initial stresses
T_last = linspace(Tb,Ts,nr)';
er_last = zeros(nr,1); % strains
et_last = zeros(nr,1);
ur_last = zeros(nr,1);      % displacement
z_last = 0; % total amount of thickening

% Set up plot
figure(1);
subplot(1,4,1); % sigma_r and sigma_t
h=plot(sigma_r_last,grid_r); hold on;
plot(sigma_r_last,grid_r,'--','Color',h.Color);
% h=legend('\sigma_r','\sigma_t','Interpreter','tex'); h.AutoUpdate=false;
title('Stress (Pa)','Interpreter','tex');
ylabel('r (m)');
subplot(1,4,2); % e_r and e_t
h=plot( sigma_r_last,grid_r); hold on;
plot( sigma_r_last,grid_r,'--','Color',h.Color); hold on;
% h=legend('r','t'); h.AutoUpdate=false;
title('Strain (-)','Interpreter','tex');
subplot(1,4,3); % temperature
plot(T_last,grid_r); hold on; title('T (K)','Interpreter','tex');
subplot(1,4,4); % radial displacement (u)
plot(ur_last,grid_r); hold on; title('u_r');
last_plot_time = 0;

time=0;
while time < t_end
    % In each timestep, we do the following
    % 1. Solve the heat equation using an implicit method
    % 2. Calculate the amount of basal freeze-on
    % 3. Solve for sigma_r
    % 4. Calculate sigma_t
    % 5. Calculate the radial displacements u(r)
    
    % 1. form discrete operators and solve the heat equation
    L = zeros(nr,nr);
    R = zeros(nr,1);
    for i=1:nr
        r = grid_r(i);
        if i==1
            drm = grid_r(i+1)-grid_r(i);
        else
            drm = grid_r(i)-grid_r(i-1);
        end
        if i==nr
            drp = drm;
        else
            drp = grid_r(i+1)-grid_r(i);
        end
        rA = r + drp/2;
        rB = r - drm/2;
        kA = k;% thermal conductivities
        kB = k;
        dr = rA-rB;
        coef_plus = -kA*rA^2/r^2/drp/dr;
        coef_center = rho_i*Cp/dt + kA*rA^2/r^2/drp/dr + kB*rB^2/r^2/drm/dr;
        coef_minus = -kB*rB^2/r^2/drm/dr;
        R(i) = rho_i*Cp/dt*T_last(i) + H;
        L(i,i) =  coef_center;
        if( i==1 )
            L(i,i+1) = coef_plus-coef_minus;
            R(i) = R(i) - 2*Tb*coef_minus;
        elseif i==nr
            L(i,i-1) = coef_minus-coef_plus;
            R(i) = R(i) - 2*Ts*coef_plus;
        else
            L(i,i-1) = coef_minus;
            L(i,i+1) = coef_plus;
        end
    end
    T = L\R;
    
    % 2. Calculate thickening
    % basal heat flux
    qb = -k*(T(2)-T(1))/(grid_r(2)-grid_r(1)); % W/m^2
    % thickening would be dx/dt = qb/(L*rho_i)
    delta_rb = dt*qb/Lf/rho_i;
    z = z_last + delta_rb;
    % mass removal from ocean
    delta_m = dt*qb/Lf;
    % calculate new ocean pressure (Manga and Wang 2007, equation 5)
    Pex_pred = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur_last(1));% ur_last because we don't yet know the uplift
    
    % re-mesh onto new grid
    new_grid_r = linspace(Ri-z,Ro,nr);
    interp_r = [new_grid_r(1) grid_r];
    T = interp1(interp_r,[Tb; T],new_grid_r)';
    T_last = interp1(interp_r,[Tb; T_last],new_grid_r)';
    sigma_r_last = interp1(interp_r,[sigma_r_last(1); sigma_r_last],new_grid_r)';
    er_last = interp1(interp_r,[er_last(1); er_last],new_grid_r)';
    et_last = interp1(interp_r,[et_last(1); et_last],new_grid_r)';
    ur_last = interp1(interp_r,[ur_last(1); ur_last],new_grid_r)';
    grid_r = new_grid_r; % end interpolation step
    
    % calculate d/dr(Tdot)
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
    for iter=1:maxiter
        if iter>1
            Pex = Pex + 0.01*(Pex_post-Pex);
        else
            Pex = Pex_pred;
        end
        % 3a. Assemble and solve equations for radial stress
        M1 = zeros(nr,nr); % coefficients on (dsigma/dt)
        M2 = zeros(nr,nr); % coefficients on (sigma_r)
        R = zeros(nr,1);
        % set ocean pressure
        P = -Pex;
        %     if time == 0
        %         P=0;
        %         dPdt = 0;
        %     else
        %         P=1e6;
        %         dPdt = 0;
        %     end
        
        for i=1:nr
            if i==1
                drm = grid_r(i+1)-grid_r(i);
            else
                drm = grid_r(i) - grid_r(i-1);
            end
            if i==nr
                drp = grid_r(i)   - grid_r(i-1);
            else
                drp = grid_r(i+1) - grid_r(i);
            end
            rA = grid_r(i)+drp/2; % half a cell +
            rB = grid_r(i)-drm/2;   % half a cell -
            drc = rA-rB;
            this_mu = mu( T(i) ); % viscosity
            % M1 - coefficients of dsigma_r/dt
            const1 = (2-nu)/E; % coefficient on d/dr term
            const2 = (1-nu)/E; % coefficient on d(r/2 d/dr)/dr term
            coef_a = rA/2/drp/drc;
            coef_b = rB/2/drm/drc;
            coef_plus   =1/dt*(  const1/drp/2    + const2*coef_a);
            coef_center =1/dt*(                   -const2*coef_a - const2*coef_b);
            coef_minus  =1/dt*( -const1/drm/2    + const2*coef_b);
            if i>1
                M1(i,i-1) = coef_minus;
            end
            M1(i,i)   = coef_center;
            if i<nr
                M1(i,i+1) = coef_plus;
            end
            %BCs
            if i==1
                M1(i,i+1) = M1(i,i+1) - coef_minus;
                R(i) = R(i) - coef_minus*2*P; % sigma_r = P at base
            elseif i==nr
                M1(i,i-1) = M1(i,i-1) - coef_plus;
                R(i) = R(i) - 2*coef_minus*0; % surface sigma_r=0
            end
            % M2 - coefficients of sigma_r
            if i==nr
                TA = T(i);
            else
                TA = (T(i+1)+T(i))/2;
            end
            if i==1
                TB=T(i);
            else
                TB = (T(i-1)+T(i))/2;
            end
            mu_A = mu(TA); % viscosity halfway between i,i+1
            mu_B = mu(TB);
            coef_plus  = 1/(4*this_mu)/(2*drc) + rA/12/mu_A/drp/drc;
            coef_center =                       -rA/12/mu_A/drp/drc - rB/12/mu_B/drm/drc;
            coef_minus =-1/(4*this_mu)/(2*drc) + rB/12/mu_B/drm/drc;
            if i>1
                M2(i,i-1) = coef_minus;
            end
            M2(i,i)   = coef_center;
            if i<nr
                M2(i,i+1) = coef_plus;
            end
            if i==1
                R(i) = R(i) - coef_minus*2*P;
                M2(i,i+1) = M2(i,i+1) - coef_minus;
            elseif i==nr
                M2(i,i-1) = M2(i,i-1) - coef_plus;
                R(i) = R(i) - coef_plus*0; % no change because sigma_r = 0 at surface
            end
            R(i) = R(i) + alpha_l*dTdotdr(i);
            %         R(i) = R(i)+alpha_l*(Tdot(i+1)-Tdot(i))/2/drc; % this term
            %         includes the coupling to the energy equation - Tdot needs
            %         to be updated
        end
        LHS = (M1+M2);
        RHS = (R+M1*sigma_r_last);
        sigma_r = LHS\RHS;
        
        % 4. calculate the tangential stress sigma_t
        % first, calculate dsr/dr
        dsrdr = zeros(size(sigma_r));
        for i=2:nr-1
            dr = grid_r(i+1)-grid_r(i-1);
            dsrdr(i) = (sigma_r(i+1)-sigma_r(i-1))/dr;
        end
        dsrdr(1) = (sigma_r(2)-sigma_r(1))/( grid_r(2)-grid_r(1));
        dsrdr(nr) = (sigma_r(nr)-sigma_r(nr-1))/(grid_r(nr)-grid_r(nr-1));
        
        sigma_t = sigma_r + (grid_r'/2).*dsrdr;
        
        % 5. Calculate the strains
        sigma_tD =  grid_r'/6.*dsrdr;     %deviatoric stresses, from Hillier and Squyres (1991) equations A8-9
        sigma_rD = -grid_r'/3.*dsrdr;
        % calculate viscosity at each node
        mu_node = zeros(nr,1);
        for i=1:nr
            mu_node(i) = mu(T(i));
        end
        dT = T-T_last;
        de_t = 1/E*(sigma_t-nu*(sigma_t+sigma_r))-alpha_l*dT + dt/2*(sigma_tD./mu_node); % change in tangential strain
        de_r = 1/E*(sigma_r-2*nu*sigma_t)-alpha_l*dT + dt/2*(sigma_rD./mu_node);         % HS91 equations A5-6
        er = er_last + de_r;
        et = et_last + de_t;
        ur = grid_r'.*et; %radial displacement
        % re-calculate excess pressure using new uplift
        Pex_post = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur(1));% ur_last because we don't yet know the uplift
        
        % check for convergence
        if abs( Pex_post-Pex )/abs(Pex) < 1e-2
            % fprintf('end of nonlinear loop. Pex_post %.2e Pex %.2e\n',Pex_post,Pex);
            fprintf('converged in %d iterations\n',iter);
            break;
        elseif iter==maxiter
            warning('Nonlinear loop failed to converge');
        end
    end%end nonlinear loop
    
    % 6. advance to next time step and plot (if needed)
    sigma_r_last = sigma_r;
    T_last = T;
    er_last = er;
    et_last = et;
    z_last = z;
    
    time = time + dt;
    if( time-last_plot_time > plot_interval)
        figure(1);
        subplot(1,4,1);
        h=plot(sigma_r,grid_r);
        plot(sigma_t,grid_r,'--','Color',h.Color);
        subplot(1,4,2);
        h=plot(er,grid_r);
        plot(et,grid_r,'--','Color',h.Color);
        subplot(1,4,3);
        plot(T,grid_r);
        subplot(1,4,4);
        plot(ur,grid_r);
        last_plot_time = time;
    end
    % stop criterion on sigma_t large
    if any(sigma_t > 3e6 )
        break;
    end
    
    
end