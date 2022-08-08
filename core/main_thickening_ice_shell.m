function result = main_thickening_ice_shell(p)
% This function solves coupled thermal and viscoelastic evolution of
% icy moons. This version incorporates the effect of ammonia (NH3)
% p is a data structure containing parameters
% Max Rudolph, March 19, 2020
% adapted for other moons by Alyssa Rhoden, 6/2021


Ro=p.Ro;
Ri = p.Ri;
Rc = p.Rc;
max_depth = Ro-Rc;
g = p.g;
Ts= p.Ts;

Qbelow = p.Qbelow;
relaxation_parameter=p.relaxation_parameter;
X0 = p.X0; % initial ammonia content.
label = p.label;

viscosity_model = 0; % 0 = Nimmo (2004), 1 = Goldsby and Kohlstedt (2001)
viscosity.d = 1e-3; % grain size in m used to calculate the viscosity
viscosity.P = 1e5; % Pressure in MPa used to calculate the viscosity

t_end = p.t_end;
nr = 512; % number of grid points in radial direction.

% label for runs.
if viscosity_model == 0
    label = [label '-nimmovisc'];
elseif viscosity_model == 1
    label = [label '-goldsbykohlstedt'];
else
    error('not implemented');
end

ifail = 1; % index into list of times at which failure occurred.

maxiter=10000;
% Define physical constants and parameters
% Physical constants
seconds_in_year = 3.1558e7;
R=8.314e-3;     % in kJ/mol/K
% Boundary conditions and internal heating
H=0; % internal heating rate.
Tb = ammonia_melting(X0);

% Elastic and Viscous properties of the ice shell
E = 5e9;        % shear modulus of ice (Pa)
nu = 0.3;       % Poisson ratio of ice (-)
beta_w = 4e-10; % Compressibility of water (1/Pa)
alpha_l = 3e-5; % coefficient of linear thermal expansion of ice ( alpha_v/3 ) (1/K)
rho_i=917;      % density of ice (kg/m^3)
%                 rho_w=1000;     % density of water (kg/m^3)
rho_w = ammonia_density(X0,rho_i*g*(Ro-Ri),Tb);
Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=p.mub;       % basal (273K) viscosity (Pa-s)
if viscosity_model == 0
    mu = @(T,stress) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
elseif viscosity_model == 1
    mu = @(T,stress) goldsby_kohlstedt(stress,T,viscosity.d,viscosity.P); % Goldsby-Kohlstedt effective viscosity
end
% Failure criterion:
tensile_strength = 3e6; % tensile strength, Pa
cohesion = 2e7;  % plastic yield strength
friction = 0.6; % friction angle for plastic yielding
% Thermal properties
Cp = 2100; %heat capacity of ice J/kg/K
Lf = 334*1000; % latent heat of fusion (J/kg)
kappa = 1e-6;% m/s/s
%k=kappa*rho_i*Cp;
k = @(T) 651./T;
%
% Basal heating model - depends on thickness and transport properties
%
%Q0 = k*(Tb-Ts)/(Ro-Ri);% time-averaged basal heat flux
[Q0,T,q] = find_steady_T(Ri,Ro,Tb,Ts,linspace(Ri,Ro,nr));
%perturbation_period = 1.0e8*seconds_in_year;
%deltaQonQ = 1.0; % fractional perturbation to Q0.

% calculate maxwell time at 100, 270
fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(100,0)/E,mu(Tb,0)/E);
fprintf('Thermal diffusion timescale %.2e\n',(4e4)^2/kappa);
% set end time and grid resolution

dtmax = 1e6*seconds_in_year;
dtmin = 3600;%*seconds_in_year;
% dt1 = 3600; % size of first timestep
% times = logspace(log10(dt1),log10(t_end+dt1),1e4)-dt1;
plot_interval = t_end;
save_interval = 1e4*seconds_in_year;
save_depths = linspace(0,max_depth,500);

nsave = ceil(t_end/save_interval) + 2;
nsave_depths = length(save_depths);
sigma_t_store = zeros(nsave_depths,nsave);

result.time = zeros(nsave,1);
result.z = zeros(nsave,1);
result.Ri = zeros(nsave,1); result.Ri(1) = Ri;
result.qb = zeros(nsave,1);
result.sigma_t = NaN*zeros(nsave_depths,nsave);
result.sigma_r = NaN*zeros(nsave_depths,nsave);
result.Pex = zeros(nsave,1);
result.XNH3 = zeros(nsave,1);
result.dTdr = zeros(nsave_depths,nsave);
result.T = zeros(nsave_depths,nsave);
result.Tb = zeros(nsave,1);
result.ur = zeros(nsave_depths,nsave);
result.failure_time = zeros(1,nsave);
result.failure_P = zeros(1,nsave);
result.failure_Pex_crit = zeros(1,nsave);
result.failure_dP = zeros(1,nsave);
result.failure_thickness = zeros(1,nsave);
result.failure_top = zeros(1,nsave);
result.failure_bottom = zeros(1,nsave);
result.failure_erupted_volume = NaN*zeros(1,nsave);
result.failure_erupted_volume_pressurechange = NaN*zeros(1,nsave);
result.failure_erupted_volume_volumechange = NaN*zeros(1,nsave);
result.failure_sigma_t = cell(1,nsave);
result.failure_sigma_r = cell(1,nsave);
result.failure_r = cell(1,nsave);
result.failure_z = zeros(1,nsave);


erupted_volume = 0;
erupted_volume_pressurechange = 0;
erupted_volume_volumechange = 0;

% set up the grid
grid_r = linspace(Ri,Ro,nr); % set up the grid

% initialize solution vectors (IC)
sigma_r_last = zeros(nr,1); % initial stresses
sigma_t_last = zeros(nr,1); % initial stresses
siiD_last = zeros(nr,1); % deviatoric stress invariant - used for viscosity
T_last = zeros(nr,1);
% Initialize T with steady numerical solution.
%T_last = solve_temperature_shell(grid_r,T_last,Tb,Ts,k,rho_i,Cp,H,Inf,0.0);
T_last(:) = solve_stefan_analytic(grid_r(end)-grid_r,k(Tb),rho_i,Cp,Lf,Tb,Ts);

er_last = zeros(nr,1); % strains
et_last = zeros(nr,1);
ur_last = zeros(nr,1);      % displacement
z_last = 0;    % total amount of thickening
dzdt_last = 0; % thickening rate
Pex_last = 0; %initial overpressure

% Set up plot
hf2=figure();

% plot_times = [0.0 0.1 0.2 0.3 0.4 0.5]*1e6*seconds_in_year; iplot=2;
plot_times = linspace(0,t_end,5); iplot=2;
hf=figure();
subplot(1,4,1); % sigma_r and sigma_t
h=plot(sigma_r_last,Ro-grid_r); hold on;
plot(sigma_t_last,Ro-grid_r,'--','Color',h.Color);
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
X = X0;

failure_mask = false(size(grid_r)); % stores whether failure occurred
failure_time = zeros(size(grid_r)); % stores the time at which failure occurred

while time < t_end && (Ri-z_last > Rc) && X <= 0.32
    % In each timestep, we do the following
    % 1. Calculate the amount of basal freeze-on and advance the mesh
    % 2. Solve the heat equation using an implicit method
    % 3. Solve for sigma_r
    % 4. Calculate sigma_t
    % 5. Calculate the radial displacements u(r)
    
    % 1. Calculate basal freeze-on and interpolate old solution onto new mesh
    % calculate heat flux
    dt = dtmax;
    Tg = Tb-(T_last(2)-Tb);
    dTdr_b_last = (T_last(2)-Tg)/2/(grid_r(2)-grid_r(1));
    qb = -k(Tb)*dTdr_b_last;
    qb_net = qb - Qbelow(time+dt); % first term is conducted heat. second term is heat supplied from below.
    
    % determine the timestep
    if abs(qb_net/Lf/rho_i*dt) > (grid_r(2)-grid_r(1))/2
        dt = abs( (grid_r(2)-grid_r(1))/2/(qb_net/Lf/rho_i) );
    end
    if dt < dtmin
        dt = dtmin;
        disp('Setting dt = dtmin');
    end
    if any(failure_mask)
        dt = dtmin;
    end
    
    % compute the derivative of Tm with respect to X
    Xp = X + 1e-6;
    Xm = X - 1e-6;
    if Xm < 0
        Xm = 0;
    end
    dTmdX = (ammonia_melting(Xp)-ammonia_melting(Xm))/(Xp-Xm);
    Vlast = 4/3*pi*((Ri-z_last)^3 - Rc^3); % compute current ocean volume
    Alast = 4*pi*(Ri-z_last)^2;
    % change in ocean internal energy per unit thickness per unit area
    % (m)*(kg/m^3)*(J/kg/K) = J/m^2/K
    dUdTm = Vlast/Alast * rho_w*ammonia_cp(X);
    % change in ammonia content per unit thickening
    % (m)^2/(m^3) -> 1/m
    dXdz = -4*pi*(Ri-z_last)^2*X/Vlast; % note that dX/dz = X0*V0/(V0+dV)^2. The expression used here is a 1st order approximation.
    % compute an effective latent heat that accounts for heat
    % added/removed from ocean to maintain Tm at base.
    L_eff = dUdTm * dTmdX * dXdz;% J/m^2/K * K * 1/m =>  J/m^3
    
    qb_net = qb - Qbelow(time+dt);
    
    % thickening would be dx/dt = qb/(L*rho_i)
    delta_rb = dt*qb_net/(Lf*rho_i + L_eff);% s * W/(J/kg/K * kg/m^3 + J/m^3)
    z = z_last + delta_rb;
    dV = -4*pi*(Ri-z)^2*delta_rb;% volume change per unit change in thickness
    X = X*Vlast/(Vlast+dV); % update nh3 content (assumes no nh3 in ice)
    % update ocean density
    rho_w = ammonia_density(X,rho_i*g*(grid_r(end)-grid_r(1)),Tb);
    beta_w = ammonia_compressibility(X,rho_i*g*(grid_r(end)-grid_r(1)),Tb);
    
    % compute the melting temperature for the new NH3 content at
    % the ocean-ice interface:
    Tmelt = ammonia_melting(X);
    Tb = Tmelt;
    dzdt = delta_rb/dt;
    
    if (Ri-z-delta_rb <= Rc)
        % code seems to get very unstable when the ocean is too
        % thin...
        break
    end
    % calculate new ocean pressure (Manga and Wang 2007, equation 5)
    Pex_pred = Pex_last + 3*(Ri-z)^2/beta_w/((Ri-z)^3-Rc^3)*(delta_rb*(rho_w-rho_i)/rho_w-ur_last(1)); % ur_last because we don't yet know the uplift
    
    new_grid_r = linspace(Ri-z,Ro,nr);
    dTdr_last = (T_last(2)-Tb)/(grid_r(2)-grid_r(1));
    [T_last,sigma_r_last,sigma_t_last,er_last,et_last] = interpolate_solution(new_grid_r,grid_r,T_last,sigma_r_last,sigma_t_last,er_last,et_last,Tb);
    grid_r = new_grid_r; % end interpolation step
    
    % 2. form discrete operators and solve the heat equation
    [T,dTdotdr] = solve_temperature_shell(grid_r,T_last,Tb,Ts,k,rho_i,Cp,H,dt,delta_rb);
    
    % 3. Nonlinear loop over pressure.
    % because the ocean pressure depends on the uplift, we make a guess
    % (above). Using this guess, we calculate stresses, strains, and
    % displacements. Then we re-calculate the pressure using the new value
    % of radial displacement. We continue until the pressure used in the
    % calculations has converged to the pressure consistent with the
    % calculated displacement;
    converged = false;
    pex_store = zeros(maxiter,1);
    pexpost_store = zeros(maxiter,1);
    for iter=1:maxiter
        if iter>10 && iter < 100
            [tmp,ind] = unique(pex_store(1:iter-1));
            Pex = interp1(pexpost_store(ind)-pex_store(ind),pex_store(ind),0,'linear','extrap');
        elseif iter>1
            Pex = Pex + relaxation_parameter*(Pex_post-Pex);
        else
            Pex = Pex_last;
        end
        
        % compute height to which water would rise
        [ptmp,Ttmp,ztmp,rhotmp] = ammonia_adiabatic_profile(X,g*rho_i*(grid_r(end)-grid_r(1)),g);
        rhobar = cumtrapz(ztmp,rhotmp);
        rhobar = 1/(ztmp(end)-ztmp(1))*rhobar(end);
        
        Pex_crit = (rhobar-rho_i)*(Ro-(Ri-z))*g;
        
        % calculate viscosity at each node
        visc_converged = false;
        visc_iter = 100;
        ivisc = 1;
        while ~visc_converged && ivisc <= visc_iter
            % compute mu for current siiD
            if ivisc == 1
                siiD = siiD_last;
            else
                siiD = siiD_post;
            end
            
            
            mu_node = zeros(nr,1);
            mu_node(:) = mu(T,siiD);
            % reduce Maxwell time in region experiencing failure
            if all(failure_mask)
                if Pex_last >= Pex_crit
                    % Calculate the volume erupted (dP)*beta*V0 + V-V0
                    pressure_contribution = (Pex_last-Pex_crit)*beta_w*(4/3*pi*((Ri-z)^3-Rc^3));
                    urelax = (Ri-z)/E*(1-2*nu)*(Pex_last-Pex_crit); % Manga and Wang (2007) equation 4
                    volume_contribution = (Ri-z)^2*urelax*4*pi; % (4*pi*R^2)*dr
                else
                    pressure_contribution = 0;
                    volume_contribution = 0;
                end
                % reset stresses and uplift
                sigma_r_last = 0*sigma_r_last;
                sigma_t_last = 0*sigma_t_last;
                er_last = 0*er_last;
                et_last = 0*et_last;
                ur_last = 0*ur_last;
                Pex=0; % force zero pressure.
                converged = true;
                Ri = Ri - z;
                z = 0;
                z_last=0;
                % move the inner radius effectively to the current position
                % of the base of the ice shell. Then set the amount of
                % freezing to zero.
            elseif( any(failure_mask) )
                minimum_viscosity_prefactor = 0; % maximum allowable fractional reduction in viscosity
                mu_node(failure_mask) = min(mu_node(failure_mask),max(minimum_viscosity_prefactor*mu_node(failure_mask),0.1*E*dt));  % timestep = 10x(maxwell time)
                above_crack = find( failure_mask,1,'last');
                above_crack_mask = false(size(failure_mask));
                above_crack_mask( above_crack+1 : end ) = true;
                %                 mu_node(above_crack_mask) = min( mu_node(above_crack_mask),100*E*dt ); % limit maximum viscosity to 100*maxwell time
                %                 for i=1:3
                %                 tmp = exp(smooth(log( mu_node )));
                %                 mu_node(~failure_mask) = tmp(~failure_mask);
                %                 end
                %                 mu_node = exp(smooth(log( mu_node )));
                if iter==1
                    Pex=0; % If failure occurs, it's better to guess that all pressure is relieved. Other choices could cause convergence problems.
                end
            end
            
            % Calculate Stresses
            [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-Pex,E,nu,dt);
            siiD_post = sqrt( 0.5*(sigma_rD.^2 + 2*sigma_tD.^2) );
            norm_change = min(norm(siiD_post-siiD)/norm(siiD),norm(siiD_post-siiD));
            %                     disp([num2str(ivisc) ' change in norm of siiD:' num2str(norm_change)]);
            if isnan(norm_change)
                keyboard
            elseif norm_change < 1e-4
                visc_converged = true;
            end
            
            ivisc = ivisc+1;
        end
        if all(failure_mask)
            erupted_volume = erupted_volume + pressure_contribution + volume_contribution;
            erupted_volume_pressurechange = erupted_volume_pressurechange + pressure_contribution;
            erupted_volume_volumechange = erupted_volume_volumechange + volume_contribution;
        end
        
        % 5. Calculate the strains
        dT = T-T_last;
        dr1 = grid_r(2)-grid_r(1);
        dr2 = grid_r(3)-grid_r(1);
        L = [0 0 1;
            dr1^2 dr1 1;
            dr2^2 dr2 1];
        R = T(1:3);
        coef = L\R;
        dTdr_b = coef(2);
        %             dTdr_b=(T(2)-Tb)/(grid_r(2)-grid_r(1));
        dT(1) = delta_rb*dTdr_b;
        
        dsigma_t = sigma_t - sigma_t_last;
        dsigma_r = sigma_r - sigma_r_last;
        %             mu_node(2:end-1) = exp(0.5*(log(mu_node(1:end-2))+log(mu_node(3:end))));
        de_t = 1/E*(dsigma_t-nu*(dsigma_t+dsigma_r))+alpha_l*dT + dt/2*(sigma_tD./mu_node); % change in tangential strain
        de_r = 1/E*(dsigma_r-2*nu*dsigma_t)         +alpha_l*dT + dt/2*(sigma_rD./mu_node); % HS91 equations A5-6
        er = er_last + de_r;
        et = et_last + de_t;
        ur = grid_r'.*et; %radial displacement
        
        ei = 2*de_t + de_r; % first invariant of strain
        de_tD = de_t - 1/3*ei;
        de_rD = de_r - 1/3*ei;
        eiiD = sqrt( 0.5*(de_rD.^2 + 2*de_tD.^2) ); % second invariant of deviatoric strain
        
        % re-calculate excess pressure using new uplift
        %             Pex_post = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur(1));
        Pex_post = Pex_last + 3*(Ri-z)^2/beta_w/((Ri-z)^3-Rc^3)*((z-z_last)*(rho_w-rho_i)/rho_w-(ur(1)-ur_last(1)));
        % Calculate the critical excess presssure necessary to
        % erupt water onto the surface.
        %fprintf('iter %d. Pex_post %.2e Pex %.2e\n',iter,Pex_post,Pex);
        
        % check for convergence
        if abs( Pex_post-Pex )/abs(Pex) < 1e-2 || abs(Pex_post-Pex) < 1e2
            fprintf('dt=%.2e yr, time=%.3e Myr, Pex_post %.6e Pex %.6e, converged in %d iterations\n',dt/seconds_in_year,(time+dt)/seconds_in_year/1e6,Pex_post,Pex,iter);
            converged = true;
        elseif iter==maxiter
            errorstring = ['Nonlinear loop failed to converge. ' sprintf('Ri=%f, X0=%f',p.Ri,p.X0)];
            error(errorstring);%'Nonlinear loop failed to converge');
            
        end
        
        pex_store(iter) = Pex;
        pexpost_store(iter) = Pex_post;
        if converged
            break;
        end
    end%end nonlinear loop
    
    
    if max(abs(diff(ur))) > 100
        % a discontinuity has developed
        %             figure();
        %             plot(1/E*(dsigma_t-nu*(dsigma_t+dsigma_r))); hold on
        %             plot(alpha_l*dT);
        %             plot(dt/2*(sigma_tD./mu_node));
        %             legend('elastic','thermal','viscous');
        %             keyboard
        %             close();
    end
    
    
    % 5. Determine whether tensile failure has occurred
    failure = tensile_failure_criterion(Ro-grid_r',sigma_t,rho_i,g,tensile_strength);
    if(any(failure)) % failure is occurring
        disp(['Failure criterion has been reached']);
        idx_shallow = find(failure,1,'last');
        idx_deep = find(failure,1,'first');
        fprintf('Shallowest, deepest failure: %f, %f\n\n',Ro-grid_r(idx_shallow),Ro-grid_r(idx_deep));
        fprintf('Failure time: %f Myr\n',time / seconds_in_year / 1e6);
        fprintf('Surface stress at failure: %f MPa\n',sigma_t(end)/1e6);
        
        % check to see if a crack could propagate to surface
        % 1. Find the midpoint of the crack
        % 2. Look upward - balance stresses on crack in upward
        % direction
        % 3. If crack reached surface, balance stresses on entire
        % crack. Otherwise balance stresses in downward direction.
        sigma_t_tot = sigma_t - rho_i*g*(Ro-grid_r');
        depth = Ro-grid_r; % depth will be in descending order, i.e. deepest first
        midpoint_depth = mean(depth([idx_shallow idx_deep]));
        [~,midpoint_ind] = max( sigma_t_tot );
        if midpoint_ind == nr
            stress_above = 0;
        else
            stress_above = cumtrapz( grid_r(midpoint_ind:end), sigma_t_tot(midpoint_ind:end) );  % integrate in upward direction
        end
        stress_below = cumtrapz( depth(midpoint_ind:-1:1), sigma_t_tot(midpoint_ind:-1:1) ); % integrate in downward direction
        if stress_above(end) >= 0
            disp('Crack reached surface');
            surface_failure = true;
            above_stress_integral = stress_above(end);
            net_tension = above_stress_integral + stress_below;
            % find depth at which crack stops
            ind = find(net_tension > 0,1,'last'); % net tension is ordered by increasing depth
            depth_tmp = depth(midpoint_ind:-1:1);
            max_depth = depth_tmp(ind); % depth at which crack stops
            min_depth = 0;
            if net_tension > 0
                disp('Crack reaches ocean!');
            end
        else
            disp('Crack cannot reach surface');
            surface_failure = false;
            % find location where integral of stress is zero in upward
            % direction
            ind = find( stress_above > 0,1,'last');
            depth_tmp = depth(midpoint_ind:end);
            min_depth = depth_tmp(ind);
            % find depth at which crack stops
            ind = find(stress_below > 0,1,'last'); % net tension is ordered by increasing depth
            depth_tmp = depth(midpoint_ind:-1:1);
            max_depth = depth_tmp(ind);
        end
        fprintf('Relieving stresses between %e-%e m\n',min_depth,max_depth);
        result.failure_thickness(ifail) = max_depth-min_depth;
        result.failure_time(ifail) = time/seconds_in_year/1e6;
        result.failure_P(ifail) = Pex;
        result.failure_Pex_crit(ifail) = Pex_crit;
        
        result.failure_top(ifail) = min_depth;
        result.failure_bottom(ifail) = max_depth;
        result.failure_sigma_t{ifail} = sigma_t;
        result.failure_sigma_r{ifail} = sigma_r;
        result.failure_r{ifail} = grid_r;
        
        
        result.failure_z(ifail) = ztmp(end);
        
        ifail = ifail + 1;
        now_failing = depth >= min_depth & depth <= max_depth;
        failure_mask = failure_mask | now_failing;
        
        %                 failure_mask = false(size(sigma_r));
        %                 failure_mask(failure) = true;
        failure_time(now_failing) = time+dt;
    else
        no_longer_failing = failure_mask & (time - failure_time) >= 10*dtmin;
        if any(failure_mask(no_longer_failing))
            result.failure_dP(ifail-1) = Pex-result.failure_P(ifail-1);
        end
        if all(failure_mask) && any(failure_mask(no_longer_failing))
            %if erupted_volume > 0
            result.failure_erupted_volume(ifail-1) = erupted_volume;
            result.failure_erupted_volume_volumechange(ifail-1) = erupted_volume_volumechange;
            result.failure_erupted_volume_pressurechange(ifail-1) = erupted_volume_pressurechange;
            %end
            erupted_volume = 0;
            erupted_volume_volumechange = 0;
            erupted_volume_pressurechange = 0;
        end
        failure_mask(no_longer_failing) = false;
    end
    yielding = eiiD > (cohesion - 1/3*ei*friction); % note that compression is negative
    if any(yielding)
        keyboard
    end
    
    
    % 6. advance to next time step and plot (if needed)
    sigma_r_last = sigma_r;
    sigma_t_last = sigma_t;
    siiD_last = siiD;
    T_last = T;
    er_last = er;
    et_last = et;
    z_last = z;
    ur_last = ur;
    Pex_last = Pex;
    
    time = time + dt;
    
    if (time >= plot_times(iplot) || time >= t_end )
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
    if (time-last_store >= save_interval || time >= t_end || any(failure_mask))
        sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
        time_store(isave) = time;
        
        result.time(isave) = time;
        result.z(isave) = z;
        result.Ri(isave) = Ri;
        result.qb(isave) = Qbelow(time);
        result.sigma_t(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
        result.sigma_r(:,isave) = interp1(Ro-grid_r,sigma_r_last,save_depths);
        result.ur(:,isave) = interp1(Ro-grid_r,ur_last,save_depths);
        result.dTdr(:,isave) = interp1(Ro-grid_r,dTdotdr*dt,save_depths);
        result.T(:,isave) = interp1(Ro-grid_r,T,save_depths);
        result.Tb(isave) = Tb;
        result.Pex(isave) = Pex;
        result.XNH3(isave) = X;
        last_store = time; isave = isave+1;
    end
end
isave = isave-1;
result.parameters = p;
result.label = label;

result.time = result.time(1:isave);
result.z = result.z(1:isave);
result.Ri = result.Ri(1:isave);
result.qb = result.qb(1:isave);
result.sigma_t = result.sigma_t(:,1:isave);
result.sigma_r = result.sigma_r(:,1:isave);
result.ur = result.ur(:,1:isave);
result.dTdr = result.dTdr(:,1:isave);
result.T = result.T(:,1:isave);
result.Tb = result.Tb(1:isave);
result.Pex = result.Pex(1:isave);
result.XNH3(isave) = X;
result.save_depths = save_depths;

ifail = ifail-1;
result.failure_time = result.failure_time(1:ifail);
result.failure_P = result.failure_P(1:ifail);
result.failure_dP = result.failure_dP(1:ifail);

result.failure_Pex_crit = result.failure_Pex_crit(1:ifail);
result.failure_dP = result.failure_dP(1:ifail);
result.failure_thickness = result.failure_thickness(1:ifail);
result.failure_top = result.failure_top(1:ifail);
result.failure_bottom = result.failure_bottom(1:ifail);
result.failure_erupted_volume= result.failure_erupted_volume(1:ifail);
result.failure_erupted_volume_pressurechange = result.failure_erupted_volume_pressurechange(1:ifail);
result.failure_erupted_volume_volumechange = result.failure_erupted_volume_volumechange(1:ifail);
result.failure_sigma_t = result.failure_sigma_t(1:ifail);
result.failure_sigma_r = result.failure_sigma_r(1:ifail);

result.failure_r = result.failure_r(1:ifail);
result.failure_z = result.failure_z(1:ifail);

