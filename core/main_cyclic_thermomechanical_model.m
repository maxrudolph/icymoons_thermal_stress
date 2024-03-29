function results = main_cyclic_thermomechanical_model(parameters)
% parameters should contain the following
Ro = parameters.Ro;     % outer radius of ice shell (m)
Ri = parameters.Ri;     % inner radius of ice shell (m)
Rc = parameters.Rc;     % core radius (m)
g = parameters.g;       % gravitational acceleration, m/s/s
k = parameters.k;       % Thermal conductivity function
tensile_strength = parameters.tensile_strength;     % tensile strength, Pa
perturbation_period = parameters.perturbation_period; % period of perturbation (s)
deltaQonQ = parameters.deltaQonQ;   % fractional perturbation to Q0.
relaxation_parameter=parameters.relaxation_parameter;   % used in nonlinear loop.
viscosity_model = parameters.viscosity_model;
Tb=parameters.Tb;       % Basal temperature (K)
Ts=parameters.Ts;       % Surface temperature (K)
t_end = parameters.end_time; % Simulation end time (s)
save_interval = parameters.save_interval; % Time between saved results (s)
save_start = parameters.save_start; % Time to begin saving (s)

viscosity.d = 1e-3; % grain size in m used to calculate the viscosity
viscosity.P = 1e5; % Pressure in MPa used to calculate the viscosity

% Numerical parameters
ifail = 1;              % index into list of times at which failure occurred.
nr = parameters.nr;               % Number of nodes in radial direction
maxiter=1000;           % Maximum iterations in nonlinear loop
% set up the grid
grid_r = linspace(Ri,Ro,nr); % set up the grid

% Basic physical properties
rho_i=900;              % density of ice (kg/m^3)
rho_w=1000;             % density of water (kg/m^3)

% Thermal properties
Cp = 2100;              % heat capacity of ice (J/kg/K)
Lf = 334*1000;          % latent heat of fusion (J/kg)
kappa = 1e-6;           % Thermal diffusivity (m/s/s)
% k=kappa*rho_i*Cp;       % Thermal conductivity (W/m/K)
R=8.314e-3;             % Ideal Gas constant (kJ/mol/K)

% Heating model
% Q0 = k*(parameters.Tb-parameters.Ts)/(parameters.Ro-parameters.Ri);
% Calculate the steady state basal heat flux, including k(T) = k1/T
% temperature dependence.
[~,T_last,qr] = find_steady_T(parameters.Ri,parameters.Ro,parameters.Tb,parameters.Ts,grid_r);
Q0 = 4*pi*parameters.Ro^2 * qr(end);
% Qtot is a function that returns the heating rate in W:
Qtot = @(time) Q0*(1+deltaQonQ*sin(-2*pi*time/perturbation_period));

% Physical constants
seconds_in_year = 3.1558e7;
% Boundary conditions and internal heating
H=0; % internal heating rate.

% Elastic and Viscous properties
E = 5e9;        % shear modulus of ice (Pa)
nu = 0.3;       % Poisson ratio of ice (-)
beta_w = 4e-10; % Compressibility of water (1/Pa)
alpha_l = 3e-5; % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)

Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
mub=1e15;       % basal viscosity (Pa-s)
mu = @(T) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
if viscosity_model == 0
    mu = @(T,stress) mub*exp(Q*(Tb-T)/R/Tb./T); % function to evaluate viscosity in Pa-s given T
elseif viscosity_model == 1
    mu = @(T,stress) goldsby_kohlstedt(stress,T,viscosity.d,viscosity.P); % Goldsby-Kohlstedt effective viscosity
end


% Failure criterion:
cohesion = 2e7;  % plastic yield strength
friction = 0.6; % friction angle for plastic yielding

% calculate maxwell time at 100, 270
fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(Ts,0)/E,mu(Tb,0)/E);
fprintf('Thermal diffusion timescale %.2e\n',(Ro-Ri)^2/kappa);

% set end time and grid resolution
dtmax = 5e4*seconds_in_year;
dtmin = 1e1*seconds_in_year;

save_depths = linspace(0,(parameters.Ro-parameters.Ri)/(1-parameters.deltaQonQ)*1.25,200);
nsave = ceil((t_end-save_start)/save_interval) + 1;
nsave_depths = length(save_depths);
sigma_t_store = zeros(nsave_depths,nsave);

results.time = zeros(nsave,1);
results.save_depths = save_depths;
results.z = zeros(nsave,1);
results.Ri = zeros(nsave,1); results.Ri(1) = Ri;
results.qb = zeros(nsave,1);
results.Qtot = zeros(nsave,1);
results.sigma_t = NaN*zeros(nsave_depths,nsave);
results.sigma_r = NaN*zeros(nsave_depths,nsave);
results.siiD = NaN*zeros(nsave_depths,nsave);
results.sigma_yield = NaN*zeros(nsave_depths,nsave);
results.Pex = zeros(nsave,1);
results.Pex_crit = zeros(nsave,1);
results.dTdr = zeros(nsave_depths,nsave);
results.T = zeros(nsave_depths,nsave);
results.ur = zeros(nsave_depths,nsave);
results.failure_time = zeros(1,nsave);
results.failure_P = zeros(1,nsave);
results.failure_Pex_crit = zeros(1,nsave);
results.failure_dP = zeros(1,nsave);
results.failure_thickness = zeros(1,nsave);
results.failure_top = zeros(1,nsave);
results.failure_bottom = zeros(1,nsave);
results.failure_initial = zeros(1,nsave);
results.failure_erupted_volume = zeros(1,nsave);
results.failure_eruption_time = zeros(1,nsave);
results.failure_erupted_volume_pressurechange = zeros(1,nsave);
results.failure_erupted_volume_volumechange = zeros(1,nsave);
erupted_volume = 0;
erupted_volume_pressurechange = 0;
erupted_volume_volumechange = 0;



% initialize solution vectors (IC)
sigma_r_last = zeros(nr,1); % initial stresses
sigma_t_last = zeros(nr,1); % initial stresses
siiD_last = zeros(nr,1); % deviatoric stress invariant - used for viscosity


% Initialize T with steady numerical solution.
for iter=1:10
    % Iterate 10 times to ensure consistency between assumed and actual k(T)
    % In practice, this doesn't matter if we start with the correct
    % analytic solution already.
    T_last = solve_temperature_shell(grid_r,T_last,Tb,Ts,k,rho_i,Cp,H,Inf,0.0);
end
er_last = zeros(nr,1); % strains
et_last = zeros(nr,1);
ur_last = zeros(nr,1); % displacement
z_last = 0;    % total amount of thickening
Pex_last = 0; %initial overpressure
siiD_last = zeros(nr,1);
grid_r_last = grid_r;

time=0; % elapsed time
dt = dtmax;

% save initial state
isave = 1;
last_store = -Inf;
if( save_start == 0 )
    sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
    time_store(isave) = time;
    last_store = time; isave = isave+1;
end
failure_mask = false(size(grid_r)); % stores whether failure occurred
failure_time = zeros(size(grid_r)); % stores the time at which failure occurred
terminate = 0;
while time < t_end && ~terminate
    % In each timestep, we do the following
    % 1. Calculate the amount of basal freeze-on and advance the mesh
    % 2. Solve the heat equation using an implicit method
    % 3. Solve for sigma_r
    % 4. Calculate sigma_t
    % 5. Calculate the radial displacements u(r)
    
    % 1. Calculate basal freeze-on and interpolate old solution onto new mesh
    % calculate heat flux
    dtlast = dt;
    dt = min(dtmax,2*dtlast); % limit amount by which timestep can increase.
    Tg = Tb-(T_last(2)-Tb);
    dTdr_b_last = (T_last(2)-Tg)/2/(grid_r(2)-grid_r(1)); % calculate dT/dr at ocean-ice interface using centered difference
    qb = -k(parameters.Tb)*dTdr_b_last;
    qbelow = Qtot(time+dt) / ( 4*pi*grid_r(1)^2 ); % This is q=Q/area
    qb_net = qb - qbelow; % first term is conducted heat. second term is heat supplied from below.
    
    % determine the timestep
    if abs(qb_net/Lf/rho_i*dt) > (grid_r(2)-grid_r(1))/2
        dt = abs( (grid_r(2)-grid_r(1))/2/(qb_net/Lf/rho_i) );
    end
    if dt > perturbation_period/200
        dt = perturbation_period/200;
    end
    if dt < dtmin
        dt = dtmin;
        disp('Setting dt = dtmin');
    end
    if any(failure_mask)
        if Pex < 0
            dt =dtmax/10;
        else
        dt = dtmin;
        end
    end
    
    qbelow = Qtot(time+dt) / ( 4*pi*grid_r(1)^2 ); % This is q=Q/area
    qb_net = qb - qbelow;
    
    % thickening would be dx/dt = qb/(L*rho_i)
    delta_rb = dt*qb_net/Lf/rho_i;
    z = z_last + delta_rb;
    dzdt = delta_rb/dt;
    
    % calculate new ocean pressure (Manga and Wang 2007, equation 5)
    Pex_pred = Pex_last + 3*(Ri-z)^2/beta_w/((Ri-z)^3-Rc^3)*(delta_rb*(rho_w-rho_i)/rho_w-ur_last(1)); % ur_last because we don't yet know the uplift
    Pex_crit = (rho_w-rho_i)*(Ro-(Ri-z))*g;
    
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
        if iter == 1
            Pex = Pex_last;
            extra_relaxation = 1.0;
        elseif iter>10 && iter < 30
            [tmp,ind] = sort(pexpost_store(1:iter-1)-pex_store(1:iter-1));
            [tmp1,ind1] = unique(tmp);
            Pex = interp1(tmp1,pex_store(ind(ind1)),0,'linear','extrap');
        elseif iter > 1
            Pex = Pex + extra_relaxation*relaxation_parameter*(Pex_post-Pex);
        end
        if ~mod(iter,100)
            Pex = Pex_last;
            extra_relaxation = extra_relaxation / 2;
        end
        if iter == maxiter
            warning('last iteration, no convergence');
        end
        
        % calculate viscosity at each node
        for ivisc=1:10
            mu_node = zeros(nr,1);
            mu_node(:) = mu(T,siiD_last);
            % reduce Maxwell time in cracked region
            if all(failure_mask) % crack reaches ocean - relieve all overpressure
                if Pex_last >= Pex_crit % allow for extrusion if P>Pex_crit
                    % Calculate the volume erupted (dP)*beta*V0 + V-V0
                    pressure_contribution = (Pex_last-Pex_crit)*beta_w*(4/3*pi*((Ri-z)^3-Rc^3));
                    urelax = (Ri-z)/E*(1-2*nu)*(Pex_last-Pex_crit); % Manga and Wang (2007) equation 4
                    volume_contribution = (Ri-z)^2*urelax*4*pi; % (4*pi*R^2)*dr
                else
                    pressure_contribution = 0;
                    volume_contribution = 0;
                end
                % reset stresses and uplift
                %             sigma_r_last = 0*sigma_r_last;
                %             sigma_t_last = 0*sigma_t_last;
                %             er_last = 0*er_last;
                %             et_last = 0*et_last;
                %             ur_last = 0*ur_last;
                minimum_viscosity_prefactor = 0; % maximum allowable fractional reduction in viscosity
                mu_node(failure_mask) = min(mu_node(failure_mask),max(minimum_viscosity_prefactor*mu_node(failure_mask),0.1*E*dt));  % timest
                Pex=0; % force zero pressure under the assumption that eruptions relieve overpressure
                converged = true;
                Ri = Ri - z;
                z = 0;
                z_last=0;
                % move the inner radius effectively to the current position
                % of the base of the ice shell. Then set the amount of
                % freezing to zero.
                % It's ok to increase the erupted volume counter inside the
                % nonlinear loop because Pex is fixed and no iteration will
                % occur.
                erupted_volume = erupted_volume + pressure_contribution + volume_contribution;
                erupted_volume_pressurechange = erupted_volume_pressurechange + pressure_contribution;
                erupted_volume_volumechange = erupted_volume_volumechange + volume_contribution;
            elseif( any(failure_mask) ) % Crack does not reach ocean - overpressure is not relieved.
                minimum_viscosity_prefactor = 0; % maximum allowable fractional reduction in viscosity
                mu_node(failure_mask) = min(mu_node(failure_mask),max(minimum_viscosity_prefactor*mu_node(failure_mask),0.1*E*dt));  % timestep = 10x(maxwell time)
                if iter==1
                    Pex=0; % If failure occurs, it's better to guess that all pressure is relieved. Other choices can cause convergence problems.
                end
            else
                %5a. Calculate the yield stress (sigma_yield) and check for
                %plastic yielding.
                si = sigma_r_last + 2*sigma_t_last; % trace of stress tensor.
                lithostatic_pressure = rho_i*g*(grid_r_last(end)-grid_r_last');
                sigma_yield = max((cohesion + friction*(1/3*si + lithostatic_pressure)),0);
                yielding = siiD_last > sigma_yield; % note that compression is negative. positive I_sigma means compression
                if any(yielding) % reduce viscosity to yielding viscosity.
                    if ivisc>1
                        siiD_last = siiD_post;
                    end
                    siiD_el = (E*dt + mu_node)./mu_node .* siiD_last; % gerya equation 13.50 - elastic stresses
                    mu_vp = E*dt*sigma_yield./(siiD_el - sigma_yield); % effective viscoplastic viscosity
                    min_visc = mu_node/1e3; % limit yielding for numerical stability
                    mu_node(yielding) = min(mu_node(yielding),mu_vp(yielding));
                    mu_node = max(mu_node,min_visc);
                end
            end
            
            % Calculate stresses.
            [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-Pex,E,nu,dt);
            siiD_post = sqrt( 0.5*(sigma_rD.^2 + 2*sigma_tD.^2) );
            if ~any(yielding)
                break;
            end
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
        
        ei = 2*de_t + de_r; % first invariant
        de_tD = de_t - 1/3*ei;
        de_rD = de_r - 1/3*ei;
        eiiD = sqrt( 0.5*(de_rD.^2 + 2*de_tD.^2) ); % second invariant of deviatoric strain
        siiD = siiD_post;
        
        
        % re-calculate excess pressure using new uplift
        %             Pex_post = 3*Ri^2/beta_w/(Ri^3-Rc^3)*(z*(rho_w-rho_i)/rho_w-ur(1));
        Pex_post = Pex_last + 3*(Ri-z)^2/beta_w/((Ri-z)^3-Rc^3)*((z-z_last)*(rho_w-rho_i)/rho_w-(ur(1)-ur_last(1)));
        %         fprintf('iter %d. Pex_post %.2e Pex %.2e\n',iter,Pex_post,Pex);
        
        % calculate an approximate derivative d_sigma_r/dPex
        %         perturb = max(1e-4,abs(1e-6*Pex));
        %         [sigma_rp,sigma_tp,~,sigma_tDp] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-(Pex+perturb),E,nu,dt);
        %         [sigma_rm,sigma_tm,~,sigma_tDm] = solve_stress_viscoelastic_shell(grid_r,mu_node,sigma_r_last,alpha_l*dTdotdr,-(Pex-perturb),E,nu,dt);
        %         dsigma_tp = sigma_tp - sigma_t_last; dsigma_rp = sigma_rp-sigma_r_last;
        %         dsigma_tm = sigma_tm - sigma_t_last; dsigma_rm = sigma_rm-sigma_r_last;
        %
        %         de_tp = 1/E*(dsigma_tp-nu*(dsigma_tp+dsigma_rp))+alpha_l*dT + dt/2*(sigma_tDp./mu_node); % change in tangential strain
        %         de_tm = 1/E*(dsigma_tm-nu*(dsigma_tm+dsigma_rm))+alpha_l*dT + dt/2*(sigma_tDm./mu_node); % change in tangential strain
        %         du_p = grid_r(1) * de_tp(1);
        %         du_m = grid_r(1) * de_tm(1);
        %         du_dp = (du_p-du_m)/(2*perturb); % d(uplift)/d(pressure)
        
        % check for convergence
        if abs( Pex_post-Pex )/abs(Pex) < 1e-2
            fprintf('dt=%.2e yr, time=%.3e Myr, Pex_post %.2e Pex %.2e, converged in %d iterations\n',dt/seconds_in_year,(time+dt)/seconds_in_year/1e6,Pex_post,Pex,iter);
            converged = true;
            
        end
        
        pex_store(iter) = Pex;
        pexpost_store(iter) = Pex_post;
        if converged
            break;
        end
    end%end nonlinear loop for pressure
    
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
    % require that the deviatoric stress in the tangential direction be consistent with vertical cracks, as required.
    % the tensile strength / 10 is meant to enforce a small positive
    % deviatoric stress
    failure = failure & sigma_tD > tensile_strength/10;
    if( any(failure) ) % failure is occurring
        disp(['Failure criterion has been reached']);
        idx_shallow = find(failure,1,'last');
        idx_deep = find(failure,1,'first');
        %         fprintf('Shallowest, deepest failure: %f, %f\n\n',Ro-grid_r(idx_shallow),Ro-grid_r(idx_deep));
        %         fprintf('Failure time: %f Myr\n',time / seconds_in_year / 1e6);
        %         fprintf('Surface stress at failure: %f MPa\n',sigma_t(end)/1e6);
        
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
        if midpoint_ind ~= 1
            stress_below = cumtrapz( depth(midpoint_ind:-1:1), sigma_t_tot(midpoint_ind:-1:1) ); % integrate in downward direction
        else
           stress_below = 0; 
        end
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
            ind = find( stress_above >= 0,1,'last');
            depth_tmp = depth(midpoint_ind:end);
            min_depth = depth_tmp(ind);
            % find depth at which crack stops
            ind = find(stress_below >= 0,1,'last'); % net tension is ordered by increasing depth
            depth_tmp = depth(midpoint_ind:-1:1);
            max_depth = depth_tmp(ind);
        end
        fprintf('Relieving stresses between %e-%e m\n',min_depth,max_depth);
        results.failure_thickness(ifail) = max_depth-min_depth;
        results.failure_time(ifail) = time;
        results.failure_P(ifail) = Pex;
        results.failure_Pex_crit(ifail) = Pex_crit;
        results.failure_top(ifail) = min_depth;
        results.failure_bottom(ifail) = max_depth;
        results.failure_initial(ifail) = midpoint_depth;
        ifail = ifail + 1;
        now_failing = depth >= min_depth & depth <= max_depth;
        failure_mask = failure_mask | now_failing;
        
        %                 failure_mask = false(size(sigma_r));
        %                 failure_mask(failure) = true;
        failure_time(now_failing) = time+dt;
    else
        no_longer_failing = failure_mask & (time - failure_time) >= 10*dtmin;
        if any(failure_mask(no_longer_failing))
            results.failure_dP(ifail-1) = Pex-results.failure_P(ifail-1);
        end
        if all(failure_mask) && any(failure_mask(no_longer_failing))
            % archive results of this eruption.
            if erupted_volume > 0
                results.failure_eruption_time(ifail-1) = time;
                results.failure_erupted_volume(ifail-1) = erupted_volume;
                results.failure_erupted_volume_volumechange(ifail-1) = erupted_volume_volumechange;
                results.failure_erupted_volume_pressurechange(ifail-1) = erupted_volume_pressurechange;
            else
                results.failure_eruption_time(ifail-1) = time;
                results.failure_erupted_volume(ifail-1) = -1;
                results.failure_erupted_volume_volumechange(ifail-1) = -1;
                results.failure_erupted_volume_pressurechange(ifail-1) = -1;
            end
            % reset erupted volume.
            erupted_volume = 0;
            erupted_volume_volumechange = 0;
            erupted_volume_pressurechange = 0;
        end
        failure_mask(no_longer_failing) = false;
    end
    
    % 6. advance to next time step and plot (if needed)
    sigma_r_last = sigma_r;
    sigma_t_last = sigma_t;
    grid_r_last = grid_r;
    siiD_last = siiD;
    T_last = T;
    er_last = er;
    et_last = et;
    z_last = z;
    ur_last = ur;
    Pex_last = Pex;
    
    time = time + dt;
    
    if time >= save_start && (time-last_store >= save_interval || time >= t_end || any(failure_mask) || any(no_longer_failing))
        sigma_t_store(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
        time_store(isave) = time;
        
        results.time(isave) = time;
        results.z(isave) = z;
        results.Ri(isave) = Ri;
        results.Qtot(isave) = Qtot(time);
        results.qb(isave) = qbelow;
        results.sigma_t(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
        results.sigma_r(:,isave) = interp1(Ro-grid_r,sigma_r_last,save_depths);
        results.siiD(:,isave) = interp1(Ro-grid_r,siiD,save_depths);
        results.sigma_yield(:,isave) = interp1(Ro-grid_r,sigma_yield,save_depths);
        
        results.ur(:,isave) = interp1(Ro-grid_r,ur_last,save_depths);
        results.dTdr(:,isave) = interp1(Ro-grid_r,dTdotdr*dt,save_depths);
        results.T(:,isave) = interp1(Ro-grid_r,T,save_depths);
        results.Pex(isave) = Pex;
        results.Pex_crit(isave) = Pex_crit;
        last_store = time; isave = isave+1;
    end
    % check for an unphysically low ocean pressure
    Ptot = Pex + rho_i*g*(grid_r(end)-grid_r(1));
    if Ptot < 612
        disp('Ocean pressure below 6 mbar. Ending calculation.');
        terminate=1;
        break;
    end
    if iter==maxiter
        for i=1:maxiter
            fprintf('iter %d. Pex_post %.12e Pex %.12e\n',i,pexpost_store(i),pex_store(i));
        end
        terminate=2;
        break;
    end
end
%% Trim results structure
ifail2 = find(results.failure_time>0,1,'last');
ifail1 = find(results.failure_time>=save_start,1,'first');
fail_mask = false(size(results.failure_time));
fail_mask(ifail1:ifail2) = true;

isave = find(results.time>0,1,'last');
% general information:
results.time = results.time(1:isave);
results.z = results.z(1:isave);
results.Ri = results.Ri(1:isave,:);
results.qb = results.qb(1:isave);
results.Qtot = results.Qtot(1:isave);
results.sigma_t = results.sigma_t(:,1:isave);
results.sigma_r = results.sigma_r(:,1:isave);
results.siiD = results.siiD(:,1:isave);
results.sigma_yield = results.sigma_yield(:,1:isave);

results.Pex = results.Pex(1:isave);
results.Pex_crit = results.Pex_crit(1:isave);
results.dTdr = results.dTdr(:,1:isave);
results.T = results.T(:,1:isave);
results.ur = results.ur(:,1:isave);
% failure information:
results.failure_time = results.failure_time(fail_mask);
results.failure_P = results.failure_P(fail_mask);
results.failure_Pex_crit = results.failure_Pex_crit(fail_mask);
% results.failure_dP = results.failure_dP(fail_mask);
results.failure_thickness = results.failure_thickness(fail_mask);
results.failure_top = results.failure_top(fail_mask);
results.failure_bottom = results.failure_bottom(fail_mask);
results.failure_initial = results.failure_initial(fail_mask);
if Ptot < 612    
    results.termination = 1;
elseif iter == maxiter
    results.termination = 2;
else
    results.termination = 0;
end
ifail2 = find(results.failure_eruption_time>0,1,'last');
ifail1 = find(results.failure_eruption_time>=save_start,1,'first');
fail_mask = false(size(results.failure_eruption_time));
fail_mask(ifail1:ifail2) = true;
results.failure_eruption_time = results.failure_eruption_time(fail_mask);
results.failure_erupted_volume = results.failure_erupted_volume(fail_mask);
results.failure_erupted_volume_pressurechange = results.failure_erupted_volume_pressurechange(fail_mask);
results.failure_erupted_volume_volumechange = results.failure_erupted_volume_volumechange(fail_mask);
