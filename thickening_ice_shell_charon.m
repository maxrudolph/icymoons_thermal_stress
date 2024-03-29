% Script to solve coupled ice shell thermal and stress evolution
% Max Rudolph, March 19, 2020
clear;
close all;
addpath core; % this is where the helper functions live.

% Numerical parameters

nrs = [512];%[512];
failure_times = 0*nrs;
failure_thickness = 0*nrs;
for isetup = 3:3
    viscosity_model = 0; % 0 = Nimmo (2004), 1 = Goldsby and Kohlstedt (2001)
    viscosity.d = 1e-3; % grain size in m used to calculate the viscosity
    viscosity.P = 1e5; % Pressure in MPa used to calculate the viscosity
    
    if isetup == 1      % Europa
        Ro = 1.561e6;          % outer radius of ice shell (m)
        Ri = Ro-1e3;         % inner radius of ice shell (m)
        Rc = Ro-1.3e5;         % core radius (m)
        g = 1.3;        % used to calculate failure, m/s/s
        relaxation_parameter=1e-4; % used in nonlinear loop.

        label = 'Europa';
    elseif isetup == 2  % Enceladus
        Ro = 2.52e5;            % outer radius of ice shell (m)
        Ri = Ro-1.0e3;          % inner radius of ice shell (m)
        Rc = Ro-1.60e5;         % core radius (m)
        g = 0.113;        % used to calculate failure, m/s/s
        relaxation_parameter=1e-3; % used in nonlinear loop.

        label = 'Enceladus';
        
    elseif isetup == 3 % Charon
        Ro = 6.06e5;
        Rc = 3.76e5;
        Ri = Ro-1.0e5;
        g = 0.279;
        relaxation_parameter = 1e-3;
        label='Charon';
    else
        error('not implemented');
    end
    if viscosity_model == 0
        label = [label '_nimmovisc'];
    elseif viscosity_model == 1
        label = [label '_goldsbykohlstedt'];
    else
        error('not implemented');
    end
    
    for inr=1:length(nrs)
        ifail = 1; % index into list of times at which failure occurred.
        nr = nrs(inr); % number of grid points
        maxiter=1000;
        % Define physical constants and parameters
        % Physical constants
        seconds_in_year = 3.1558e7;
        R=8.314e-3;     % in kJ/mol/K
        % Boundary conditions and internal heating
        H=0; % internal heating rate.
        Tb=270;
        Ts=40;
        
        % Elastic and Viscous properties
        E = 5e9;        % shear modulus of ice (Pa)
        nu = 0.3;       % Poisson ratio of ice (-)
        beta_w = 4e-10; % Compressibility of water (1/Pa)
        alpha_l = 3e-5; % coefficient of linear thermal expansion ( alpha_v/3 ) (1/K)
        rho_i=900;      % density of ice (kg/m^3)
        rho_w=1000;     % density of water (kg/m^3)
        Q=40;           % activation energy, kJ/mol, Nimmo 2004 (kJ/mol)
        mub=1e15;       % basal viscosity (Pa-s)
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
%         k=kappa*rho_i*Cp;
        k = @(T) 651./T;
        %
        % Basal heating model - depends on thickness and transport properties
        %
%         Q0 = k*(Tb-Ts)/(Ro-Ri);% time-averaged basal heat flux
        [Q0,T,q] = find_steady_T(Ri,Ro,Tb,Ts,linspace(Ri,Ro,nr));
        perturbation_period = 1.0e8*seconds_in_year;
        deltaQonQ = 1.0; % fractional perturbation to Q0.
        
        Qbelow = @(time) 0;
        % calculate maxwell time at 100, 270
        fprintf('Maxwell time at surface, base %.2e %.2e\n',mu(100,0)/E,mu(Tb,0)/E);
        fprintf('Thermal diffusion timescale %.2e\n',(4e4)^2/kappa);
        % set end time and grid resolution
        t_end = 1e8*seconds_in_year;%  3*perturbation_period;
        % dt = 1e4*seconds_in_year; % time step in seconds
        dtmax = 1e4*seconds_in_year;
        dtmin = 3600;%*seconds_in_year;
        % dt1 = 3600; % size of first timestep
        % times = logspace(log10(dt1),log10(t_end+dt1),1e4)-dt1;
        plot_interval = t_end;
        save_interval = 1e2*seconds_in_year;        
        save_depths = linspace(0,Ro-Rc,500);
        
        nsave = ceil(t_end/save_interval) + 1;
        nsave_depths = length(save_depths);
        sigma_t_store = zeros(nsave_depths,nsave);
        
        results.time = zeros(nsave,1);
        results.z = zeros(nsave,1);
        results.Ri = zeros(nsave,1); results.Ri(1) = Ri;
        results.qb = zeros(nsave,1);
        results.sigma_t = NaN*zeros(nsave_depths,nsave);
        results.sigma_r = NaN*zeros(nsave_depths,nsave);
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
        results.failure_erupted_volume = zeros(1,nsave);
        results.failure_erupted_volume_pressurechange = zeros(1,nsave);
        results.failure_erupted_volume_volumechange = zeros(1,nsave);
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
        
        failure_mask = false(size(grid_r)); % stores whether failure occurred
        failure_time = zeros(size(grid_r)); % stores the time at which failure occurred
        
        while time < t_end
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
            
            qb_net = qb - Qbelow(time+dt);
            
            % thickening would be dx/dt = qb/(L*rho_i)
            delta_rb = dt*qb_net/Lf/rho_i;
            z = z_last + delta_rb;
            dzdt = delta_rb/dt;
            
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
                if iter>10
                    [tmp,ind] = unique(pex_store(1:iter-1));
                    Pex = interp1(pexpost_store(ind)-pex_store(ind),pex_store(ind),0,'linear','extrap');
                elseif iter>1
                    Pex = Pex + relaxation_parameter*(Pex_post-Pex);
                else
                    Pex = Pex_last;
                end
                Pex_crit = (rho_w-rho_i)*(Ro-(Ri-z))*g;

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
                if abs( Pex_post-Pex )/abs(Pex) < 1e-3 || abs(Pex_post-Pex) < 1e2
                    fprintf('dt=%.2e yr, time=%.3e Myr, Pex_post %.6e Pex %.6e, converged in %d iterations\n',dt/seconds_in_year,(time+dt)/seconds_in_year/1e6,Pex_post,Pex,iter);
                    converged = true;
                elseif iter==maxiter
                    error('Nonlinear loop failed to converge');
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
                results.failure_thickness(ifail) = max_depth-min_depth;
                results.failure_time(ifail) = time/seconds_in_year/1e6;
                results.failure_P(ifail) = Pex;
                results.failure_Pex_crit(ifail) = Pex_crit;

                results.failure_top(ifail) = min_depth;
                results.failure_bottom(ifail) = max_depth;
                results.failure_sigma_t{ifail} = sigma_t;
                results.failure_sigma_r{ifail} = sigma_r;
                results.failure_r{ifail} = grid_r;
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
                    if erupted_volume > 0
                        results.failure_erupted_volume(ifail-1) = erupted_volume;
                        results.failure_erupted_volume_volumechange(ifail-1) = erupted_volume_volumechange;
                        results.failure_erupted_volume_pressurechange(ifail-1) = erupted_volume_pressurechange;
                    end
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
                
                results.time(isave) = time;
                results.z(isave) = z;
                results.Ri(isave) = Ri;
                results.qb(isave) = Qbelow(time);
                results.sigma_t(:,isave) = interp1(Ro-grid_r,sigma_t_last,save_depths);
                results.sigma_r(:,isave) = interp1(Ro-grid_r,sigma_r_last,save_depths);
                results.ur(:,isave) = interp1(Ro-grid_r,ur_last,save_depths);
                results.dTdr(:,isave) = interp1(Ro-grid_r,dTdotdr*dt,save_depths);
                results.T(:,isave) = interp1(Ro-grid_r,T,save_depths);
                results.Pex(isave) = Pex;
                results.Pex_crit(isave) = Pex_crit;
                last_store = time; isave = isave+1;
            end
        end
        %%
        mask = 1:(isave-1);
        figure();
        
        subplot(3,1,1);
        plot(results.time(mask)/seconds_in_year,results.z(mask))
        ylabel('Amount of freezing (m)');
        subplot(3,1,2);
        plot(results.time(mask)/seconds_in_year,results.sigma_t(1,mask),'DisplayName',sprintf('%.02f km',save_depths(1)/1000));
        hold on
        plot(results.time(mask)/seconds_in_year,results.sigma_t(10,mask),'DisplayName',sprintf('%.02f km',save_depths(10)/1000));
        plot(results.time(mask)/seconds_in_year,results.sigma_t(20,mask),'DisplayName',sprintf('%.02f km',save_depths(20)/1000));
        legend();
        ylabel('\sigma_t (Pa)')
        subplot(3,1,3);
        plot(results.time(mask)/seconds_in_year,results.Pex(mask));
        ylabel('Overpressure (Pa)');
        
        %% Pseudocolor stress plot
        figure();
        xscale = 'log';
        t=tiledlayout(3,1,'TileSpacing','none','Padding','none');
        t.Units = 'centimeters';
        t.OuterPosition = [1 1 9.5 9];
        nexttile
        
        contourf(results.time(mask)/seconds_in_year,save_depths/1000,results.sigma_t(:,mask),64,'Color','none'); %shading flat;
        hold on
        plot(results.time(mask)/seconds_in_year,((Ro-results.Ri(mask))+results.z(mask))/1000,'Color','k','LineWidth',1);
        set(gca,'YLim',[0 ceil(1+max(((Ro-results.Ri(mask))+results.z(mask))/1000))]);
        set(gca,'YDir','reverse');
        ax1 = gca();
        ax1.FontSize=8;
        hcb = colorbar();
        hcb.Label.String = 'Tensile Stress (Pa)';
        text(0.025,0.85,char('A'+(isetup-1)*2),'FontSize',8,'Units','normalized');
        xlabel('Time (years)');
        title(label);
        ylabel('Depth (km)');
        set(gca,'XScale',xscale);
        hold on;
        for i=1:ifail-1
            plot(results.failure_time(i)*1e6*[1 1],[results.failure_top(i) results.failure_bottom(i)]/1e3,'r');
        end
        nexttile
        plot(results.time(mask)/seconds_in_year,results.Pex(mask));
        plot(results.time(mask)/seconds_in_year,results.Pex_crit(mask));
        ylabel('P_{ex} (Pa)');
        set(gca,'XScale',xscale);
        ax2 = gca();
        ax2.Position(3) = ax1.Position(3);
        ax2.XLim = ax1.XLim;
        ax2.FontSize=8;
        hold on
        plot(results.failure_time(1:ifail-1)*1e6,results.failure_P(1:ifail-1),'r.');
        end_color = [0 0.9 0];
        plot(results.failure_time(1:ifail-1)*1e6,(results.failure_P(1:ifail-1)+results.failure_dP(1:ifail-1)),'LineStyle','none','Color',end_color,'Marker','o','MarkerFaceColor',end_color,'MarkerSize',2);
        text(0.025,0.85,char('B'+(isetup-1)*2),'FontSize',8,'Units','normalized');

        xlabel('Time (years)');
        nexttile
        hold on;
        for i=1:ifail-1
            plot(results.failure_time(i)*1e6*[1 1],results.failure_erupted_volume(i)/(4*pi*Ro^2)*[0 1],'b');
            %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_volumechange(i)/(4*pi*Ro^2),'go');
            %         plot(results.failure_time(i)*1e6,results.failure_erupted_volume_pressurechange(i)/(4*pi*Ro^2),'rx');
        end
        ylabel('Erupted vol. (m)');
        xlabel('Time (years)');
        set(gca,'XScale',xscale);
        ax3=gca();
        ax3.XLim = ax1.XLim;
        ax3.Position(3) = ax1.Position(3);
        ax3.Box = 'on';
        ax3.FontSize=8;
        text(0.025,0.85,char('C'+(isetup-1)*3),'FontSize',8,'Units','normalized');

        fig = gcf();
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition(3) = 6.00;
        fig.Color = 'w';
        %exportgraphics(gcf,[label '_thickening.eps'],'ContentType','vector');
        
        %% Make a figure showing tensile stress at the time of failure
        figure();
        hold on
        for i=1:ifail-1
            plot(results.failure_sigma_t{i}/1e6,(results.failure_r{i}(end)-results.failure_r{i})/1000);
        end
        set(gca,'YDir','reverse');
        xlabel('Tensile Stress (MPa)');
        ylabel('Depth (km)');%
        %% write failure stresses to a file
        odir = [label '_failure_stress'];
        mkdir(odir);
        for i=1:ifail-1
            filename = [odir '/' sprintf('failure_stresses_%d_%.2fkm.csv',i,results.failure_thickness(i)/1000)];
            fh = fopen(filename,'w');
            fprintf(fh,'# Failure event %d, crack thickness %d m , ice shell thickness %d m\n',i,round(results.failure_thickness(i)),round(results.failure_r{i}(end)-results.failure_r{i}(1)));
            fprintf(fh,'# Radius (m), sigma_t (Pa), sigma_r (Pa)\n');
            for j=1:nr
                fprintf(fh,'%e,%e,%e\n',results.failure_r{i}(j),results.failure_sigma_t{i}(j),results.failure_sigma_r{i}(j));
            end
            fclose(fh);
        end
    end           
end
    
